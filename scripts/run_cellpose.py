#!/usr/bin/env python3
"""Phase 2: Cellpose 3D inference on BA-TOS crops, one fixed YAML config per model (no per-crop tuning).

Runs in the cellpose env (C:\\BA_TOS_DATA\\envs\\cellpose, Python 3.12, torch cu128):

    C:\\BA_TOS_DATA\\envs\\cellpose\\Scripts\\python.exe -I scripts\\run_cellpose.py \\
        --config configs\\cellpose_cpsam_v2.yaml --set dev|test|all|<crop_id> [<crop_id> ...]

Reads <crop>/gray.tif (zyx), runs CellposeModel.eval with the config's `eval` arguments (everything else at the
cellpose defaults), writes <runs_root>/<tag>/<crop>/labels.tif and run.json (model, versions, device, weights
sha256, resolved eval arguments, runtime, peak VRAM, peak host RSS, git commit, command). Prediction only:
evaluation is done separately with batos_v2.evaluate_v2 in the batos env. A crop whose run.json exists is
skipped (resume); a crash writes FAILED.json and the driver continues. GT is never read here.

`--override key=value` (repeatable) changes an eval argument for a labelled sensitivity run; the run tag
then defaults to <model>__<key>-<value>. Values are parsed as YAML. The special value `diameter=2r_granulo`
sets diameter = 2 * r_granulo per crop, read from the frozen v3 run record of that crop (GT-free scale estimate;
dev: v3_scale_test/runs/v3/x1, test: test_protocol/v3); a missing record is an error.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import platform
import subprocess
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import psutil
import tifffile
import yaml

REPO = Path(__file__).resolve().parents[1]
MANIFESTS = {"dev": REPO / "data" / "metadata" / "crops_manifest.csv",
             "test": REPO / "data" / "metadata" / "crops_manifest_test.csv"}
V3_RUNS = {"dev": Path(r"C:\BA_TOS_DATA\runs\v3_scale_test\runs\v3\x1"),
           "test": Path(r"C:\BA_TOS_DATA\runs\test_protocol\v3")}


def r_granulo(crop: str, split: str) -> float:
    j = V3_RUNS[split] / crop / "v3_run.json"
    if not j.exists():
        raise SystemExit(f"no v3 run record for {crop}: {j}")
    r = json.loads(j.read_text(encoding="utf-8")).get("r_granulo")
    if r is None:
        raise SystemExit(f"r_granulo missing in {j}")
    return float(r)


def crops_for(selection: list[str], extra_manifest: Path | None = None) -> list[tuple[str, Path, str]]:
    rows = {}
    manifests = dict(MANIFESTS)
    if extra_manifest is not None:  # e.g. degraded inputs for the robustness test (split name = "extra")
        manifests["extra"] = extra_manifest
    for key, path in manifests.items():
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                rows[r["crop_id"]] = (key, Path(r["gray_path"]))
    out = []
    for s in selection:
        if s in ("dev", "test", "extra"):
            out += [(c, p, k) for c, (k, p) in rows.items() if k == s]
        elif s == "all":
            out += [(c, p, k) for c, (k, p) in rows.items()]
        elif s in rows:
            out.append((s, rows[s][1], rows[s][0]))
        else:
            raise SystemExit(f"unknown crop or set: {s}")
    seen, uniq = set(), []
    for c, p, k in out:
        if c not in seen:
            seen.add(c)
            uniq.append((c, p, k))
    return uniq


def git_info() -> dict:
    def g(*a):
        return subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True).stdout.strip()
    return {"commit": g("rev-parse", "HEAD"), "dirty": bool(g("status", "--porcelain", "--untracked-files=no"))}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


class RssPeak:
    def __init__(self):
        import threading
        self.proc, self.peak, self.stop = psutil.Process(), 0, threading.Event()
        self.t = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self.stop.is_set():
            self.peak = max(self.peak, self.proc.memory_info().rss)
            time.sleep(0.2)

    def __enter__(self):
        self.t.start()
        return self

    def __exit__(self, *a):
        self.stop.set()
        self.t.join()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--set", nargs="+", required=True, dest="selection", help="dev | test | all | crop ids")
    ap.add_argument("--override", action="append", default=[], help="eval key=value (labelled sensitivity run)")
    ap.add_argument("--tag", help="output sub-folder (default: model, or model__overrides)")
    ap.add_argument("--manifest", type=Path, help="extra CSV (crop_id, gray_path); its crops form the set 'extra'")
    args = ap.parse_args(argv)

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    eval_kw = dict(cfg["eval"])
    overrides = {}
    for o in args.override:
        k, v = o.split("=", 1)
        overrides[k] = yaml.safe_load(v)
    eval_kw.update(overrides)
    tag = args.tag or cfg["model"] + "".join(f"__{k}-{v}" for k, v in overrides.items())
    root = Path(cfg["runs_root"]) / tag
    crops = crops_for(args.selection, args.manifest)

    import torch
    import cellpose
    from cellpose import models

    if cfg["model"] not in models.MODEL_NAMES:
        raise SystemExit(f"model {cfg['model']!r} not in cellpose MODEL_NAMES {models.MODEL_NAMES}")
    if cfg.get("gpu") and not torch.cuda.is_available():
        raise SystemExit("config requests GPU but torch.cuda.is_available() is False")
    dev = torch.device("cuda" if cfg.get("gpu") else "cpu")
    model = models.CellposeModel(gpu=bool(cfg.get("gpu")), pretrained_model=cfg["model"])
    weights = Path(models.cache_model_path(cfg["model"]))
    sig = inspect.signature(model.eval)
    resolved = {k: (eval_kw[k] if k in eval_kw else p.default) for k, p in sig.parameters.items()
                if k not in ("x",) and p.default is not inspect.Parameter.empty}
    unknown = set(eval_kw) - set(resolved)
    if unknown:
        raise SystemExit(f"unknown eval arguments: {unknown}")

    env = {"python": sys.version.split()[0], "platform": platform.platform(), "cellpose": cellpose.version
           if hasattr(cellpose, "version") else __import__("importlib.metadata").metadata.version("cellpose"),
           "torch": torch.__version__, "cuda": torch.version.cuda, "numpy": np.__version__,
           "device": str(dev), "device_name": torch.cuda.get_device_name(0) if dev.type == "cuda" else platform.processor()}
    common = {"method": cfg.get("method", "cellpose"), "model": cfg["model"], "tag": tag, "config": str(args.config),
              "config_sha256": sha256(args.config), "overrides": overrides, "eval_resolved": resolved,
              "weights": str(weights), "weights_sha256": sha256(weights), "environment": env, "git": git_info(),
              "command": " ".join([sys.executable, *sys.argv])}
    print(f"[{datetime.now():%H:%M:%S}] {tag}: {len(crops)} crops -> {root}", flush=True)

    for i, (crop, gray_path, split) in enumerate(crops, 1):
        out = root / crop
        if (out / "run.json").exists():
            print(f"  [{i}/{len(crops)}] {crop}: done, skipped", flush=True)
            continue
        out.mkdir(parents=True, exist_ok=True)
        (out / "FAILED.json").unlink(missing_ok=True)
        start = datetime.now()
        try:
            kw = dict(eval_kw)
            if kw.get("diameter") == "2r_granulo":
                kw["diameter"] = 2.0 * r_granulo(crop, split)
            img = tifffile.imread(gray_path)
            if dev.type == "cuda":
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
            with RssPeak() as rss:
                t0 = time.perf_counter()
                masks, _flows, _styles = model.eval(img, **kw)
                runtime = time.perf_counter() - t0
            masks = np.asarray(masks)
            if masks.shape != img.shape:
                raise RuntimeError(f"mask shape {masks.shape} != image shape {img.shape}")
            n = int(masks.max())
            dtype = np.uint16 if n < 65536 else np.uint32
            tifffile.imwrite(out / "labels.tif", masks.astype(dtype), compression="zlib")
            ids, sizes = np.unique(masks[masks > 0], return_counts=True)
            rec = common | {
                "crop": crop, "input": str(gray_path), "input_shape": list(img.shape), "input_dtype": str(img.dtype),
                "split": split, "eval_used": kw, "output": str(out / "labels.tif"), "started": start.isoformat(timespec="seconds"),
                "runtime_seconds": round(runtime, 2),
                "peak_vram_allocated_mb": round(torch.cuda.max_memory_allocated() / 2**20, 1) if dev.type == "cuda" else None,
                "peak_vram_reserved_mb": round(torch.cuda.max_memory_reserved() / 2**20, 1) if dev.type == "cuda" else None,
                "peak_host_rss_mb": round(rss.peak / 2**20, 1),
                "n_labels": int(ids.size),
                "size_voxels": {q: float(np.percentile(sizes, p)) if sizes.size else None
                                for q, p in (("min", 0), ("p05", 5), ("median", 50), ("p95", 95), ("max", 100))},
                "n_labels_lt500": int((sizes < 500).sum()),
            }
            (out / "run.json").write_text(json.dumps(rec, indent=2), encoding="utf-8")
            print(f"  [{i}/{len(crops)}] {crop}: {ids.size} labels, median {rec['size_voxels']['median']}, "
                  f"{runtime:.0f} s, VRAM {rec['peak_vram_reserved_mb']} MB", flush=True)
        except Exception as e:  # recorded, scored later as an empty segmentation
            (out / "FAILED.json").write_text(json.dumps(common | {"crop": crop, "error": repr(e),
                                                                  "traceback": traceback.format_exc()}, indent=2),
                                             encoding="utf-8")
            print(f"  [{i}/{len(crops)}] {crop}: FAILED {e!r}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
