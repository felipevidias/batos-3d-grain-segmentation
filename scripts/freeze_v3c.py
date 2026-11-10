#!/usr/bin/env python3
"""Phase 5i: write configs/batos_paper_v3c.yaml (status frozen_for_test) with the SHA-256 of every file the test-set-2
run executes or imports (v3c code, the frozen v3 / v2-fixflat code it wraps, the test-2 definition / cut / run /
aggregation drivers) plus the dev / synthetic evidence. Run once, only if both acceptance criteria of
docs/V3C_CANONICAL.md hold; then commit the config and tag batos-paper-v3c. scripts/run_test2_v3c.py refuses to start if
any recorded hash differs.

Hashes are taken over LF-normalised bytes (CRLF -> LF): with core.autocrlf the checkout line endings are not part of the
code (the v3 freeze recorded LF bytes and fails a raw-byte check on a CRLF checkout of identical content).
"""

from __future__ import annotations

import datetime as dt
import hashlib
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.frozen import SCRIPTS  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402

_spec = importlib.util.spec_from_file_location("freeze_v3", REPO / "scripts" / "freeze_v3.py")
fv3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fv3)

OWN = fv3.OWN + [
    "src/batos_v3c/__init__.py", "src/batos_v3c/canonical.py", "src/batos_v3c/eval_lean.py",
    "scripts/run_v3c_dev.py", "scripts/check_eval_lean.py", "scripts/aggregate_v3c_dev.py",
    "scripts/make_synthetic_canon.py", "scripts/define_test_crops.py", "scripts/define_test_crops2.py",
    "scripts/extract_and_cut.py", "scripts/make_crops.py", "scripts/cut_test2.py", "scripts/cut_test2_v3c.py",
    "scripts/run_test2_v3c.py", "scripts/aggregate_test2_v3c.py", "configs/batos_paper_v3_scalefree.yaml",
    "docs/V3C_CANONICAL.md",
]
CONFIG = REPO / "configs" / "batos_paper_v3c.yaml"


def sha_lf(p: Path) -> str:
    return hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main() -> int:
    if CONFIG.exists():
        raise FileExistsError(f"{CONFIG} exists; v3c is frozen once")
    g = git_state()
    if g["branch"] != "journal-extension-v32" or g.get("dirty"):
        raise RuntimeError(f"freeze from a clean journal-extension-v32 tree: {g}")
    acc = pd.read_csv(REPO / "results" / "v3c_dev" / "acceptance.csv")
    if len(acc) != 4 or not (acc["pass"].astype(str) == "True").all():
        raise RuntimeError(f"acceptance not met; v3c is not frozen:\n{acc}")
    cfg = {
        "preset": "batos-paper-v3c",
        "git_tag": "batos-paper-v3c",
        "git_commit_code": g["commit"],
        "frozen_on": dt.date.today().isoformat(),
        "status": "frozen_for_test",
        "method": {
            "idea": "canonical-scale normalisation around the frozen batos-paper-v2-fixflat pipeline (docs/V3C_CANONICAL.md)",
            "scale": "r_granulo = batos_v3.scale.estimate on the 200^3 probe (only scale input); f = 10 / r_granulo",
            "window": "native side round(200 r_granulo / 10) at the probe centre (canonical field of view 200^3)",
            "resample_in": "scipy.ndimage.zoom order 3, grid_mode, grid-mirror; clipped, rounded, input dtype; none at f = 1",
            "canonical_run": "batos_v3.pipeline.run_crop preset fixflat, e13 replaced by r = 10, openrec 15; all v2 constants",
            "size_filter_own": "500 canonical voxels",
            "map_back": "labels nearest neighbour; markers by voxel-centre coordinate map",
            "identity": "f = 1 reproduces batos-paper-v2-fixflat bit-for-bit (EFRGP01_00_crop04, both routes)",
        },
        "evaluation": {"evaluator": "batos_v3c.eval_lean (exactly equal to batos_v2.evaluate_v2; scripts/check_eval_lean.py)",
                       "protocol": "markers_dominant (seed-ignore) primary; no_markers_fraction050 secondary; IoU .25/.50/.75",
                       "gt_valid": ">= round(1000 (r_granulo/10)^3) native voxels (1000 canonical), not touching a face"},
        "acceptance_dev": acc.to_dict(orient="records"),
        "design_data": {"dev": "12 SIBGRAPI EFR crops", "synthetic": "36 canonical-FOV windows f 0.7/1.4/2.0 "
                        "(data/metadata/crops_synth_canon.json)", "evidence": "results/v3c_dev/tables.md",
                        "test_set_1_used_for": "nothing in the v3c design", "test_set_2_defined_after_freeze": True},
        "test_set_2": {"definition": "python -I scripts/define_test_crops2.py (seed 20261109)",
                       "cut": "python -I scripts/cut_test2_v3c.py", "run": "python -I scripts/run_test2_v3c.py",
                       "aggregate": "python -I scripts/aggregate_test2_v3c.py", "methods": ["v3", "v3c"]},
        "hash_convention": "sha256 of LF-normalised bytes",
        "code_sha256_at_freeze": {p: sha_lf(REPO / p) for p in OWN},
        "recovered_scripts_sha256": {k: {"path": str(SCRIPTS[k]), "sha256": sha_lf(Path(SCRIPTS[k]))}
                                     for k in fv3.RECOVERED + ["openrec_lowmem"] if k in SCRIPTS},
    }
    head = ("# batos-paper-v3c — FROZEN method under test (Phase 5i).\n"
            "# Designed on the 12 EFR dev crops + canonical-FOV synthetic scale set only. Test set 2 is defined after\n"
            "# this freeze. scripts/run_test2_v3c.py checks every hash below (LF-normalised) before running.\n")
    CONFIG.write_text(head + yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True, width=120), encoding="utf-8")
    print(CONFIG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
