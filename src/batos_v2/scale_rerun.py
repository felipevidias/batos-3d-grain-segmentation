"""Phase 5b follow-up B: frozen pipeline re-run with openrec radius = max(8, round(1.5 * r_granulo)).

Everything is the frozen pipeline (batos_v2.sensitivity.CropRun, all stages from openrec downstream re-run; frozen
flags; evaluate_v2 with markers, seed-ignore, no size filter) except:
  * the scale estimate: r_granulo (results/scale/<crop>_granulo.json, key --scale-key, default r_ps_peak) replaces
    e13's r_est in the openrec radius formula max(8, round(1.5 r));
  * the openrec stage runs batos_v2.openrec_lowmem, verified bit-identical to the frozen script at r = 15 on
    EFRGP01_00_crop04 (all four outputs), so radii >= 20 fit in < 1 GB instead of ~18 GB.
All stages are re-run even when the radius equals the frozen one (identity check of the substitution).
Outputs: C:\\BA_TOS_DATA\\runs\\scale_granulo\\<key>\\<crop>\\, results/scale/rerun_summary.csv, results/all_runs.csv.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from .common import REPO_ROOT
from .extinction_run import TrackedRun, summarize, upsert
from .run_crop import git_state
from .sensitivity import RamBudget

ROOT = Path(r"C:\BA_TOS_DATA\runs\scale_granulo")


class ScaleRun(TrackedRun):
    r_granulo: float = 0.0

    def frozen_paths(self):
        paths = super().frozen_paths()
        self.frozen_r_est = float(self.auto["estimated_grain_radius"])
        self.auto = dict(self.auto, estimated_grain_radius=float(self.r_granulo))
        self.frozen_radius_original = self.frozen_radius
        self.frozen_radius = -1          # force the full re-run from openrec downstream
        return paths

    def _openrec(self, args):
        return [sys.executable, "-I", "-m", "batos_v2.openrec_lowmem", *map(str, args)], "batos_v2.openrec_lowmem"

    replace = {"openrec": _openrec}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="frozen pipeline with granulometric openrec radius")
    ap.add_argument("--crops", nargs="*", default=["EFRGP01_00_crop04", "EFRGP01_00_crop05", "EFRGP01_01_crop05"])
    ap.add_argument("--scale-key", default="r_ps_peak")
    ap.add_argument("--ram-gb", type=float, default=6.0)
    args = ap.parse_args(argv)
    git = git_state()
    rows, runs = [], []
    for crop in args.crops:
        r = float(json.loads((REPO_ROOT / "results" / "scale" / f"{crop}_granulo.json").read_text(
            encoding="utf-8"))[args.scale_key])
        out = ROOT / args.scale_key / crop
        label = f"frozen_openrec_from_{args.scale_key}"
        if not (out / "DONE").exists():
            if out.exists():
                raise RuntimeError(f"partial run folder exists, inspect and remove it first: {out}")
            param = {"name": "openrec_radius_factor", "group": "scale", "kind": "scale", "frozen": 1.5,
                     "type": "float"}
            run = ScaleRun(crop, param, "r_granulo", 1.5, out, RamBudget(args.ram_gb), keep_tifs=True)
            run.r_granulo = r
            run.note = (f"openrec radius from {args.scale_key}={r} (frozen e13 r_est replaced); "
                        "openrec via batos_v2.openrec_lowmem (bit-identical to frozen at r=15)")
            run.execute()
            rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
            rec |= {"git": git, "r_granulo": r, "frozen_r_est": run.frozen_r_est,
                    "frozen_openrec_radius": run.frozen_radius_original}
            (out / "run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
        rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
        r1, r2 = summarize(out, crop, label, None, git)
        for x in r1:
            x |= {"r_granulo": rec["r_granulo"], "frozen_r_est": rec["frozen_r_est"],
                  "frozen_openrec_radius": rec["frozen_openrec_radius"]}
        rows += r1
        runs += r2
    upsert(REPO_ROOT / "results" / "scale" / "rerun_summary.csv", rows, ["crop_id", "variant", "route"])
    upsert(REPO_ROOT / "results" / "all_runs.csv", runs, ["experiment_id", "method"])
    print(pd.DataFrame(rows)[["crop_id", "route", "r_granulo", "openrec_radius", "frozen_openrec_radius",
                              "n_markers", "grains_seeded", "grains_split", "gate_triggered", "gate_mv", "gate_qv",
                              "f1_050", "f1_075", "max_stage_peak_mb"]].to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
