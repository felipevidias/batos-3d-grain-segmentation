"""Root-cause fix test: re-run the frozen pipeline from the flatten stage onward with batos_v2.flatten_fixed.

Upstream (scale, openrec) is read from the frozen run. The flatten stage is the copy batos_v2.flatten_fixed
(variant A: bg_ref = median(blur >= 1); variant B: median of the blur over e13's Otsu foreground of the raw crop).
Every later stage (legacy ToS seed, background, legacy markers, true_dual, e10, e11 markers, EDT-Otsu, watershed, E23c
gate, selective watershed, evaluate_v2) is the unmodified frozen call. Outputs:
C:\\BA_TOS_DATA\\runs\\flatten_fix\\<variant>\\<crop>\\, results/flatten_fix/summary.csv, results/all_runs.csv.
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
from .sensitivity import FROZEN, RamBudget

ROOT = Path(r"C:\BA_TOS_DATA\runs\flatten_fix")


class FlattenFixRun(TrackedRun):
    raw_gray: str = ""

    def _flatten(self, args):
        cmd = [sys.executable, "-I", "-m", "batos_v2.flatten_fixed", *map(str, args)]
        if self.value == "B":
            cmd += ["--raw-gray", self.raw_gray]
        return cmd, "batos_v2.flatten_fixed"

    replace = {"flatten": _flatten}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="flatten bg_ref fix test")
    ap.add_argument("--crops", nargs="*", default=["EFRGP01_00_crop05", "EFRGP01_01_crop05", "EFRGP01_00_crop04"])
    ap.add_argument("--variants", nargs="*", default=["A", "B"], choices=["A", "B", "frozen"])
    ap.add_argument("--ram-gb", type=float, default=6.0)
    args = ap.parse_args(argv)
    git = git_state()
    rows, runs = [], []
    for v in args.variants:
        for crop in args.crops:
            out = ROOT / v / crop
            label = f"flatten_fix_{v}"
            if not (out / "DONE").exists():
                if out.exists():
                    raise RuntimeError(f"partial run folder exists, inspect and remove it first: {out}")
                rec0 = json.loads((FROZEN / crop / "frozen_run.json").read_text(encoding="utf-8"))
                param = {"name": label, "group": "flatten_fix", "kind": "flag", "stage": "flatten",
                         "flag": "--variant", "frozen": "frozen", "type": "str"}
                run = FlattenFixRun(crop, param, "fix", v, out, RamBudget(args.ram_gb), keep_tifs=True)
                run.raw_gray = rec0["inputs"]["gray"]
                run.note = f"flatten stage = batos_v2.flatten_fixed --variant {v}; everything else frozen"
                run.execute()
                rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
                info = json.loads((out / "pipe" / "01_preprocess" / "flatten_fixed_info.json").read_text(
                    encoding="utf-8"))
                rec |= {"git": git, "bg_ref": info["bg_ref"]}
                (out / "run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
            rec = json.loads((out / "run.json").read_text(encoding="utf-8"))
            r1, r2 = summarize(out, crop, label, None, git)
            for x in r1:
                x["bg_ref"] = rec.get("bg_ref")
            rows += r1
            runs += r2
    upsert(REPO_ROOT / "results" / "flatten_fix" / "summary.csv", rows, ["crop_id", "variant", "route"])
    upsert(REPO_ROOT / "results" / "all_runs.csv", runs, ["experiment_id", "method"])
    print(pd.DataFrame(rows).reindex(columns=["crop_id", "variant", "route", "bg_ref", "n_markers", "grains_seeded",
                                              "grains_split", "gate_triggered", "gate_mv", "gate_qv", "f1_050",
                                              "f1_075", "max_stage_peak_mb"]).to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
