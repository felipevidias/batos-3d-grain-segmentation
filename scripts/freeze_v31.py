#!/usr/bin/env python3
"""Phase 5f step 3: write configs/batos_paper_v31.yaml (status frozen_for_test) with the SHA-256 of every file the
test-set-2 run executes or imports (v3.1 code, the frozen v3 code it wraps, the test-1 drivers it re-points, the
test-2 drivers and aggregator) plus the dev / synthetic evidence. Run once after the step-2 commit; then commit the
config and tag batos-paper-v3.1. scripts/run_test2.py refuses to start if any recorded hash differs.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import importlib.util
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.frozen import SCRIPTS  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402
from batos_v31 import constants as c31  # noqa: E402

_spec = importlib.util.spec_from_file_location("freeze_v3", REPO / "scripts" / "freeze_v3.py")
fv3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fv3)

OWN = fv3.OWN + [
    "src/batos_v31/__init__.py", "src/batos_v31/constants.py", "src/batos_v31/pipeline.py",
    "src/batos_v2/flatten_fix_run.py", "src/batos_v2/scale_granulo.py",
    "scripts/run_test_protocol.py", "scripts/define_test_crops.py", "scripts/extract_and_cut.py", "scripts/make_crops.py",
    "scripts/define_test_crops2.py", "scripts/cut_test2.py", "scripts/run_test2.py", "scripts/aggregate_test2.py",
    "scripts/freeze_v31.py", "configs/batos_paper_v3_scalefree.yaml",
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    g = git_state()
    if g["branch"] != "journal-extension":
        raise RuntimeError(f"freeze from journal-extension: {g}")
    # The tree may be dirty only through the user's uncommitted sweep table (not code); the SHA-256 list guards the code.
    cfg = {
        "preset": "batos-paper-v3.1",
        "git_tag": "batos-paper-v3.1",
        "git_commit_code": g["commit"],
        "frozen_on": dt.date.today().isoformat(),
        "status": "frozen_for_test",
        "method": {
            "base": "batos-paper-v3-scalefree (configs/batos_paper_v3_scalefree.yaml), unchanged except kind N",
            "change": "count constants (kind N) = v2 * max(1, s^-3) instead of v2 / s^3 (never below the v2 capacity)",
            "n_constants": {f"{c.stage} {c.flag}": c.v2 for c in c31.N_CONSTANTS},
            "n_values_example": {f"s={s:g}": {f"{c.stage} {c.flag}": c31.scaled(c, s) for c in c31.N_CONSTANTS}
                                 for s in (0.7, 1.0, 1.4, 1.5, 2.0)},
            "identical_to_v3_for": "s <= 1 (all dev crops; synthetic f = 0.7)",
            "mechanism": "batos_v31.pipeline.run_crop = batos_v3.pipeline.run_crop with batos_v31.constants",
        },
        "design_data": {"dev": "12 SIBGRAPI EFR crops", "synthetic": "36 rescaled dev crops f 0.7/1.4/2.0",
                        "evidence": "results/v31_dev/tables.md, REPORT.md", "test_set_1_used_for": "diagnosis only "
                        "(results/test_v3/REPORT.md section 7b); no parameter fitted on it",
                        "test_set_2_defined_after_freeze": True},
        "test_set_2": {"definition": "python -I scripts/define_test_crops2.py (seed 20261109)",
                       "cut": "python -I scripts/cut_test2.py", "run": "python -I scripts/run_test2.py",
                       "aggregate": "python -I scripts/aggregate_test2.py",
                       "methods": ["frozen", "fixflat", "v3", "v3.1"]},
        "constants_at_s1": c31.stage_flags(1.0),
        "code_sha256_at_freeze": {p: sha(REPO / p) for p in OWN if p != "scripts/freeze_v31.py"},
        "recovered_scripts_sha256": {k: {"path": str(SCRIPTS[k]), "sha256": sha(Path(SCRIPTS[k]))}
                                     for k in fv3.RECOVERED + ["e13_estimate_scale"]},
    }
    head = ("# batos-paper-v3.1 — FROZEN method under test (Phase 5f step 3).\n"
            "# Designed on the 12 EFR dev crops + synthetic scale set only. Test set 2 is defined after this freeze.\n"
            "# scripts/run_test2.py checks every hash below (and the v3 config's) before running.\n")
    out = REPO / "configs" / "batos_paper_v31.yaml"
    out.write_text(head + yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True, width=120), encoding="utf-8")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
