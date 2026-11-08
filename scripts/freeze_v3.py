#!/usr/bin/env python3
"""Phase 5e step 4: write configs/batos_paper_v3_scalefree.yaml (status frozen_for_test) with the SHA-256 of every file
the v3 preset executes (own code + recovered frozen scripts), the dev / synthetic evidence, and the v3 constants table.
Run once, after the step-3 commit; then commit the config and tag batos-paper-v3-scalefree. run_test_v3.py refuses
to start if any recorded hash differs.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from batos_v2.frozen import SCRIPTS  # noqa: E402
from batos_v2.run_crop import git_state  # noqa: E402
from batos_v3 import constants  # noqa: E402

OWN = ["src/batos_v3/__init__.py", "src/batos_v3/constants.py", "src/batos_v3/scale.py", "src/batos_v3/pipeline.py",
       "src/batos_v3/make_patched.py", "src/batos_v3/patched/e11_build_control_markers_v3.py",
       "src/batos_v3/patched/e11_run_fixed_mask_watershed_v3.py",
       "src/batos_v3/patched/e23c_build_selective_marker_route_v3.py", "src/batos_v3/patched/patched_manifest.json",
       "src/batos_v2/flatten_fixed.py", "src/batos_v2/openrec_lowmem.py", "src/batos_v2/extinction_run.py",
       "src/batos_v2/frozen.py", "src/batos_v2/sensitivity.py", "src/batos_v2/run_crop.py",
       "src/batos_v2/evaluate_v2.py", "scripts/run_test_v3.py", "scripts/aggregate_test_v3.py"]
RECOVERED = ["openrec_v2", "batos_tos_seed", "batos_background", "batos_markers", "true_dual", "spatialize",
             "e11_markers"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    g = git_state()
    if g["branch"] != "journal-extension":
        raise RuntimeError(g)
    cfg = {
        "preset": "batos-paper-v3-scalefree",
        "git_tag": "batos-paper-v3-scalefree",
        "git_commit_code": g["commit"],
        "frozen_on": dt.date.today().isoformat(),
        "status": "frozen_for_test",
        "method": {
            "base": "batos-paper-v2-fixflat (stage order and frozen scripts of batos_v2.sensitivity.CropRun)",
            "scale": "r_granulo = argmax of the specimen-only binary pattern spectrum, r = 2..40 (batos_v3.scale)",
            "s": "r_granulo / R_REF",
            "R_REF": constants.R_REF,
            "openrec_radius": "max(8, round(1.5 r_granulo)); openrec_lowmem (bit-identical) above 15",
            "stage03_dark_threshold": "Otsu of the selector on specimen voxels (--dark-threshold)",
            "empty_rule": "e11 'no MIN candidate inside internal_dark' -> recorded empty result (F1 0)",
            "gate": "E23c median(region_voxels) <= 3 max(1, s^3) and q75 <= 4 max(1, s^3) -> EDT-Otsu",
            "size_filter": "500 s^3 (reported alongside filter off and 500)",
            "constants": "src/batos_v3/constants.py TABLE (kinds L/V/VF/N/D/Q/R/X/T); docs/V3_SCALE_FREE.md",
        },
        "evaluation": {"protocol": "evaluate_v2 markers_dominant (seed-ignore, paper Sec. III-C), primary",
                       "secondary": "no_markers_fraction050", "iou": [0.25, 0.5, 0.75],
                       "gt_valid": ">= 1000 voxels, not touching a crop face",
                       "size_filter": ["none (primary)", "500", "own (500 s^3)"], "routes": ["BA-TOS_ToS", "BA-TOS_sel"]},
        "design_data": {"dev": "12 SIBGRAPI EFR crops (data/metadata/crops_sibgrapi.json)",
                        "synthetic": "36 rescaled dev crops f 0.7/1.4/2.0 (data/metadata/crops_synth_scale.json)",
                        "evidence": "results/v3_scale_test/REPORT.md",
                        "test_untouched_until_freeze": True},
        "drivers": {"test": "python -I scripts/run_test_v3.py", "aggregate": "python -I scripts/aggregate_test_v3.py",
                    "single_crop": "python -I -m batos_v3.pipeline --preset v3 --crop-id <id> --gray <tif> --out <dir>"},
        "constants_at_s1": constants.stage_flags(1.0),
        "code_sha256_at_freeze": {p: sha(REPO / p) for p in OWN},
        "recovered_scripts_sha256": {k: {"path": str(SCRIPTS[k]), "sha256": sha(Path(SCRIPTS[k]))} for k in RECOVERED},
    }
    head = ("# batos-paper-v3-scalefree — FROZEN method under test (Phase 5e step 4).\n"
            "# Designed on the 12 EFR dev crops + synthetic scale test only. No parameter, code or threshold may change\n"
            "# before the 52-crop test results are reported. scripts/run_test_v3.py checks every hash below.\n")
    out = REPO / "configs" / "batos_paper_v3_scalefree.yaml"
    out.write_text(head + yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True, width=120), encoding="utf-8")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
