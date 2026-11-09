"""BA-TOS v3.1 single-crop driver: the frozen v3 pipeline (batos_v3.pipeline.run_crop, preset v3) with the v3.1
count rule (batos_v31.constants). Nothing else differs; no v3 file is modified (v3 is SHA-256 frozen).

Mechanism: batos_v3.pipeline reads its flags through its module-level name `constants`; during the call that name is
bound to batos_v31.constants and restored afterwards. The run record (v3_run.json) is then re-written with
preset 'v3.1' and the N values used.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import batos_v3.pipeline as P
from batos_v31 import constants as C31
from batos_v3.pipeline import StageError, outputs  # noqa: F401  (re-exported for drivers)


def run_crop(crop_id: str, gray: Path, out: Path) -> dict:
    orig = P.constants
    P.constants = C31
    try:
        rec = P.run_crop(crop_id, gray, out, "v3")
    finally:
        P.constants = orig
    rec["preset"] = "v3.1"
    rec["n_rule"] = "v2 * max(1, s^-3)"
    rec["n_values"] = {f"{c.stage} {c.flag}": C31.scaled(c, rec["s"]) for c in C31.N_CONSTANTS}
    (Path(out) / "v3_run.json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--crop-id", required=True)
    ap.add_argument("--gray", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    rec = run_crop(a.crop_id, a.gray, a.out)
    print(json.dumps({k: rec.get(k) for k in ("crop_id", "preset", "s", "n_values", "gate", "tos_markers",
                                              "empty_result", "seconds_total")}, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
