"""BA-TOS v3.1 constants: identical to batos_v3.constants except kind N (count of objects per crop).

v3   N: v2 / s^3            (shrinks the candidate / marker capacity for large grains)
v3.1 N: v2 * max(1, s^-3)   (never below the v2 capacity; grows for small grains as in v3)

Why (Phase 5f, docs/EXPERIMENT_LOG.md 2026-11-09): on held-out test set 1 the MIN candidate count hit the v3 cap
1800/s^3 on every crop, and EO interior crops (s = 1.5) were left with 1.19 candidates / 0.35 ToS markers per grain.
Grain counts per crop do not fall as s^-3: on the synthetic scale set they fall as ~s^-2.4 (face-touching grains), and
grain shape changes the count at a given r_granulo. A capacity floor removes the assumption instead of re-fitting it.
Dev evidence that extra capacity is not harmful at s = 1: K = 7200 on EFRGP01_00_crop04 gave F1@0.50 0.998 (Table IV).

N constants (all three changed): tos_seed --max-selected-nodes 900; markers_legacy --max-markers 420;
true_dual --max-selected-nodes 3600 (= K, 1800 MIN + 1800 MAX: one flag covers both pools).
For s <= 1 (all dev crops, s = 1; synthetic f = 0.7) every flag equals v3's.
"""

from __future__ import annotations

from batos_v3 import constants as v3
from batos_v3.constants import INT_FLAGS, R_REF, TABLE, C  # noqa: F401  (re-exported)

N_CONSTANTS = [c for c in TABLE if c.kind == "N"]


def scaled(c: C, s: float):
    if c.kind != "N":
        return v3.scaled(c, s)
    return max(1, int(round(float(c.v2) * max(1.0, s ** -3))))


def stage_flags(s: float) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for c in TABLE:
        if c.kind == "T" or c.flag.startswith("(") or c.stage in ("scale", "openrec", "post"):
            continue
        out.setdefault(c.stage, []).extend([c.flag, str(scaled(c, s))])
    return out


def size_filter(s: float) -> int:
    return v3.size_filter(s)


if __name__ == "__main__":
    for x in (0.7, 1.0, 1.4, 1.5, 2.0):
        print(x, {f"{c.stage} {c.flag}": (v3.scaled(c, x), scaled(c, x)) for c in N_CONSTANTS})
