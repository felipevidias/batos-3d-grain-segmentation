"""BA-TOS v3 constants: every v2 constant of the pipeline, its dimension, and its v3 expression (single source of truth;
docs/V3_SCALE_FREE.md is generated from this table by `python -m batos_v3.constants --markdown`).

s = r_granulo / R_REF, with R_REF = 10 voxels = r_granulo of every one of the 12 dev crops (GT-free, step 1). So on the
dev crops s = 1 and every scaled constant equals its v2 value; constants are re-expressed in grain units, not re-tuned.

kind  meaning                                              v3 value
L     length (voxels)                                      v2 * s        (int flags: max(1, round(.)); 0 stays 0)
V     volume / ToS area (voxel count)                      v2 * s^3      (int flags: max(1, round(.)); 0 stays 0)
VF    voxel-degeneracy limit                               v2 * max(1, s^3)  (never below the voxel grid)
N     count of objects per crop (~ crop volume / grain)    v2 / s^3      (max(1, round(.)))
D     dimensionless (ratio, weight, score, fraction, topology)   kept
Q     rank / percentile of the image's own histogram       kept (already a statistic of the image)
R     intensity on the robust p0.5-p99.5 -> 0..255 scale  kept (already relative to the image's own robust range)
X     intensity guard removed (set to its non-binding limit; never or almost never binding on 115 dev+test runs)
T     replaced by a rule (see note)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

R_REF = 10.0


@dataclass(frozen=True)
class C:
    stage: str
    flag: str
    v2: object
    kind: str
    note: str = ""
    v3_fixed: object = None  # value for kind X


TABLE = [
    # 0 scale / 1 openrec
    C("scale", "(e13 r_est)", "median of top-1500 EDT 5^3-maxima", "T",
      "replaced by r_granulo (batos_v3.scale): argmax of the specimen-only pattern spectrum, r = 2..40"),
    C("openrec", "--radius", "max(8, round(1.5 r_e13))", "T", "max(8, round(1.5 r_granulo)) (floor 8 kept, user spec)"),
    C("openrec", "(robust p1/p99)", "1 / 99", "Q", "hard-coded rank normalisation; kept"),
    C("openrec", "(implementation)", "openrec_v2", "D", "openrec_lowmem (bit-identical) when radius > 15"),
    # 2 flatten (fixflat stage)
    C("flatten", "--lambda-xy", 1.0, "L", "per-slice Gaussian sigma of the background blur"),
    C("flatten", "--variant", "A", "D", "bg_ref = median(blur >= 1): one grey level, quantisation guard; kept"),
    C("flatten", "(robust p0.5/p99.5)", "0.5 / 99.5", "Q", "hard-coded; kept"),
    # 3 legacy ToS seed (run_batos_pipeline.sh step 2)
    C("tos_seed", "--height-factor", 0.012, "D", "fraction of the image range"),
    C("tos_seed", "--area-factor", 22, "D", "multiple of the mean node area"),
    C("tos_seed", "--area-min", 2, "V"),
    C("tos_seed", "--area-max", 2200, "V"),
    C("tos_seed", "--delta-parent-min", 0, "D"),
    C("tos_seed", "--dark-percentile-max", 78, "Q"),
    C("tos_seed", "--min-score", 0.35, "D"),
    C("tos_seed", "--max-selected-nodes", 900, "N"),
    C("tos_seed", "--support-area-min", 20, "V"),
    C("tos_seed", "--support-area-max", 45000, "V"),
    C("tos_seed", "--support-max-steps", 16, "D", "tree steps (intensity levels), not a length"),
    C("tos_seed", "--component-min", 2, "V"),
    C("tos_seed", "--component-max", 1500, "V"),
    # 4 background rejection (step 3)
    C("background", "--dark-percentile", 40, "T", "dropped: dark threshold = Otsu of the selector on specimen voxels"),
    C("background", "--otsu-offset", -10, "T", "dropped (idem)"),
    C("background", "--dark-threshold-max", 135, "T", "dropped (idem); passed as --dark-threshold"),
    C("background", "--dark-connectivity", 1, "D"),
    C("background", "--seed-area-min", 1, "V"),
    C("background", "--seed-area-max", 3000, "V"),
    C("background", "--max-external-overlap", 0.10, "D"),
    C("background", "--min-internal-overlap", 0.00, "D"),
    C("background", "--internal-search-radius", 6, "L"),
    C("background", "--center-radius", 2, "L", "script default, now passed"),
    C("background", "--ring-inner", 4, "L", "script default, now passed"),
    C("background", "--ring-outer", 10, "L", "script default, now passed"),
    C("background", "--center-dark-max", 210.0, "X", "0 rejections in 115 runs", 255.0),
    C("background", "--min-ring-contrast", -5.0, "X", "0 rejections in 115 runs", -255.0),
    C("background", "--max-ring-external-coverage", 0.70, "D"),
    C("background", "--nms-min-dist", 5.0, "L"),
    # 5 legacy markers (step 4; used only for the watershed mask union)
    C("markers_legacy", "--connectivity", 1, "D"),
    C("markers_legacy", "--area-min", 3, "V"),
    C("markers_legacy", "--area-max", 900, "V"),
    C("markers_legacy", "--max-dim", 30, "L"),
    C("markers_legacy", "--max-elongation", 8.0, "D"),
    C("markers_legacy", "--min-fill", 0.005, "D"),
    C("markers_legacy", "--center-radius", 2, "L", "script default, now passed"),
    C("markers_legacy", "--ring-inner", 4, "L", "script default, now passed"),
    C("markers_legacy", "--ring-outer", 10, "L", "script default, now passed"),
    C("markers_legacy", "--center-dark-max", 205.0, "X", "3 rejections in 1.93 M components (115 runs)", 255.0),
    C("markers_legacy", "--min-ring-contrast", -2.0, "R", "binding (130k rejections): sign test with a 2-level "
      "(0.8 % of the robust range) tolerance; kept"),
    C("markers_legacy", "--nms-min-dist", 7.0, "L"),
    C("markers_legacy", "--max-markers", 420, "N"),
    # 6 true-dual selector (K = 3600)
    C("true_dual", "--height-factor", 0.012, "D"),
    C("true_dual", "--area-factor", 22.0, "D"),
    C("true_dual", "--area-min", 2, "V"),
    C("true_dual", "--area-max", 2200, "V"),
    C("true_dual", "--delta-parent-min", 0, "D"),
    C("true_dual", "--dark-percentile-max", 78.0, "Q"),
    C("true_dual", "--support-area-min", 20, "V"),
    C("true_dual", "--support-area-max", 45000, "V"),
    C("true_dual", "--support-max-steps", 16, "D"),
    C("true_dual", "--support-area-ratio-preferred", 18.0, "D", "volume ratio"),
    C("true_dual", "--support-area-ratio-width", 80.0, "D", "volume ratio"),
    C("true_dual", "--support-contrast-norm", 80.0, "R", "score normaliser in grey levels"),
    C("true_dual", "--delta-norm", 20.0, "R", "score normaliser in grey levels"),
    C("true_dual", "--area-preferred-position", 0.42, "D", "position on log(area_min..area_max), both scaled"),
    C("true_dual", "--area-preferred-width", 0.55, "D"),
    C("true_dual", "--w-height", 0.22, "D"),
    C("true_dual", "--w-delta", 0.12, "D"),
    C("true_dual", "--w-area", 0.18, "D"),
    C("true_dual", "--w-dark", 0.18, "D"),
    C("true_dual", "--w-support", 0.3, "D"),
    C("true_dual", "--min-score", 0.35, "D", "tau_s"),
    C("true_dual", "--cleanup-open-radius", 0, "L"),
    C("true_dual", "--cleanup-close-radius", 0, "L"),
    C("true_dual", "--component-min", 2, "V"),
    C("true_dual", "--component-max", 1500, "V"),
    C("true_dual", "--max-selected-nodes", 3600, "N", "K"),
    # 7 e11 ToS-MIN markers
    C("markers_tos", "--polarity", "min", "D"),
    C("markers_tos", "--nms-radius", 6.0, "L"),
    C("markers_tos", "--edt-priority-weight", 0.20, "D", "weight of log1p(EDT)/log1p(max EDT)"),
    C("markers_tos", "--min-priority-bonus", 0.05, "D"),
    C("markers_tos", "(no MIN candidate in internal_dark)", "RuntimeError", "T", "recorded empty result (F1 0)"),
    # 8 EDT-Otsu control (gate fallback)
    C("edt_otsu", "--mode", "otsu_edt", "D"),
    C("edt_otsu", "--nms-radius", 6.0, "L"),
    C("edt_otsu", "--v3-min-component", 500, "V", "hard-coded in the frozen script; patched copy"),
    C("edt_otsu", "--v3-closing-iterations", 1, "L", "hard-coded (6-conn closing); patched copy"),
    # 9 fixed-mask watershed (ToS route and selective route)
    C("ws_tos", "--outside-policy", "drop", "D"),
    C("ws_tos", "--v3-min-component", 500, "V", "mask components >= 500; hard-coded; patched copy"),
    C("ws_tos", "--v3-closing-iterations", 1, "L", "mask closing; hard-coded; patched copy"),
    C("ws_sel", "--outside-policy", "drop", "D"),
    C("ws_sel", "--v3-min-component", 500, "V", "idem"),
    C("ws_sel", "--v3-closing-iterations", 1, "L", "idem"),
    # 10 E23c gate
    C("route", "--v3-median-limit", 3.0, "VF", "gate median(region_voxels) <= 3; hard-coded; patched copy. A collapsed "
      "support is a 1-3 voxel discretisation degeneracy: grows with s^3 for s > 1, never below 3"),
    C("route", "--v3-q75-limit", 4.0, "VF", "gate q75(region_voxels) <= 4; idem"),
    # 11 post-processing
    C("post", "size_filter", 500, "V", "prediction size filter (reported alongside filter off)"),
]

INT_FLAGS = {"--area-min", "--area-max", "--max-selected-nodes", "--support-area-min", "--support-area-max",
             "--component-min", "--component-max", "--seed-area-min", "--seed-area-max", "--internal-search-radius",
             "--center-radius", "--ring-inner", "--ring-outer", "--max-dim", "--max-markers", "--cleanup-open-radius",
             "--cleanup-close-radius", "--v3-min-component", "--v3-closing-iterations", "size_filter"}


def scaled(c: C, s: float):
    if c.kind == "X":
        return c.v3_fixed
    if c.kind in ("D", "Q", "R", "T"):
        return c.v2
    v = float(c.v2)
    if v == 0:
        return c.v2
    x = {"L": v * s, "V": v * s ** 3, "VF": v * max(1.0, s ** 3), "N": v / s ** 3}[c.kind]
    if c.flag in INT_FLAGS or c.kind == "N":
        return max(1, int(round(x)))
    return round(x, 4)


def stage_flags(s: float) -> dict[str, list[str]]:
    """v3 CLI flags per stage (flags of kind T and pseudo-flags '(...)' are handled by the driver)."""
    out: dict[str, list[str]] = {}
    for c in TABLE:
        if c.kind == "T" or c.flag.startswith("(") or c.stage in ("scale", "openrec", "post"):
            continue
        out.setdefault(c.stage, []).extend([c.flag, str(scaled(c, s))])
    return out


def size_filter(s: float) -> int:
    return scaled(next(c for c in TABLE if c.flag == "size_filter"), s)


def markdown(example_s=(0.7, 1.0, 1.4, 2.0)) -> str:
    head = "| stage | constant | v2 | kind | v3 rule | " + " | ".join(f"s={x:g}" for x in example_s) + " | note |"
    lines = [head, "|" + "---|" * (6 + len(example_s))]
    rule = {"L": "v2·s", "V": "v2·s³", "VF": "v2·max(1, s³)", "N": "v2/s³", "D": "kept", "Q": "kept", "R": "kept", "X": "removed", "T": "rule"}
    for c in TABLE:
        vals = " | ".join(str(scaled(c, x)) if c.kind not in ("T",) else "—" for x in example_s)
        lines.append(f"| {c.stage} | `{c.flag}` | {c.v2} | {c.kind} | {rule[c.kind]} | {vals} | {c.note} |")
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--s", type=float, default=1.0)
    a = ap.parse_args()
    print(markdown() if a.markdown else stage_flags(a.s))
