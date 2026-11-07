#!/usr/bin/env python3
"""Phase 1f forensics: per-stage GT-free quantities for every test crop and method, from the CSV/JSON intermediates
kept after pruning (C:\\BA_TOS_DATA\\runs\\test_protocol\\{frozen,fixflat}\\<crop>). GT enters only via the F1 already
computed by evaluate_v2 and the GT inscribed radius (diagnosis only; nothing here feeds back into the method).

Stages, in pipeline order:
  0 scale     e13 r_est, openrec radius (frozen_run.json) vs granulometric r_ps_peak and GT inscribed radius
  1 flatten   bg_ref (fixflat only)
  2 tos_seed  number of seed components (bgreject_report rows)
  3 background seed components kept (reason == keep) vs rejected as external_background; from the stage log:
              otsu, p40, dark_threshold = min(otsu-10, p40, 135), dark / external / internal voxels, threshold branch
  4 e10       internal_components (spatialization_summary.json)
  5 e11       ToS-MIN markers, EDT-Otsu markers
  6 gate      mv / qv / triggered
  7 result    F1@0.50 ToS / sel (markers/dominant, no filter)
Output: results/test_protocol/forensics_stages.csv (one row per crop x method).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
ROOT = Path(r"C:\BA_TOS_DATA\runs\test_protocol")
OUT = REPO / "results" / "test_protocol"


def one(base: Path, pattern: str) -> Path | None:
    hits = sorted(base.glob(pattern))
    return hits[0] if len(hits) == 1 else None


def nrows(p: Path | None) -> int | None:
    return None if p is None or not p.exists() else int(len(pd.read_csv(p)))


def stages(method: str, crop: str) -> dict:
    base = ROOT / method / crop
    pipe = base / ("batos_openrec_auto_minparam" if method == "frozen" else "pipe")
    r = {"crop_id": crop, "method": method, "crashed": (base / "FAILED.json").exists()}
    bg = one(pipe / "03_background", "*_bgreject_report.csv")
    if bg is not None:
        d = pd.read_csv(bg)
        r |= {"seed_components": len(d), "bg_keep": int((d.reason == "keep").sum()),
              "bg_external": int((d.reason == "external_background").sum()),
              "bg_other": int((~d.reason.isin(["keep", "external_background"])).sum()),
              "seed_area_median": float(d.seed_area.median())}
    logf = base / ("frozen_commands.log" if method == "frozen" else "commands.log")
    if logf.exists():
        txt = logf.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"otsu=(\d+), percentile40\.0=([\d.]+), dark_threshold=(\d+)", txt)
        if m:
            otsu, p40, thr = int(m[1]), float(m[2]), int(m[3])
            r |= {"bg_otsu": otsu, "bg_p40": p40, "bg_dark_thr": thr,
                  "bg_thr_branch": "otsu-10" if thr == otsu - 10 else ("p40" if thr == int(p40) else "cap135")}
            for k in ("dark voxels", "external dark", "internal dark"):
                v = re.search(rf"{k}\s+=\s+(\d+)", txt)
                r[f"bg_{k.replace(' ', '_')}"] = int(v[1]) if v else None
            if r.get("bg_dark_voxels"):
                r["bg_dark_fraction"] = r["bg_dark_voxels"] / 200 ** 3
                r["bg_external_share"] = r["bg_external_dark"] / r["bg_dark_voxels"]
    s = base / "e11" / "spatial" / "spatialization_summary.json"
    if s.exists():
        r["internal_components"] = json.loads(s.read_text(encoding="utf-8"))["internal_components"]
    diag = base / "e11" / "spatial" / "candidate_diagnostics_no_gt.csv"
    if diag.exists():
        d = pd.read_csv(diag, usecols=["polarity", "point_internal_component", "region_voxels"])
        m = d[d.polarity == "min"]
        r |= {"min_candidates": len(m), "min_in_internal": int((m.point_internal_component > 0).sum()),
              "min_region_voxels_median": float(m.region_voxels.median())}
    r["tos_markers"] = nrows(base / "e11" / "markers_tos_min" / "markers.csv")
    r["edt_otsu_markers"] = nrows(base / "e11" / "markers_edt_otsu" / "markers.csv")
    rs = base / "e23c" / "route" / "route_summary.csv"
    if rs.exists():
        g = pd.read_csv(rs).iloc[0]
        r |= {"gate_mv": g.region_voxels_median, "gate_qv": g.region_voxels_q75, "gate": bool(g.gate_triggered),
              "selected_source": g.selected_source}
    return r


def main() -> int:
    info = pd.read_csv(OUT / "crop_info.csv").set_index("crop_id")
    prim = pd.read_csv(OUT / "per_crop_primary.csv").set_index("crop_id")
    rows = []
    for crop in info.index:
        for method in ("frozen", "fixflat"):
            r = stages(method, crop)
            i = info.loc[crop]
            r = {"specimen": i.specimen, "test_class": i.test_class, "f": i.specimen_fraction_sub4,
                 "r_est_e13": i.r_est_e13, "openrec_radius": i.openrec_radius, "r_ps_peak": i.r_ps_peak,
                 "gt_r_in": i.gt_inscribed_r_median, "bg_ref": i.bg_ref_fixflat if method == "fixflat" else None,
                 **r, "f1_tos": prim.loc[crop, f"{method}_ToS_f1_050"], "f1_sel": prim.loc[crop, f"{method}_sel_f1_050"]}
            rows.append(r)
    df = pd.DataFrame(rows)
    df["r_ratio_e13_over_gt"] = df.r_est_e13 / df.gt_r_in
    df.to_csv(OUT / "forensics_stages.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    cols = ["crop_id", "method", "f", "r_est_e13", "openrec_radius", "bg_otsu", "bg_p40", "bg_dark_thr",
            "bg_thr_branch", "bg_dark_fraction", "bg_external_share", "bg_keep", "min_in_internal", "tos_markers",
            "gate_mv", "gate_qv", "crashed", "f1_tos", "f1_sel"]
    print(df[cols].round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
