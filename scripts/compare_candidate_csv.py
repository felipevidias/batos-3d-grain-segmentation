#!/usr/bin/env python3
"""Compare a candidate_diagnostics_no_gt.csv (or *_tosconop_autodual_candidates.csv) with a reference.

Reports: point-coordinate overlap (crop identification), exact row equality on (node, polarity, point),
per-column agreement for shared numeric columns (max |diff|, fraction equal), and the first rank
where the two lists diverge. Read-only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def compare(ours: Path, ref: Path) -> dict:
    a, b = pd.read_csv(ours), pd.read_csv(ref)
    out: dict = {"ours": str(ours), "reference": str(ref), "rows_ours": len(a), "rows_reference": len(b)}
    if {"point_z", "point_y", "point_x"} <= set(a.columns) & set(b.columns):
        pa = set(map(tuple, a[["point_z", "point_y", "point_x"]].dropna().astype(int).to_numpy()))
        pb = set(map(tuple, b[["point_z", "point_y", "point_x"]].dropna().astype(int).to_numpy()))
        out["point_overlap_fraction_of_reference"] = len(pa & pb) / max(1, len(pb))
    key = [c for c in ("node", "polarity", "point_z", "point_y", "point_x") if c in a.columns and c in b.columns]
    n = min(len(a), len(b))
    same_rows = (a[key].iloc[:n].astype(str).to_numpy() == b[key].iloc[:n].astype(str).to_numpy()).all(axis=1)
    out["rowwise_key_columns"] = key
    out["rowwise_identical_fraction"] = float(same_rows.mean()) if n else 0.0
    out["first_divergent_rank"] = int(np.argmin(same_rows)) + 1 if not same_rows.all() else None
    out["node_set_overlap_fraction"] = (len(set(a["node"]) & set(b["node"])) / max(1, len(set(b["node"])))
                                        if "node" in key else None)
    columns = {}
    for col in [c for c in b.columns if c in a.columns]:
        x, y = a[col].iloc[:n], b[col].iloc[:n]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            d = np.abs(x.to_numpy(float) - y.to_numpy(float))
            columns[col] = {"max_abs_diff": float(np.nanmax(d)) if n else None,
                            "fraction_equal": float(np.mean(np.isclose(x, y, equal_nan=True)))}
        else:
            columns[col] = {"fraction_equal": float((x.astype(str).to_numpy() == y.astype(str).to_numpy()).mean())}
    out["columns"] = columns
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ours", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    result = compare(args.ours, args.reference)
    text = json.dumps(result, indent=1)
    if args.out:
        args.out.write_text(text, encoding="utf-8")
    summary = {k: v for k, v in result.items() if k != "columns"}
    print(json.dumps(summary, indent=1))
    bad = {c: v for c, v in result["columns"].items() if v.get("fraction_equal", 1) < 1}
    print(f"columns compared: {len(result['columns'])}; not fully equal: {len(bad)}")
    for c, v in list(bad.items())[:20]:
        print(f"  {c}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
