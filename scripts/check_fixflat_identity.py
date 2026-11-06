#!/usr/bin/env python3
"""batos-paper-v2-fixflat vs frozen: selector (lambda=1 flatten) and BA-TOS_ToS segmentation identity per crop.
Writes results/batos_v2_fixflat/identity_vs_frozen.csv."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

REPO = Path(__file__).resolve().parents[1]
F = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_frozen")
R = Path(r"C:\BA_TOS_DATA\runs\batos_paper_v2_fixflat\A")
SEG = "e11/watershed/segmentations/tos_e11_min_distance_c1.tif"


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def sel(b):
    return tifffile.imread(next(b.rglob("*_flatten_lambda1p0_8u.tif")))


rows = []
for s in (0, 1):
    for c in range(2, 8):
        k = f"EFRGP01_0{s}_crop0{c}"
        a, b = sel(F / k / "batos_openrec_auto_minparam" / "01_preprocess"), sel(R / k / "pipe" / "01_preprocess")
        d = np.abs(a.astype(np.int16) - b)
        sa, sb = tifffile.imread(F / k / SEG), tifffile.imread(R / k / SEG)
        rows.append({"crop_id": k, "selector_identical": sha(a) == sha(b), "selector_max_absdiff": int(d.max()),
                     "selector_frac_voxels_changed": float((d > 0).mean()),
                     "tos_segmentation_identical": sha(sa) == sha(sb),
                     "tos_segmentation_frac_voxels_changed": float((sa != sb).mean())})
df = pd.DataFrame(rows)
df.to_csv(REPO / "results" / "batos_v2_fixflat" / "identity_vs_frozen.csv", index=False)
print(df.to_string(index=False))
