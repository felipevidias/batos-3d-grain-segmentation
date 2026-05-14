#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import tifffile
from skimage.measure import label, regionprops

if len(sys.argv) < 2:
    print("Uso: python3 count_tiff_labels.py arquivo.tif [--binary]")
    sys.exit(1)

path = Path(sys.argv[1])
binary = "--binary" in sys.argv[2:]
arr = tifffile.imread(path)
if binary:
    lab = label(arr > 0, background=0, connectivity=1)
else:
    lab = arr
u = np.unique(lab)
u = u[u > 0]
print(f"file = {path}")
print(f"shape = {arr.shape}, dtype = {arr.dtype}, min = {arr.min()}, max = {arr.max()}, mean = {arr.mean():.6f}")
print(f"labels/components = {len(u)}")
if binary:
    areas = [int(r.area) for r in regionprops(lab)]
    if areas:
        print(f"area min/mean/max = {min(areas)} / {np.mean(areas):.2f} / {max(areas)}")
