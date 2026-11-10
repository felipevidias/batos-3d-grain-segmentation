"""Export the exclusive ToS supports of the e10 spatialization as a label volume (Phase 5h, H3). GT is never read.

Runs the frozen e10 script's own functions (imported, not copied): exact 3D Tree of Shapes of the selector,
tree/CSV consistency check, propagation of the nearest selected ancestor. Voxel value = candidate_label (row order of
the true_dual candidate CSV, as in candidate_diagnostics_no_gt.csv); exclusive support of candidate c = {leaf == c}.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from batos_v2.frozen import SCRIPTS


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--image", type=Path, required=True)
    ap.add_argument("--candidate-csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    spec = importlib.util.spec_from_file_location("e10_frozen", SCRIPTS["spatialize"])
    e10 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(e10)
    image = e10.read3d(a.image, "Imagem")
    df = pd.read_csv(a.candidate_csv)
    e10.validate_candidate_csv(df)
    df["polarity"] = df["polarity"].astype(str).str.lower().str.strip()
    tree, alt = e10.hg.component_tree_tree_of_shapes_image3d(image)
    alt = np.asarray(alt, dtype=np.int64)
    ext = np.asarray(e10.hg.attribute_extrema(tree, alt), dtype=bool)
    checked = e10.add_tree_checks(df, tree, alt, ext)
    bad = int((~checked["polarity_matches_tree"]).sum() + (~checked["csv_alt_matches_tree"]).sum()
              + (~checked["csv_parent_alt_matches_tree"]).sum())
    if bad:
        raise RuntimeError(f"candidate CSV does not match the tree ({bad} mismatches)")
    leaf = e10.propagate_nearest_selected_ancestor(tree, checked["node"].astype(np.int64).to_numpy(), image.shape)
    del tree, alt, ext
    dtype = np.uint16 if len(df) < 65535 else np.uint32
    a.out.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(a.out, leaf.astype(dtype), compression="zlib")
    print(f"supports: {len(df)} candidates -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
