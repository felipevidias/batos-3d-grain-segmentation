#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile as tiff
import json
import math

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/batos_local_split_large_labels/local_split_manifest.csv"
OUT = ROOT / "results_batos/paper_tables/batos_ls_auto_params"
OUT.mkdir(parents=True, exist_ok=True)


def read_labels(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise ValueError(f"Expected 3D label volume, got {arr.shape}: {path}")
    return arr.astype(np.int32, copy=False)


def label_sizes(labels):
    ids, counts = np.unique(labels, return_counts=True)
    sizes = np.array([c for i, c in zip(ids, counts) if i > 0], dtype=np.float64)
    return sizes


def robust_grain_scale_from_sizes(sizes):
    if len(sizes) == 0:
        raise RuntimeError("No labels found.")

    # Remove extremos fortes: fragmentos pequenos e labels fundidos muito grandes.
    q10, q90 = np.percentile(sizes, [10, 90])
    clean = sizes[(sizes >= q10) & (sizes <= q90)]

    if len(clean) < 5:
        clean = sizes

    v_med = float(np.median(clean))
    v_q25 = float(np.percentile(clean, 25))
    v_q75 = float(np.percentile(clean, 75))

    r_eq = float((3.0 * v_med / (4.0 * math.pi)) ** (1.0 / 3.0))

    params = {
        "grain_volume_median": v_med,
        "grain_volume_q25": v_q25,
        "grain_volume_q75": v_q75,
        "equivalent_radius": r_eq,

        # Parâmetros automáticos do local split
        "MERGE_MIN_VOLUME": int(round(1.8 * v_med)),
        "MIN_SPLIT_PART_SIZE": int(round(0.25 * v_med)),
        "PEAK_MIN_DISTANCE": int(max(3, round(0.45 * r_eq))),
        "PEAK_THRESHOLD": float(max(2.0, 0.25 * r_eq)),

        # Mantém limite de segurança, mas o split pode usar expected_peaks por label.
        "MAX_PEAKS_GLOBAL_CAP": 6,
    }

    return params


def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(
            f"Manifesto não encontrado: {MANIFEST}\n"
            "Rode a partir da raiz do projeto."
        )

    df = pd.read_csv(MANIFEST)

    all_sizes = []
    per_crop_rows = []

    for _, r in df.iterrows():
        scan = str(r["scan"])
        crop_id = str(r["crop_id"])
        pred_path = Path(str(r["original_pred"]))

        labels = read_labels(pred_path)
        sizes = label_sizes(labels)

        if len(sizes) == 0:
            continue

        all_sizes.extend(sizes.tolist())

        p = robust_grain_scale_from_sizes(sizes)

        row = {
            "scan": scan,
            "crop_id": crop_id,
            "n_labels": int(len(sizes)),
            **p,
        }
        per_crop_rows.append(row)

    per_crop = pd.DataFrame(per_crop_rows)
    per_crop.to_csv(OUT / "auto_params_per_crop.csv", index=False)

    global_params = robust_grain_scale_from_sizes(np.array(all_sizes, dtype=np.float64))

    with open(OUT / "auto_params_global.json", "w") as f:
        json.dump(global_params, f, indent=2)

    print("=== AUTO PARAMS GLOBAL ===")
    for k, v in global_params.items():
        print(f"{k}: {v}")

    print()
    print("wrote:", OUT / "auto_params_per_crop.csv")
    print("wrote:", OUT / "auto_params_global.json")


if __name__ == "__main__":
    main()
