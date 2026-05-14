#!/usr/bin/env python3
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import tifffile as tiff

ROOT = Path(".")
MANIFEST = ROOT / "results_batos/paper_tables/batos_local_split_large_labels/local_split_manifest.csv"
OUT_ROOT = ROOT / "results_batos/gmic_3d_view_batos_ls"
OUT_ROOT.mkdir(parents=True, exist_ok=True)


def check(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")
    return path


def read_label(path):
    arr = tiff.imread(str(path))
    if arr.ndim != 3:
        raise ValueError(f"Esperava label 3D escalar, recebi shape={arr.shape}: {path}")
    return arr.astype(np.int32)


def label_to_rgb_volume(labels):
    labels = labels.astype(np.int64)
    rgb = np.zeros(labels.shape + (3,), dtype=np.uint8)

    ids = np.unique(labels)
    ids = ids[ids > 0]

    for lab in ids:
        r = int((lab * 37) % 255)
        g = int((lab * 91) % 255)
        b = int((lab * 53) % 255)

        r = max(r, 60)
        g = max(g, 60)
        b = max(b, 60)

        rgb[labels == lab] = (r, g, b)

    return rgb


def prepare_one(row):
    scan = str(row["scan"])
    crop_id = str(row["crop_id"])
    name = f"{scan}_{crop_id}"

    gt_path = check(row["gt"])
    batos_path = check(row["original_pred"])
    ls_path = check(row["split_pred"])

    out = OUT_ROOT / name
    out.mkdir(parents=True, exist_ok=True)

    # 1) Cópias dos labels 3D escalares reais
    raw_gt = out / "01_GT_labels_3D_scalar.tif"
    raw_batos = out / "02_BA_TOS_labels_3D_scalar.tif"
    raw_ls = out / "03_BA_TOS_LS_labels_3D_scalar.tif"

    shutil.copy(gt_path, raw_gt)
    shutil.copy(batos_path, raw_batos)
    shutil.copy(ls_path, raw_ls)

    # 2) Versões coloridas 3D em RGB para visualizar melhor no G'MIC
    gt = read_label(gt_path)
    batos = read_label(batos_path)
    ls = read_label(ls_path)

    gt_rgb = label_to_rgb_volume(gt)
    batos_rgb = label_to_rgb_volume(batos)
    ls_rgb = label_to_rgb_volume(ls)

    rgb_gt = out / "01_GT_labels_3D_RGB.tif"
    rgb_batos = out / "02_BA_TOS_labels_3D_RGB.tif"
    rgb_ls = out / "03_BA_TOS_LS_labels_3D_RGB.tif"

    tiff.imwrite(str(rgb_gt), gt_rgb, photometric="rgb")
    tiff.imwrite(str(rgb_batos), batos_rgb, photometric="rgb")
    tiff.imwrite(str(rgb_ls), ls_rgb, photometric="rgb")

    # 3) Um volume 3D lado a lado: GT | BA-TOS | BA-TOS-LS
    z, y, x, c = gt_rgb.shape
    gap = np.ones((z, y, 8, c), dtype=np.uint8) * 255
    side = np.concatenate([gt_rgb, gap, batos_rgb, gap, ls_rgb], axis=2)

    side_out = out / "00_GT_BA_TOS_BA_TOS_LS_side_by_side_3D_RGB.tif"
    tiff.imwrite(str(side_out), side, photometric="rgb")

    print()
    print("====", name, "====")
    print("RAW GT       :", raw_gt)
    print("RAW BA-TOS   :", raw_batos)
    print("RAW BA-TOS-LS:", raw_ls)
    print("RGB GT       :", rgb_gt)
    print("RGB BA-TOS   :", rgb_batos)
    print("RGB BA-TOS-LS:", rgb_ls)
    print("SIDE-BY-SIDE :", side_out)


def main():
    if not MANIFEST.exists():
        raise FileNotFoundError(
            f"Manifesto não encontrado: {MANIFEST}\n"
            "Você provavelmente não está na pasta raiz do projeto.\n"
            "Rode antes:\n"
            "cd ~/Pesquisa-Grain_Seg/Filtro_Conexo/ToSConOpCpp/CodeReport_LysandreMacke"
        )

    df = pd.read_csv(MANIFEST)

    # Prepara todos os crops
    for _, row in df.iterrows():
        prepare_one(row)

    print()
    print("Tudo pronto em:", OUT_ROOT)


if __name__ == "__main__":
    main()
