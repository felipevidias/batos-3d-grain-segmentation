#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff


def read_crop_tiff(path, z0, y0, x0, size):
    """
    Lê apenas o crop necessário.

    Primeiro tenta memmap. Se o TIFF não permitir memmap,
    cai para leitura por páginas usando key=range(z0,z1).
    """
    path = Path(path)
    z1 = z0 + size
    y1 = y0 + size
    x1 = x0 + size

    try:
        arr = tiff.memmap(str(path))
        print(f"[crop] memmap OK: {path}")
        print(f"[crop] shape={arr.shape}, dtype={arr.dtype}")

        if arr.ndim != 3:
            raise RuntimeError(f"Esperava volume 3D, shape={arr.shape}")

        crop = np.asarray(arr[z0:z1, y0:y1, x0:x1])
        return crop

    except Exception as e:
        print(f"[crop] memmap falhou para {path}: {e}")
        print("[crop] tentando leitura por páginas...")

        pages = tiff.imread(str(path), key=range(z0, z1))
        pages = np.asarray(pages)

        if pages.ndim != 3:
            raise RuntimeError(f"Leitura por páginas não retornou 3D, shape={pages.shape}")

        print(f"[crop] pages shape={pages.shape}, dtype={pages.dtype}")
        crop = pages[:, y0:y1, x0:x1]
        return crop


def robust_u8(vol):
    if vol.dtype == np.uint8:
        return vol.copy()

    v = vol.astype(np.float32)
    nz = v[v > 0]
    if nz.size == 0:
        nz = v.ravel()

    lo = np.percentile(nz, 0.5)
    hi = np.percentile(nz, 99.5)

    if hi <= lo:
        hi = lo + 1.0

    out = (v - lo) * 255.0 / (hi - lo)
    return np.clip(out, 0, 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gray", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--x", type=int, required=True)
    ap.add_argument("--y", type=int, required=True)
    ap.add_argument("--z", type=int, required=True)
    ap.add_argument("--size", type=int, default=200)
    ap.add_argument("--prefix", default="paired_crop")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[crop] gray :", args.gray)
    print("[crop] label:", args.label)
    print(f"[crop] origin: x={args.x}, y={args.y}, z={args.z}, size={args.size}")

    gray_crop = read_crop_tiff(args.gray, args.z, args.y, args.x, args.size)
    label_crop = read_crop_tiff(args.label, args.z, args.y, args.x, args.size)

    print("[crop] gray crop shape :", gray_crop.shape, gray_crop.dtype)
    print("[crop] label crop shape:", label_crop.shape, label_crop.dtype)

    if gray_crop.shape != label_crop.shape:
        raise RuntimeError(f"Shape mismatch: gray={gray_crop.shape}, label={label_crop.shape}")

    gray_crop_8u = robust_u8(gray_crop)

    gray_out = out_dir / f"{args.prefix}_gray_crop.tif"
    gray8_out = out_dir / f"{args.prefix}_gray_crop_8u.tif"
    label_out = out_dir / f"{args.prefix}_label_crop.tif"

    tiff.imwrite(gray_out, gray_crop)
    tiff.imwrite(gray8_out, gray_crop_8u)
    tiff.imwrite(label_out, label_crop)

    labels = np.unique(label_crop)
    labels_fg = labels[labels > 0]

    print("[crop] wrote gray     :", gray_out)
    print("[crop] wrote gray 8u  :", gray8_out)
    print("[crop] wrote label    :", label_out)
    print("[crop] gray min/max   :", int(gray_crop.min()), int(gray_crop.max()))
    print("[crop] label min/max  :", int(label_crop.min()), int(label_crop.max()))
    print("[crop] foreground IDs :", len(labels_fg))

    if len(labels_fg) > 0:
        print("[crop] first IDs      :", labels_fg[:20])


if __name__ == "__main__":
    main()
