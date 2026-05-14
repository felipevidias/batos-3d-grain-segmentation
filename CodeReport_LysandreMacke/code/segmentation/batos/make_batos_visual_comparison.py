#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt


def read3d(path):
    arr = np.asarray(tiff.imread(str(path)))
    if arr.ndim == 4 and arr.shape[-1] in (3, 4):
        arr = arr[..., 0]
    if arr.ndim != 3:
        raise RuntimeError(f"Expected 3D image, got {arr.shape} for {path}")
    return arr


def robust_u8(img):
    if img.dtype == np.uint8:
        return img
    v = img.astype(np.float32)
    lo, hi = np.percentile(v, [0.5, 99.5])
    if hi <= lo:
        hi = lo + 1
    return np.clip((v - lo) * 255.0 / (hi - lo), 0, 255).astype(np.uint8)


def colorize_labels_2d(lbl):
    lbl = lbl.astype(np.int64)
    out = np.zeros(lbl.shape + (3,), dtype=np.uint8)
    ids = np.unique(lbl)
    ids = ids[ids > 0]

    for lab in ids:
        r = (lab * 37 + 53) % 255
        g = (lab * 73 + 91) % 255
        b = (lab * 109 + 17) % 255
        out[lbl == lab] = (r, g, b)

    return out


def overlay(gray2d, lbl2d, alpha=0.55):
    gray = robust_u8(gray2d)
    rgb = np.repeat(gray[..., None], 3, axis=-1).astype(np.float32)
    col = colorize_labels_2d(lbl2d).astype(np.float32)
    mask = lbl2d > 0

    out = rgb.copy()
    out[mask] = (1 - alpha) * rgb[mask] + alpha * col[mask]
    return np.clip(out, 0, 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gray", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--otsu", required=True)
    ap.add_argument("--batos-raw", required=True)
    ap.add_argument("--batos-final", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="comparison")
    ap.add_argument("--z", type=int, default=-1)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gray = read3d(args.gray)
    gt = read3d(args.gt)
    otsu = read3d(args.otsu)
    batos_raw = read3d(args.batos_raw)
    batos_final = read3d(args.batos_final)

    if args.z < 0:
        z = gray.shape[0] // 2
    else:
        z = args.z

    panels = [
        ("Grayscale", robust_u8(gray[z])),
        ("Reference labels", colorize_labels_2d(gt[z])),
        ("Otsu + WS", overlay(gray[z], otsu[z])),
        ("BA-TOS raw", overlay(gray[z], batos_raw[z])),
        ("BA-TOS min500", overlay(gray[z], batos_final[z])),
    ]

    fig, axes = plt.subplots(1, len(panels), figsize=(18, 4))

    for ax, (title, img) in zip(axes, panels):
        if img.ndim == 2:
            ax.imshow(img, cmap="gray")
        else:
            ax.imshow(img)
        ax.set_title(title, fontsize=10)
        ax.axis("off")

    fig.suptitle(f"{args.prefix} | z={z}", fontsize=12)
    fig.tight_layout()

    png_out = out_dir / f"{args.prefix}_z{z}_comparison.png"
    pdf_out = out_dir / f"{args.prefix}_z{z}_comparison.pdf"

    fig.savefig(png_out, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_out, bbox_inches="tight")
    plt.close(fig)

    # Also save individual panels as PNG.
    for title, img in panels:
        safe = title.lower().replace(" ", "_").replace("+", "plus").replace("-", "_")
        p = out_dir / f"{args.prefix}_z{z}_{safe}.png"
        plt.figure(figsize=(4, 4))
        if img.ndim == 2:
            plt.imshow(img, cmap="gray")
        else:
            plt.imshow(img)
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(p, dpi=300, bbox_inches="tight", pad_inches=0)
        plt.close()

    print("wrote:", png_out)
    print("wrote:", pdf_out)
    print("out-dir:", out_dir)


if __name__ == "__main__":
    main()
