from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    import tifffile
except ImportError:
    tifffile = None


INPUT_DIR = Path("results_batos/paper_figures/batos_ls_color_qc")
OUTPUT_DIR = Path("results_batos/paper_figures/batos_ls_stable_panels")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def normalize_gray(arr):
    arr = arr.astype(np.float32)
    p1, p99 = np.percentile(arr, [1, 99])
    if p99 <= p1:
        p1, p99 = arr.min(), arr.max()
    arr = (arr - p1) / (p99 - p1 + 1e-8)
    arr = np.clip(arr, 0, 1)
    return (arr * 255).astype(np.uint8)


def load_gray_mip(gray_path):
    if tifffile is None:
        raise RuntimeError("Instale tifffile: pip install tifffile")

    vol = tifffile.imread(str(gray_path))

    if vol.ndim == 3:
        mip = vol.max(axis=0)
    elif vol.ndim == 2:
        mip = vol
    else:
        raise ValueError(f"Formato inesperado para gray: {gray_path}, shape={vol.shape}")

    mip = normalize_gray(mip)
    return Image.fromarray(mip).convert("RGB")


def load_rgb_png(path):
    return Image.open(path).convert("RGB")


def resize_to(img, size):
    return img.resize(size, Image.Resampling.NEAREST)


def add_title(img, title):
    title_h = 32
    canvas = Image.new("RGB", (img.width, img.height + title_h), "white")
    canvas.paste(img, (0, title_h))

    draw = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 16)
    except Exception:
        font = ImageFont.load_default()

    bbox = draw.textbbox((0, 0), title, font=font)
    tw = bbox[2] - bbox[0]
    draw.text(((img.width - tw) // 2, 8), title, fill="black", font=font)
    return canvas


def make_panel(crop_dir):
    crop_id = crop_dir.name

    gray_path = crop_dir / f"{crop_id}_gray.tif"
    gt_path = crop_dir / f"{crop_id}_gt_overlay_MIP.png"
    batos_path = crop_dir / f"{crop_id}_batos_overlay_MIP.png"
    batos_ls_path = crop_dir / f"{crop_id}_batos_ls_overlay_MIP.png"

    required = [gray_path, gt_path, batos_path, batos_ls_path]
    missing = [p for p in required if not p.exists()]
    if missing:
        print(f"[SKIP] {crop_id}: faltando {missing}")
        return None

    gray = load_gray_mip(gray_path)
    gt = load_rgb_png(gt_path)
    batos = load_rgb_png(batos_path)
    batos_ls = load_rgb_png(batos_ls_path)

    target_size = gt.size
    gray = resize_to(gray, target_size)
    batos = resize_to(batos, target_size)
    batos_ls = resize_to(batos_ls, target_size)

    panels = [
        add_title(gray, "Grayscale MIP"),
        add_title(gt, "Reference GT"),
        add_title(batos, "BA-TOS"),
        add_title(batos_ls, "BA-TOS-LS"),
    ]

    gap = 12
    total_w = sum(p.width for p in panels) + gap * (len(panels) - 1)
    total_h = max(p.height for p in panels)

    out = Image.new("RGB", (total_w, total_h), "white")

    x = 0
    for p in panels:
        out.paste(p, (x, 0))
        x += p.width + gap

    out_path = OUTPUT_DIR / f"{crop_id}_stable_panel.png"
    out.save(out_path, dpi=(300, 300))
    print(f"[OK] {out_path}")
    return out_path


def main():
    crop_dirs = sorted([p for p in INPUT_DIR.iterdir() if p.is_dir()])

    if not crop_dirs:
        raise RuntimeError(f"Nenhum diretório de crop encontrado em {INPUT_DIR}")

    generated = []
    for crop_dir in crop_dirs:
        out = make_panel(crop_dir)
        if out is not None:
            generated.append(out)

    manifest = OUTPUT_DIR / "stable_panel_manifest.txt"
    with open(manifest, "w") as f:
        for p in generated:
            f.write(str(p) + "\n")

    print()
    print(f"Total de painéis gerados: {len(generated)}")
    print(f"Manifesto: {manifest}")


if __name__ == "__main__":
    main()
