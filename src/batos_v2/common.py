"""Shared helpers for batos_v2: IO, config, hashing, robust 8-bit conversion."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
import yaml

CONNECTIVITY_6 = np.array(
    [[[0, 0, 0], [0, 1, 0], [0, 0, 0]],
     [[0, 1, 0], [1, 1, 1], [0, 1, 0]],
     [[0, 0, 0], [0, 1, 0], [0, 0, 0]]], dtype=bool)

REPO_ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.is_absolute():
        path = (Path.cwd() / path) if (Path.cwd() / path).exists() else REPO_ROOT / path
    with open(path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict) or config.get("preset") != "batos-paper-v2":
        raise ValueError(f"{path} is not a batos-paper-v2 config")
    return config


def config_hash(config: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(config, sort_keys=True).encode("utf-8")).hexdigest()


def read3d(path: str | Path, description: str = "volume") -> np.ndarray:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"{description}: {path}")
    array = np.squeeze(np.asarray(tifffile.imread(path)))
    if array.ndim == 4 and array.shape[-1] in (3, 4):
        array = array[..., 0]
    if array.ndim != 3:
        raise ValueError(f"{description} must be 3D, got {array.shape}: {path}")
    return array


def write3d(path: str | Path, array: np.ndarray, compress: bool = True) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(path, array, compression="zlib" if compress else None)
    return path


def array_sha256(array: np.ndarray) -> str:
    data = np.ascontiguousarray(array)
    data = data.astype(data.dtype.newbyteorder("<"), copy=False)
    return hashlib.sha256(data.tobytes()).hexdigest()


def robust_u8(volume: np.ndarray, p_low: float = 0.5, p_high: float = 99.5) -> np.ndarray:
    """Exact port of src/batos robust_u8: percentiles over >0 voxels; uint8 input is copied."""
    volume = np.asarray(volume)
    if volume.dtype == np.uint8:
        return volume.copy()
    values = volume.astype(np.float32)
    nonzero = values[values > 0]
    if nonzero.size == 0:
        nonzero = values.ravel()
    low = np.percentile(nonzero, p_low)
    high = np.percentile(nonzero, p_high)
    if high <= low:
        high = low + 1.0
    out = (values - low) * 255.0 / (high - low)
    return np.clip(out, 0, 255).astype(np.uint8)


def normalize_p1_p99(image: np.ndarray) -> np.ndarray:
    """Port of e11 normalize_image: p1-p99 over finite voxels -> [0, 1] float32."""
    work = image.astype(np.float32, copy=False)
    finite = np.isfinite(work)
    values = work[finite]
    low, high = np.percentile(values, [1.0, 99.0])
    if high <= low:
        low, high = float(values.min()), float(values.max())
    if high <= low:
        return np.zeros(work.shape, dtype=np.float32)
    result = np.clip((work - low) / (high - low), 0.0, 1.0)
    result[~finite] = 0.0
    return result.astype(np.float32, copy=False)


def write_json(path: str | Path, payload: dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=_json_default), encoding="utf-8")
    return path


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Not JSON serializable: {type(value)}")


def synthetic_grains(shape=(64, 64, 64), radius=6.0, spacing=15, seed=0,
                     bright=45000.0, dark=15000.0, noise=1500.0):
    """Synthetic test volume: bright balls on a dark matrix, each ball with a small dark core.

    Returns (uint16 gray, int32 labels). Used only by unit tests.
    """
    rng = np.random.default_rng(seed)
    zz, yy, xx = np.indices(shape)
    gray = np.full(shape, dark, dtype=np.float32)
    labels = np.zeros(shape, dtype=np.int32)
    label = 0
    offset = spacing // 2 + 2
    for z in range(offset, shape[0] - 4, spacing):
        for y in range(offset, shape[1] - 4, spacing):
            for x in range(offset, shape[2] - 4, spacing):
                cz, cy, cx = (np.array([z, y, x]) + rng.integers(-1, 2, 3)).tolist()
                d2 = (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2
                ball = d2 <= radius ** 2
                label += 1
                labels[ball] = label
                gray[ball] = bright
                gray[d2 <= 2] = bright - 12000.0  # dark core (internal minimum)
    gray += rng.normal(0, noise, shape).astype(np.float32)
    return np.clip(gray, 0, 65535).astype(np.uint16), labels
