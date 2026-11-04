"""Run recovered linux_patches scripts verbatim (unmodified) as pipeline stages.

The E13 README lists the frozen E11 settings, and the E11 README states that
e11_build_markers_ablation.py in MIN mode reproduces the E10 markers voxel by voxel.
So these recovered scripts ARE the frozen marker code. They are called as subprocesses from
paths.e11_audit_dir; nothing is copied or edited.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def _script(config: dict[str, Any], name: str) -> Path:
    path = Path(config["paths"]["e11_audit_dir"]) / name
    if not path.is_file():
        raise FileNotFoundError(f"Recovered script not found: {path}")
    return path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _call(script: Path, args: list[str]) -> str:
    command = [sys.executable, "-I", str(script), *args]
    done = subprocess.run(command, capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(f"{script.name} failed ({done.returncode}):\n{done.stderr[-2000:]}")
    return done.stdout


def e11_markers(config: dict[str, Any], diagnostics_csv: Path, internal_tif: Path, out_dir: Path) -> dict[str, Any]:
    """e11_build_markers_ablation.py --polarity/--nms-radius/--edt-priority-weight/--min-priority-bonus."""
    cfg = config["markers"]
    script = _script(config, "e11_build_markers_ablation.py")
    _call(script, ["--diagnostic", str(diagnostics_csv), "--internal-source", str(internal_tif),
                   "--polarity", cfg["polarity"], "--nms-radius", str(cfg["nms_radius"]),
                   "--edt-priority-weight", str(cfg["edt_priority_weight"]),
                   "--min-priority-bonus", str(cfg["min_priority_bonus"]), "--out-dir", str(out_dir)])
    summary = json.loads((Path(out_dir) / "marker_config.json").read_text(encoding="utf-8"))
    return {"implementation": "e11_verbatim", "script": str(script), "script_sha256": _sha256(script),
            "markers_csv": str(Path(out_dir) / "markers.csv"), "markers_tif": str(Path(out_dir) / "markers_uint16.tif"),
            "marker_count": summary["marker_count"], "markers_before_nms": summary["markers_before_nms"],
            "candidates_inside_internal": summary["candidates_inside_internal"],
            "input_spatial_candidates": summary["input_spatial_candidates"]}


def e11_edt_otsu(config: dict[str, Any], flat_image_tif: Path, match_count_from: Path, out_dir: Path) -> dict[str, Any]:
    """e11_build_control_markers.py --mode otsu_edt --match-count-from <ToS markers>."""
    script = _script(config, "e11_build_control_markers.py")
    _call(script, ["--image", str(flat_image_tif), "--mode", "otsu_edt",
                   "--nms-radius", str(config["gate"]["edt_nms_radius"]),
                   "--match-count-from", str(match_count_from), "--out-dir", str(out_dir)])
    summary = json.loads((Path(out_dir) / "control_config.json").read_text(encoding="utf-8"))
    return {"implementation": "e11_verbatim", "script": str(script), "script_sha256": _sha256(script),
            "markers_tif": str(Path(out_dir) / "markers_uint16.tif"), **summary}
