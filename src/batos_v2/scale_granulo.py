"""Phase 5b side study: GT-free granulometric grain-radius estimators (compared with e13's r_est).

Foreground = e13's own Otsu foreground (robust p1-p99 uint8 -> Otsu >= -> closing ball(1) -> fill holes), imported
read-only from the frozen e13_estimate_scale.py so that it is identical.

1. Pattern spectrum (binary granulometry). Opening by a ball of radius r (skimage ball: |d|^2 <= r^2), computed
   exactly with Euclidean distance transforms: erosion = {EDT_fg > r}, dilation of that set = {EDT to it <= r}.
   The crop outside is NOT treated as background (grains cut by a face are not shrunk by the face). V(r) = opened
   volume; PS(r) = V(r) - V(r+1), r = 2..30. Estimates: r_ps_peak = argmax PS; r_ps_median = radius at which half
   of the volume removed between r = 2 and r = 31 has been removed.
2. EDT local maxima: median EDT value over the voxels that are local maxima of the EDT in their 3x3x3 neighbourhood
   inside the foreground (no percentile cut, no top-N truncation; e13 uses a 5^3 window, keeps peaks >= EDT p75 and
   the 1500 largest).
GT is never read here.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import scipy.ndimage as ndi
import tifffile

E13 = Path(r"C:\BA_TOS_DATA\linux_import\linux_patches\BA_TOS_E13_E16_conditional_acceptance\audit\e13_estimate_scale.py")
RADII = range(2, 31)


def load_e13():
    spec = importlib.util.spec_from_file_location("e13_estimate_scale_ro", E13)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def opened_volume(edt_fg: np.ndarray, r: int) -> int:
    core = edt_fg > r
    if not core.any():
        return 0
    return int((ndi.distance_transform_edt(~core) <= r).sum())


def estimate(gray: np.ndarray) -> dict:
    e13 = load_e13()
    fg, thr = e13.make_foreground(e13.robust_u8(gray))
    edt = ndi.distance_transform_edt(fg)
    vols = {r: opened_volume(edt, r) for r in range(2, 32)}
    ps = {r: vols[r] - vols[r + 1] for r in RADII}
    radii = np.array(list(ps))
    spec = np.array([ps[r] for r in radii], dtype=float)
    cum = np.cumsum(spec) / max(1.0, spec.sum())
    r_peak = int(radii[int(np.argmax(spec))])
    r_med = float(radii[int(np.searchsorted(cum, 0.5))])
    loc = (edt == ndi.maximum_filter(edt, size=3)) & fg & (edt > 0)
    r_edt = float(np.median(edt[loc]))
    e13_r, method, npk = e13.estimate_radius(fg)
    return {"otsu_threshold_u8": thr, "foreground_fraction": float(fg.mean()),
            "r_ps_peak": r_peak, "r_ps_median": r_med, "r_edt_localmax_median": r_edt,
            "n_edt_localmax": int(loc.sum()), "e13_r_est_recomputed": float(e13_r), "e13_method": method,
            "pattern_spectrum": {int(r): int(v) for r, v in ps.items()},
            "opened_volume": {int(r): int(v) for r, v in vols.items()}}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="GT-free granulometric grain radius")
    ap.add_argument("gray", type=Path)
    ap.add_argument("--out-json", type=Path, required=True)
    args = ap.parse_args(argv)
    res = estimate(np.squeeze(tifffile.imread(args.gray)))
    res["gray"] = str(args.gray)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k not in ("pattern_spectrum", "opened_volume")}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
