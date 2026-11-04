import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from batos_v2 import gate
from batos_v2.common import load_config, read3d, synthetic_grains


def _frames(region_voxels):
    nodes = np.arange(100, 100 + len(region_voxels))
    diag = pd.DataFrame({"node": nodes, "polarity": "min", "score": 0.5,
                         "has_reconstructed_region": True, "region_voxels": region_voxels})
    return diag, pd.DataFrame({"node": nodes})


class TestGate(unittest.TestCase):
    def setUp(self):
        self.cfg = load_config("configs/batos_paper_v2.yaml")

    def test_collapsed_supports_trigger(self):
        diag, markers = _frames([1, 2, 2, 2, 3, 3, 3, 4, 40])  # median 3, q75 3
        self.assertTrue(gate.compute_gate(diag, markers, self.cfg["gate"])["gate_true"])

    def test_normal_supports_do_not_trigger(self):
        diag, markers = _frames([17, 20, 43, 60, 100])
        self.assertFalse(gate.compute_gate(diag, markers, self.cfg["gate"])["gate_true"])

    def test_q75_condition_matters(self):
        diag, markers = _frames([1, 2, 3, 3, 3, 10, 10, 10])  # median 3, q75 10
        g = gate.compute_gate(diag, markers, self.cfg["gate"])
        self.assertEqual(g["region_voxels_median"], 3.0)
        self.assertFalse(g["gate_true"])

    def test_missing_node_fails_loudly(self):
        diag, markers = _frames([5, 6])
        markers.loc[0, "node"] = 999
        with self.assertRaises(RuntimeError):
            gate.compute_gate(diag, markers, self.cfg["gate"])

    def test_edt_otsu_matched_count_and_routing(self):
        gray, labels = synthetic_grains()
        flat = (gray // 256).astype(np.uint8)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tos = np.zeros(gray.shape, np.uint16)
            pts = np.argwhere(labels > 0)[::997][:5]
            for i, p in enumerate(pts, 1):
                tos[tuple(p)] = i
            from batos_v2.common import write3d
            write3d(tmp / "tos.tif", tos)
            diag, markers = _frames([2] * 5)
            diag.to_csv(tmp / "diag.csv", index=False)
            markers.to_csv(tmp / "markers.csv", index=False)
            route = gate.run(tmp / "diag.csv", tmp / "markers.csv", tmp / "tos.tif", flat, self.cfg, tmp / "out")
            self.assertTrue(route["gate"]["gate_true"])
            self.assertEqual(route["selected_source"], "edt_otsu")
            self.assertEqual(route["marker_counts"]["edt_otsu"], 5)
            sel = read3d(route["selective_markers"])
            self.assertEqual(int((sel > 0).sum()), 5)


if __name__ == "__main__":
    unittest.main()
