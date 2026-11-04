import tempfile
import unittest
from pathlib import Path

import numpy as np

from batos_v2 import watershed
from batos_v2.common import load_config, read3d, synthetic_grains, write3d


class TestWatershed(unittest.TestCase):
    def test_size_filter_keeps_ids(self):
        labels = np.zeros((10, 10, 10), np.int32)
        labels[:5] = 3
        labels[9, 9, 9] = 7
        out = watershed.size_filter(labels, 2)
        self.assertEqual(set(np.unique(out)), {0, 3})

    def test_one_basin_per_grain_on_synthetic(self):
        gray, labels = synthetic_grains()
        flat = (gray // 256).astype(np.uint8)
        cfg = load_config("configs/batos_paper_v2.yaml")
        markers = np.zeros(gray.shape, np.uint16)
        for i, lab in enumerate(range(1, labels.max() + 1), 1):
            c = np.rint(np.argwhere(labels == lab).mean(axis=0)).astype(int)
            markers[tuple(c)] = i
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write3d(tmp / "tos.tif", markers)
            internal = np.zeros(gray.shape, np.uint8)
            info = watershed.run(flat, flat, internal, {"tos_e11_min": tmp / "tos.tif"}, cfg, tmp / "ws")
            rows = {r["size_filter"]: r for r in info["rows"]}
            self.assertEqual(set(rows), {"500", "off"})
            seg = read3d(rows["off"]["segmentation"])
            self.assertEqual(rows["off"]["prediction_count"], labels.max())
            # each GT grain is dominated by one basin
            for lab in range(1, labels.max() + 1):
                vals, counts = np.unique(seg[labels == lab], return_counts=True)
                self.assertGreater(counts.max() / counts.sum(), 0.9)
            # synthetic balls (~900 voxels) survive the 500-voxel filter
            self.assertEqual(rows["500"]["prediction_count"], labels.max())


if __name__ == "__main__":
    unittest.main()
