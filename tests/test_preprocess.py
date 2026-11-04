import tempfile
import unittest
from pathlib import Path

import numpy as np

from batos_v2 import preprocess
from batos_v2.common import load_config, read3d, synthetic_grains


class TestPreprocess(unittest.TestCase):
    def setUp(self):
        self.gray, self.labels = synthetic_grains()
        self.config = load_config("configs/batos_paper_v2.yaml")

    def test_scale_estimate_near_true_radius(self):
        scale = preprocess.estimate_scale(self.gray)
        self.assertFalse(scale["gt_used"])
        self.assertGreater(scale["estimated_grain_radius"], 3.0)
        self.assertLess(scale["estimated_grain_radius"], 9.0)

    def test_flatten_lambda_matches_legacy_formula(self):
        vol = (np.arange(8 * 9 * 10) % 251).reshape(8, 9, 10).astype(np.uint8)
        out = preprocess.flatten_volume_u8(vol, 1.0, 0.0)
        self.assertEqual(out.dtype, np.uint8)
        self.assertEqual(out.shape, vol.shape)

    def test_run_writes_e13_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            info = preprocess.run(self.gray, self.config, Path(tmp), "syn")
            openrec = read3d(info["openrec_image"])
            selector = read3d(info["selector_image"])
            self.assertTrue(Path(info["openrec_image"]).name.startswith("syn_flatten_openrec_r"))
            self.assertTrue(Path(info["selector_image"]).name.endswith("_flatten_lambda1p0_8u.tif"))
            self.assertEqual(openrec.dtype, np.uint8)
            self.assertEqual(selector.shape, self.gray.shape)
            # grains stay brighter than the matrix after flattening
            self.assertGreater(openrec[self.labels > 0].mean(), openrec[self.labels == 0].mean() + 50)


    def test_all_openrec_variants_run(self):
        variants = [{"mode": "tophat"}, {"mode": "divide_background"}, {"presmooth_sigma": 1.0},
                    {"input": "robust_u8"}, {"dimensionality": "2d_per_slice"}]
        for override in variants:
            cfg = dict(self.config["openrec"], **override)
            flat, stats = preprocess.openrec_flatten(self.gray, 8, cfg)
            self.assertEqual(flat.dtype, np.uint8, override)
            self.assertEqual(flat.shape, self.gray.shape, override)
        with self.assertRaises(ValueError):
            preprocess.openrec_flatten(self.gray, 8, dict(self.config["openrec"], mode="nope"))


if __name__ == "__main__":
    unittest.main()
