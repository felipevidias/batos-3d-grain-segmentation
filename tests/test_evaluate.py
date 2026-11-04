import unittest

import numpy as np
import pandas as pd

from batos_v2 import evaluate
from batos_v2.common import load_config, synthetic_grains

CFG = load_config("configs/batos_paper_v2.yaml")["eval"]
CFG_SMALL = dict(CFG, min_gt_voxels=100)


def _seeds_for(labels):
    seeds = np.zeros(labels.shape, np.int32)
    for lab in np.unique(labels):
        if lab:
            c = np.argwhere(labels == lab)[0]
            seeds[tuple(c)] = lab
    return seeds


class TestEvaluate(unittest.TestCase):
    def setUp(self):
        _, self.gt = synthetic_grains()

    def test_perfect_prediction(self):
        rows = evaluate.evaluate_crop(self.gt, _seeds_for(self.gt), self.gt, CFG_SMALL)
        seed = rows[0]
        self.assertEqual(seed["evaluation_protocol"], "seed_ignore")
        for s in ("025", "050", "075"):
            self.assertEqual(seed[f"f1_iou{s}"], 1.0)
        self.assertEqual(seed["n_pred_ignored"], 0)

    def test_prediction_seeded_in_excluded_gt_is_ignored(self):
        gt = self.gt.copy()
        gt[0, 0, 0] = 999                      # tiny face-touching label -> excluded
        gt[0:3, 0:3, 0:3] = 999
        pred = gt.copy()
        rows = evaluate.evaluate_crop(pred, _seeds_for(pred), gt, CFG_SMALL)
        seed, strict = rows
        self.assertEqual(seed["n_pred_ignored"], 1)
        self.assertEqual(seed["fp_iou050"], 0)
        self.assertEqual(seed["f1_iou050"], 1.0)
        self.assertEqual(strict["fp_iou050"], 0)  # fully inside excluded -> not in strict domain

    def test_background_seeded_prediction_is_fp(self):
        pred = self.gt.copy()
        pred[pred == 0] = 0
        pred[30:33, 0:3, 30:33] = 777          # in background (gt==0)
        self.assertTrue((self.gt[30:33, 0:3, 30:33] == 0).all())
        rows = evaluate.evaluate_crop(pred, _seeds_for(pred), self.gt, CFG_SMALL)
        self.assertEqual(rows[0]["fp_iou050"], 1)
        self.assertEqual(rows[0]["n_pred_ignored"], 0)

    def test_merge_gives_fn(self):
        pred = self.gt.copy()
        pred[pred == 2] = 1                    # merge two grains
        rows = evaluate.evaluate_crop(pred, _seeds_for(pred), self.gt, CFG_SMALL)
        n = int(self.gt.max())
        # the merged blob has IoU ~0.5 with each grain: it may match one of them at 0.50, none at 0.75
        self.assertIn(rows[0]["tp_iou050"], (n - 2, n - 1))
        self.assertEqual(rows[0]["tp_iou075"], n - 2)
        self.assertEqual(rows[0]["fn_iou075"], 2)

    def test_summaries_are_labelled(self):
        rows = evaluate.evaluate_crop(self.gt, _seeds_for(self.gt), self.gt, CFG_SMALL)
        frame = pd.DataFrame([dict(crop_id=c, method="m", size_filter="off", **r) for c in ("a", "b") for r in rows])
        pooled, macro = evaluate.summarize(frame, CFG["iou_thresholds"])
        self.assertTrue((pooled["summary_type"] == "pooled_object").all())
        self.assertTrue((macro["summary_type"] == "crop_macro").all())
        self.assertEqual(float(macro.loc[macro.evaluation_protocol == "seed_ignore", "f1_iou050_mean"].iloc[0]), 1.0)


if __name__ == "__main__":
    unittest.main()
