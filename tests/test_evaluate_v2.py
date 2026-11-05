import importlib.util
import unittest
from pathlib import Path

import numpy as np

from batos_v2 import evaluate_v2 as ev
from batos_v2.common import synthetic_grains

E09B = Path(r"C:\BA_TOS_DATA\linux_import\chat_recovered\e13_crop05_metadata_full_20260725_182408\audit"
            r"\e09b_ignore_aware_audit.py")
KW = dict(min_gt_voxels=100)


def seeds_for(labels, offset=0):
    """One-voxel marker at the first voxel of each label (marker label = label + offset)."""
    seeds = np.zeros(labels.shape, np.int32)
    for lab in np.unique(labels):
        if lab:
            seeds[tuple(np.argwhere(labels == lab)[len(np.argwhere(labels == lab)) // 2])] = lab + offset
    return seeds


def with_excluded(gt):
    """Add a face-touching GT grain (label 900) -> excluded; returns (gt, its mask)."""
    gt = gt.copy()
    gt[0:6, 0:6, 0:6] = 900
    return gt, gt == 900


class TestEvaluateV2(unittest.TestCase):
    def setUp(self):
        _, self.gt = synthetic_grains()

    def test_perfect_prediction(self):
        m, table, morph = ev.evaluate(self.gt, self.gt, seeds_for(self.gt), **KW)
        for s in ("025", "050", "075"):
            self.assertEqual(m[f"f1_iou{s}"], 1.0)
            self.assertEqual(m[f"ap_iou{s}"], 1.0)
            self.assertEqual(m[f"mean_matched_iou{s}"], 1.0)
        self.assertEqual(m["n_pred_ignored"], 0)
        self.assertEqual(len(morph), m["valid_gt_count"])
        for col in ("volume_rel_error", "d_eq_error", "centroid_error", "sphericity_error", "elongation_error"):
            self.assertTrue(np.allclose(morph[col], 0.0), col)

    def test_seeded_in_excluded_gt_is_ignored_with_and_without_markers(self):
        gt, ex = with_excluded(self.gt)
        pred = gt.copy()
        n_valid = len(np.unique(self.gt)) - 1
        with_m, table, _ = ev.evaluate(pred, gt, seeds_for(pred), **KW)
        self.assertEqual(with_m["n_pred_ignored"], 1)
        self.assertTrue(table.set_index("prediction_label").loc[900, "ignored"])
        self.assertEqual(with_m["fp_iou050"], 0)
        self.assertEqual(with_m["f1_iou050"], 1.0)
        no_m, table, _ = ev.evaluate(pred, gt, None, **KW)
        self.assertEqual(no_m["support"], "own_voxels")
        self.assertEqual(no_m["rule"], "fraction")
        self.assertEqual(no_m["n_pred_ignored"], 1)
        self.assertEqual(no_m["tp_iou050"], n_valid)
        self.assertEqual(no_m["fp_iou050"], 0)

    def test_background_seeded_prediction_stays_fp(self):
        pred = self.gt.copy()
        pred[30:36, 0:2, 30:36] = 777   # pure background blob (synthetic balls never reach y < 2)
        self.assertTrue((self.gt[pred == 777] == 0).all())
        for markers in (seeds_for(pred), None):
            m, table, _ = ev.evaluate(pred, self.gt, markers, **KW)
            self.assertEqual(table.set_index("prediction_label").loc[777, "support_category"], "background")
            self.assertFalse(table.set_index("prediction_label").loc[777, "ignored"])
            self.assertEqual(m["fp_iou050"], 1)
            self.assertEqual(m["n_pred_background_support"], 1)

    def test_label_invariance(self):
        """Prediction and marker labels unrelated to each other (Cellpose-like / relabelled) -> same result."""
        gt, _ = with_excluded(self.gt)
        pred = gt.copy()
        pred[30:36, 0:2, 30:36] = 777
        base, _, _ = ev.evaluate(pred, gt, seeds_for(pred), **KW)
        rng = np.random.default_rng(0)
        perm = np.concatenate([[0], rng.permutation(int(pred.max())) + 5000])
        permuted, _, _ = ev.evaluate(perm[pred], gt, seeds_for(pred, offset=40000), **KW)
        for key, value in base.items():
            if isinstance(value, float) and np.isnan(value):
                self.assertTrue(np.isnan(permuted[key]), key)
            else:
                self.assertAlmostEqual(value, permuted[key], places=12, msg=key)

    def test_dominant_vs_fraction_rule(self):
        gt, ex = with_excluded(self.gt)
        pred = gt.copy()
        markers = seeds_for(pred)
        markers[pred == 900] = 0
        # seed of prediction 900: 3 voxels in excluded GT, 2 in a valid grain -> dominant: ignored;
        # fraction 0.6 >= 0.5 ignored; >= 0.7 not ignored
        valid_lab = int(gt[9, 9, 9]) or int(np.unique(gt[gt > 0])[0])
        vox = np.argwhere(gt == valid_lab)[:2]
        pred[tuple(vox.T)] = 900
        markers[tuple(vox.T)] = 1
        markers[0, 0, 0:3] = 1
        prot = ev.reference_protocol(gt, 100)
        t = ev.classify_predictions(pred, prot, markers, "dominant").set_index("prediction_label")
        self.assertTrue(t.loc[900, "ignored"])
        self.assertAlmostEqual(t.loc[900, "support_excluded_fraction"], 0.6)
        t = ev.classify_predictions(pred, prot, markers, "fraction", 0.5).set_index("prediction_label")
        self.assertTrue(t.loc[900, "ignored"])
        t = ev.classify_predictions(pred, prot, markers, "fraction", 0.7).set_index("prediction_label")
        self.assertFalse(t.loc[900, "ignored"])

    def test_unseeded_prediction_is_evaluated(self):
        pred = self.gt.copy()
        markers = seeds_for(pred)
        lab = int(np.unique(pred[pred > 0])[0])
        markers[pred == lab] = 0
        m, table, _ = ev.evaluate(pred, self.gt, markers, **KW)
        self.assertEqual(table.set_index("prediction_label").loc[lab, "support_category"], "unseeded")
        self.assertEqual(m["n_pred_unseeded"], 1)
        self.assertEqual(m["f1_iou050"], 1.0)

    def test_matching_max_cardinality_then_max_iou(self):
        # two maximum matchings: {A1, B2} total 0.6 and {A2, B1} total 0.8 -> the latter
        pp, gg = np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1])
        iou = np.array([0.30, 0.40, 0.40, 0.30])
        rp, cg, v = ev._max_card_matching(pp, gg, iou, 2, 2, 0.25)
        self.assertEqual(sorted(zip(rp.tolist(), cg.tolist())), [(0, 1), (1, 0)])
        self.assertAlmostEqual(v.sum(), 0.8)
        # cardinality beats IoU: A-1 0.9 alone vs {A2 0.3, B1 0.3}
        pp, gg, iou = np.array([0, 0, 1]), np.array([0, 1, 0]), np.array([0.9, 0.3, 0.3])
        rp, cg, v = ev._max_card_matching(pp, gg, iou, 2, 2, 0.25)
        self.assertEqual(len(v), 2)

    def test_coco_ap(self):
        pp, gg, iou = np.array([1, 2]), np.array([0, 1]), np.array([0.9, 0.8])
        self.assertAlmostEqual(ev.coco_ap(pp, gg, iou, np.array([0.9, 0.5, 0.4]), 2, 0.5), 2 / 3, places=6)
        self.assertAlmostEqual(ev.coco_ap(pp, gg, iou, np.array([0.1, 0.9, 0.8]), 2, 0.5), 1.0)

    def test_scores_give_coco_ap_and_none_give_nan(self):
        m, _, _ = ev.evaluate(self.gt, self.gt, None, **KW)
        self.assertTrue(np.isnan(m["coco_ap_iou050"]))
        scores = {int(p): 1.0 for p in np.unique(self.gt) if p}
        m, _, _ = ev.evaluate(self.gt, self.gt, None, scores=scores, **KW)
        self.assertEqual(m["coco_ap_iou050"], 1.0)

    def test_morphology_descriptors(self):
        z, y, x = np.indices((40, 40, 60))
        ball = (z - 20) ** 2 + (y - 20) ** 2 + (x - 30) ** 2 <= 10 ** 2
        d = ev._descriptors(ball, (0, 0, 0), 1.0)
        self.assertAlmostEqual(d["d_eq"], 20.0, delta=0.5)
        self.assertGreater(d["sphericity"], 0.85)
        self.assertLessEqual(d["sphericity"], 1.0)
        self.assertAlmostEqual(d["elongation"], 1.0, delta=0.02)
        self.assertEqual((d["centroid_z"], d["centroid_y"], d["centroid_x"]), (20.0, 20.0, 30.0))
        ell = ((z - 20) / 6.0) ** 2 + ((y - 20) / 6.0) ** 2 + ((x - 30) / 18.0) ** 2 <= 1
        self.assertAlmostEqual(ev._descriptors(ell, (0, 0, 0), 1.0)["elongation"], 3.0, delta=0.15)
        self.assertAlmostEqual(ev._descriptors(ball, (0, 0, 0), 2.0)["volume"], ball.sum() * 8.0)

    def test_shape_mismatch_fails_loudly(self):
        with self.assertRaises(ValueError):
            ev.evaluate(self.gt[:-1], self.gt, None, **KW)

    @unittest.skipUnless(E09B.is_file(), "frozen e09b not available")
    def test_matches_frozen_e09b_seed_ignore(self):
        spec = importlib.util.spec_from_file_location("e09b", E09B)
        e09b = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(e09b)
        gt, _ = with_excluded(self.gt)
        pred = gt.copy()
        pred[30:36, 0:2, 30:36] = 777
        # watershed-like: prediction label == 6-connected marker component label in raster order
        markers = e09b.normalize_markers(seeds_for(pred), gt.shape)
        relabel = np.zeros(int(pred.max()) + 1, np.int64)
        for k in range(1, int(markers.max()) + 1):
            relabel[int(pred[markers == k][0])] = k
        pred = relabel[pred]
        prot = e09b.prepare_gt_protocol(gt.astype(np.int64), 100)
        mt = e09b.classify_marker_labels(markers, gt.astype(np.int64), prot)
        ignored = set(mt.loc[mt.marker_category == "ignored_gt", "marker_label"].astype(int))
        ref = e09b.compute_instance_metrics(pred.astype(np.int32), prot["valid_gt"], prot["ignore_mask"],
                                            ignored, (0.25, 0.50, 0.75))
        mine, _, _ = ev.evaluate(pred, gt, markers, **KW)
        for key in ("prediction_count_evaluable", "valid_gt_count", "tp_iou050", "fp_iou050", "fn_iou050",
                    "f1_iou025", "f1_iou050", "f1_iou075", "mean_matched_iou050", "mean_matched_iou075"):
            self.assertAlmostEqual(mine[key], ref[key], places=12, msg=key)


if __name__ == "__main__":
    unittest.main()
