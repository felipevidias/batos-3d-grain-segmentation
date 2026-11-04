import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from batos_v2 import run_crop
from batos_v2.common import load_config, synthetic_grains, write3d


class TestRunCrop(unittest.TestCase):
    def test_end_to_end_on_synthetic_crop(self):
        gray, labels = synthetic_grains()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            write3d(tmp / "crops" / "SYN_crop01" / "gray.tif", gray)
            write3d(tmp / "crops" / "SYN_crop01" / "labels.tif", labels)
            meta = {"crop_convention": {"slicing_status": "verified"},
                    "crops": [{"crop_id": "SYN_crop01", "scan": "SYN", "x": 0, "y": 0, "z": 0, "size": 64}]}
            (tmp / "meta.json").write_text(json.dumps(meta))
            cfg = load_config("configs/batos_paper_v2.yaml")
            cfg["paths"].update({"crops_metadata": str(tmp / "meta.json"), "crops_root": str(tmp / "crops"),
                            "runs_root": str(tmp / "runs"), "all_runs": str(tmp / "all_runs.csv")})
            cfg["eval"]["min_gt_voxels"] = 100
            (tmp / "cfg.yaml").write_text(yaml.safe_dump(cfg))
            self.assertEqual(run_crop.main(["--crop-id", "SYN_crop01", "--config", str(tmp / "cfg.yaml")]), 0)
            record = json.loads((tmp / "runs" / "SYN_crop01" / "run.json").read_text())
            for key in ("command", "timestamp", "git", "config_sha256", "inputs", "outputs", "metrics"):
                self.assertIn(key, record)
            metrics = pd.DataFrame(record["metrics"])
            self.assertEqual(set(metrics["method"]), {"batos_tos", "batos_sel", "edt_otsu"})
            self.assertEqual(set(metrics["size_filter"]), {"500", "off"})
            all_runs = pd.read_csv(tmp / "all_runs.csv")
            self.assertEqual(len(all_runs), len(metrics))
            # refuses to overwrite without --force
            with self.assertRaises(FileExistsError):
                run_crop.main(["--crop-id", "SYN_crop01", "--config", str(tmp / "cfg.yaml")])


if __name__ == "__main__":
    unittest.main()
