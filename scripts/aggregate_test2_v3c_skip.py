#!/usr/bin/env python3
"""Run the frozen scripts/aggregate_test2_v3c.py unchanged, with N_CROPS = number of crops in data/metadata/crops_test2.json
(test set 2 was defined by scripts/define_test_crops2_skip.py, which can skip empty stratum cells; the frozen aggregator
hard-codes 52). Nothing else is changed.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("agg_t2", REPO / "scripts" / "aggregate_test2_v3c.py")
agg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(agg)

if __name__ == "__main__":
    agg.N_CROPS = len(json.loads((REPO / "data" / "metadata" / "crops_test2.json").read_text(encoding="utf-8"))["crops"])
    print(f"[agg] N_CROPS = {agg.N_CROPS}", file=sys.stderr)
    raise SystemExit(agg.main())
