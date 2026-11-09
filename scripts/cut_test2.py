#!/usr/bin/env python3
"""Phase 5f step 4: cut test set 2 with the Phase-1f cutter (scripts/extract_and_cut.py, unchanged) re-pointed at
data/metadata/crops_test2.json -> C:\\BA_TOS_DATA\\crops_test2, manifest data/metadata/crops_manifest_test2.csv.
Arguments are passed through (--specimens, --keep-extracted).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("extract_and_cut", REPO / "scripts" / "extract_and_cut.py")
eac = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(eac)

eac.META = REPO / "data" / "metadata" / "crops_test2.json"
eac.MANIFEST = REPO / "data" / "metadata" / "crops_manifest_test2.csv"
eac.OUT = Path(r"C:\BA_TOS_DATA\crops_test2")

if __name__ == "__main__":
    raise SystemExit(eac.main())
