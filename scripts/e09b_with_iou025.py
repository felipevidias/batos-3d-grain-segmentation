#!/usr/bin/env python3
"""Run the frozen e09b_ignore_aware_audit.py unmodified, with IoU 0.25 added to its thresholds.

e09b computes IoU 0.50/0.75 only, via the default `thresholds=(0.50, 0.75)` of its own
compute_instance_metrics(); main() never passes thresholds. This wrapper imports the frozen module from
its path, sets that default to (0.25, 0.50, 0.75), and calls its main() with the given CLI. No other
change. The 0.50/0.75 columns must equal a plain e09b run (checked by aggregate_sibgrapi.py).

Usage: python -I scripts/e09b_with_iou025.py <path/to/e09b_ignore_aware_audit.py> <e09b args...>
"""
import importlib.util
import sys

path = sys.argv[1]
spec = importlib.util.spec_from_file_location("e09b_frozen", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.compute_instance_metrics.__defaults__ = ((0.25, 0.50, 0.75),)
sys.argv = [path, *sys.argv[2:]]
raise SystemExit(module.main())
