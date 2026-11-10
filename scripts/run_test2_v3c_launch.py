#!/usr/bin/env python3
"""Launch the frozen scripts/run_test2_v3c.py unchanged, with one correction to its v3 freeze check.

configs/batos_paper_v3_scalefree.yaml records a MIXED hash convention: of its 25 entries, 13 match the raw bytes of the
current CRLF working copy and 12 match the LF-normalised bytes (each file matches exactly one; checked 2026-11-10). The
frozen run_test2_v3c.check_cfg accepts only LF-normalised hashes and so refused to start (nothing was run). This
launcher accepts a recorded hash if it equals EITHER the raw or the LF-normalised SHA-256 of the file, for both configs
(for batos_paper_v3c.yaml, recorded LF-normalised, this is the same test), and keeps the `git diff <tag>` check. Then it
calls the frozen main().
"""

from __future__ import annotations

import hashlib
import importlib.util
import subprocess
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("rt2", REPO / "scripts" / "run_test2_v3c.py")
rt2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rt2)


def _ok(p: Path, v: str) -> bool:
    b = p.read_bytes()
    return v in (hashlib.sha256(b).hexdigest(), hashlib.sha256(b.replace(b"\r\n", b"\n")).hexdigest())


def check_cfg(cfg_path: Path, tag: str | None) -> dict:
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if cfg.get("status") != "frozen_for_test":
        raise RuntimeError(f"{cfg_path} is not frozen_for_test")
    bad = [k for k, v in cfg["code_sha256_at_freeze"].items() if not _ok(REPO / k, v)]
    bad += [k for k, v in cfg["recovered_scripts_sha256"].items() if not _ok(Path(v["path"]), v["sha256"])]
    if tag:
        d = subprocess.run(["git", "diff", "--name-only", tag, "--", *cfg["code_sha256_at_freeze"]], cwd=REPO,
                           capture_output=True, text=True, check=True).stdout.split()
        bad += [f"git diff vs {tag}: {x}" for x in d]
    if bad:
        raise RuntimeError(f"code changed since the freeze recorded in {cfg_path.name}: {bad}")
    rt2.log(f"[launch] {cfg_path.name}: {len(cfg['code_sha256_at_freeze'])} + "
            f"{len(cfg['recovered_scripts_sha256'])} hashes OK (raw or LF-normalised)")
    return cfg


if __name__ == "__main__":
    rt2.check_cfg = check_cfg
    raise SystemExit(rt2.main())
