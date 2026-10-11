"""Where data and run outputs live.

Everything bulky goes under ASTRO_DATA (default: ./data in the repository),
which is git-ignored. Set ASTRO_DATA to an external drive on the MacBook if
needed.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def data_root() -> Path:
    root = Path(os.environ.get("ASTRO_DATA", REPO_ROOT / "data")).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    return root


def sub(*parts: str) -> Path:
    """Return (and create) a directory under the data root."""
    p = data_root().joinpath(*parts)
    p.mkdir(parents=True, exist_ok=True)
    return p


def config_dir() -> Path:
    return REPO_ROOT / "config"
