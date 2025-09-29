"""Basic smoke tests for configuration loading and interaction registry helpers."""
from pathlib import Path
import sys

# Ensure src is on path
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.append(str(SRC))

from utils.config import load_config, get_interaction_types, get_interaction_display_names  # type: ignore


def test_load_config():
    cfg = load_config()
    assert cfg.app_name.lower().startswith("protein"), "Unexpected app name"
    assert cfg.interactions.hbond_distance_cutoff > 0


def test_interaction_type_lists():
    types = get_interaction_types()
    assert isinstance(types, list) and len(types) > 5
    display = get_interaction_display_names()
    for t in types:
        assert t in display, f"Missing display name for {t}"
