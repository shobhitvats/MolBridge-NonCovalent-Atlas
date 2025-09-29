"""Test configuration ensuring src package discoverability & settings reset helpers."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # points to src/
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

def reset_settings_cache():  # convenience for tests toggling env flags
    try:
        from utils.settings import get_settings
        get_settings.cache_clear()  # type: ignore
    except Exception:
        pass