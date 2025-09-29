"""Tests for scenario profile loading with optional PyYAML dependency handling."""
from types import SimpleNamespace

def test_scenario_profile_accessor(monkeypatch, tmp_path):
    # Create temp templates directory with scenario_profiles.yaml if yaml available
    try:
        import yaml  # type: ignore
    except Exception:
        yaml = None  # type: ignore
    from utils.config import AppConfig
    base = tmp_path
    (base / 'templates').mkdir()
    content = {
        'Electrostatic Focus': {
            'interactions': ['hydrogenbond'],
            'preset': 'default'
        }
    }
    if yaml is not None:
        (base / 'templates' / 'scenario_profiles.yaml').write_text(yaml.safe_dump(content))
    cfg = AppConfig(base_dir=base, presets={'default': {'hydrogenbond_distance_cutoff': 3.5}}, interactions=SimpleNamespace(hydrogenbond_distance_cutoff=3.5))
    prof = cfg.get_scenario_profile('electrostatic_focus')
    if yaml is None:
        assert prof is None
    else:
        assert prof is not None
        assert prof['interactions'] == ['hydrogenbond']
