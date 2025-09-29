"""Tests for command palette core generation and execution safety.
These are lightweight and avoid Streamlit heavy runtime by isolating helper logic.
"""
from types import SimpleNamespace

def build_interface(monkeypatch):
    from ui.main_interface import MainInterface
    from utils.config import load_config
    cfg = load_config()
    # Inject a minimal preset & profile for testing
    cfg.presets['default'] = {'hydrogenbond_distance_cutoff': getattr(cfg.interactions, 'hydrogenbond_distance_cutoff', 3.5)}
    cfg.scenario_profiles = {'electrostatic_focus': {'interactions': ['hydrogenbond'], 'preset': 'default'}}
    mi = MainInterface(cfg)
    return mi

def test_command_items_groups(monkeypatch):
    mi = build_interface(monkeypatch)
    items = mi._command_items()
    labels = [it['label'] for it in items]
    assert any(l.startswith('Apply Preset:') for l in labels)
    assert any(l.startswith('Apply Profile:') for l in labels)
    assert any('Save Layout Snapshot' in l for l in labels)


def test_execute_preset(monkeypatch):
    mi = build_interface(monkeypatch)
    # Change value then apply preset
    mi.config.interactions.hydrogenbond_distance_cutoff = 4.0
    mi._execute_command_palette_action('preset:default')
    assert mi.config.interactions.hydrogenbond_distance_cutoff == 3.5


def test_execute_profile(monkeypatch):
    import streamlit as st
    mi = build_interface(monkeypatch)
    st.session_state.selected_interactions = {'hydrogenbond': False, 'ionicinteraction': False}
    mi._execute_command_palette_action('profile:electrostatic_focus')
    assert st.session_state.selected_interactions['hydrogenbond'] is True


def test_layout_snapshot_version_field(monkeypatch, tmp_path):
    import json, streamlit as st
    mi = build_interface(monkeypatch)
    st.session_state.selected_interactions = {'hydrogenbond': True}
    st.session_state.individual_strength_filters = {'hydrogenbond': 'all'}
    # simulate save
    d = tmp_path / 'cache'
    d.mkdir()
    import time
    snap_path = d / 'layout_123.json'
    snap_path.write_text(json.dumps({'version': mi.LAYOUT_SNAPSHOT_VERSION, 'selected_interactions': {'hydrogenbond': False}, 'individual_strength_filters': {'hydrogenbond': 'all'}}))
    # monkeypatch Path('cache') to point to tmp
    from pathlib import Path
    monkeypatch.setattr('ui.main_interface.Path', lambda p='cache': Path(d))
    mi._execute_command_palette_action('layout:restore_last')
    # interactions should update from snapshot
    assert st.session_state.selected_interactions['hydrogenbond'] is False
