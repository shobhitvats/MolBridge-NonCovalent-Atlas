"""Preset Manager Component

Provides a lightweight interface to list, apply, rename and delete saved
parameter presets stored in Streamlit session state. Each preset entry:
  name -> { values: {param: float}, created_at: iso, diff: optional }

This file intentionally keeps zero business logic about parameter semantics;
the main interface owns the authoritative parameter application.
"""
from __future__ import annotations
from typing import Dict, Any
import streamlit as st


def _compute_diff(current: Dict[str, float], baseline: Dict[str, float]):
    diff = {}
    for k, v in current.items():
        b = baseline.get(k)
        if b is not None and abs(float(v) - float(b)) > 1e-9:
            diff[k] = {'from': b, 'to': v}
    return diff


def render_preset_manager(apply_callback, get_current_params_callback):
    """Render presets expander.

    Parameters
    ----------
    apply_callback: callable(name:str, values:dict) -> None
        Function invoked when user clicks Apply.
    get_current_params_callback: callable() -> dict
        Returns flattened current parameter map for diff display.
    """
    if 'saved_parameter_presets' not in st.session_state:
        return
    presets = st.session_state.get('saved_parameter_presets', {})
    if not presets:
        return
    with st.expander("💾 Preset Manager", expanded=False):
        st.caption("Apply / rename / delete saved parameter presets.")
        current = get_current_params_callback() or {}
        for name, meta in list(presets.items()):
            values = meta.get('values', {})
            diff_map = _compute_diff(values, current)
            cols = st.columns([3, 2, 2, 2, 1])
            with cols[0]:
                new_name = st.text_input(
                    f"Preset Name {name}", value=name, key=f"preset_rename_{name}",
                    label_visibility="collapsed"
                )
                if new_name != name and new_name.strip():
                    # rename
                    st.session_state.saved_parameter_presets[new_name] = presets.pop(name)
                    st.experimental_rerun()
            with cols[1]:
                st.write(meta.get('created_at', ''))
            with cols[2]:
                st.write(f"{len(values)} params")
            with cols[3]:
                if diff_map:
                    st.write(f"Δ {len(diff_map)}")
                else:
                    st.write("—")
            with cols[4]:
                if st.button("Apply", key=f"apply_preset_{name}"):
                    apply_callback(name, values)
                    st.success(f"Applied {name}")
            # second row for delete
            d_cols = st.columns([6,1])
            with d_cols[1]:
                if st.button("✖", key=f"delete_preset_{name}", help="Delete preset"):
                    presets.pop(name, None)
                    st.experimental_rerun()
        st.caption("Δ = differing parameters vs current working set")
