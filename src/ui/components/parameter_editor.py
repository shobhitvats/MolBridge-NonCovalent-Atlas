"""Parameter Editor Component

Encapsulates rendering of sliders for a single interaction's parameters in a
compact 2-column responsive layout with min/max badges and optional help
overlays. The main interface supplies metadata list (from registry).
"""
from __future__ import annotations
from typing import List, Dict, Any
import streamlit as st


def render_parameter_sliders(itype: str, params: List[Dict[str, Any]], config, preset_values: Dict[str, float]):
    changed_any = False
    if not params:
        st.info("No tunable parameters exposed.")
        return changed_any
    # responsive columns (2) per row
    for meta in params:
        # decide column container
        cols = st.columns(2)
        with cols[0]:
            label = meta['label']
            help_txt = meta.get('help')
        # single value
        if 'fields' in meta:
            f_min, f_max = meta['fields']
            cur_min = getattr(config.interactions, f_min, meta['min'])
            cur_max = getattr(config.interactions, f_max, meta['max'])
            if cur_min > cur_max:
                cur_min, cur_max = cur_max, cur_min
            new_min, new_max = st.slider(
                label,
                meta['min'], meta['max'], (float(cur_min), float(cur_max)), meta['step'],
                key=f"param_{itype}_{f_min}_{f_max}"
            )
            if (new_min, new_max) != (cur_min, cur_max):
                if hasattr(config.interactions, f_min):
                    setattr(config.interactions, f_min, float(new_min))
                if hasattr(config.interactions, f_max):
                    setattr(config.interactions, f_max, float(new_max))
                st.session_state['parameter_mismatch'] = True
            baseline_min = preset_values.get(f_min)
            baseline_max = preset_values.get(f_max)
            if (baseline_min is not None and abs(float(baseline_min) - float(new_min)) > 1e-9) or \
               (baseline_max is not None and abs(float(baseline_max) - float(new_max)) > 1e-9):
                changed_any = True
            if help_txt:
                with st.expander("ⓘ Details", expanded=False):
                    st.markdown(help_txt)
        else:
            field = meta['field']
            current_val = getattr(config.interactions, field, meta['min'])
            new_val = st.slider(
                meta['label'], meta['min'], meta['max'], float(current_val), meta['step'],
                key=f"param_{itype}_{field}"
            )
            if new_val != current_val:
                if hasattr(config.interactions, field):
                    setattr(config.interactions, field, float(new_val))
                st.session_state['parameter_mismatch'] = True
            baseline_val = preset_values.get(field)
            if baseline_val is not None and abs(float(baseline_val) - float(new_val)) > 1e-9:
                changed_any = True
            if meta.get('help'):
                with st.expander("ⓘ Details", expanded=False):
                    st.markdown(meta['help'])
    return changed_any
