"""Centralized color palette & legend utilities.

Provides a colorblind-friendly categorical palette shared across
all interaction visualizations. Also renders a lightweight legend.

Future Enhancements:
 - Load palette overrides from user config or plugin manifests
 - Provide palette validation (contrast ratios)
 - Add toggle for alternative (dark/light) variants
"""
from __future__ import annotations

import streamlit as st

# Colorblind-friendly qualitative palette (derived from multiple accessible schemes)
# Keys are normalized (lowercase, no spaces) interaction identifiers
COLORBLIND_PALETTE = {
    'hydrogenbonds': '#1b9e77',
    'halogenbonds': '#d95f02',
    'pipi': '#7570b3',
    'ionicinteractions': '#e7298a',
    'hydrophobiccontacts': '#66a61e',
    'chpi': '#e6ab02',
    'chalcogenbonds': '#a6761d',
    'pnictogenbonds': '#666666',
    'tetrelbonds': '#1f78b4',
    'anionpi': '#b2df8a',
    'npistar': '#fb9a99',
    'londondispersion': '#fdbf6f',
    'cationpi': '#cab2d6',
    'saltbridges': '#ffff99',
    'sulfurpi': '#b15928',
    'metalcoordination': '#8dd3c7'
}

def normalize_key(name: str) -> str:
    return name.lower().replace(' ', '').replace('-', '')

def get_color(interaction_key: str) -> str:
    return COLORBLIND_PALETTE.get(normalize_key(interaction_key), '#999999')

def render_interaction_legend(selected=None, title: str = "Interaction Legend"):
    """Render a compact flexbox legend of interaction categories.

    Args:
        selected: Optional iterable of keys to emphasize (others dimmed)
        title: Legend title displayed above chips
    """
    selected = set(normalize_key(s) for s in (selected or []))
    st.markdown(f"#### {title}")
    st.markdown('<div style="display:flex;flex-wrap:wrap;gap:6px;">', unsafe_allow_html=True)
    for key, color in COLORBLIND_PALETTE.items():
        opacity = 1.0 if not selected or key in selected else 0.35
        pill = (
            f"<div style='background:{color};opacity:{opacity};padding:4px 9px;"\
            "border-radius:14px;font-size:0.65rem;font-weight:600;"\
            "font-family:system-ui,Segoe UI,Roboto,sans-serif;color:#000;"\
            "border:1px solid rgba(0,0,0,0.15);'>"\
            f"{key}</div>"
        )
        st.markdown(pill, unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)
