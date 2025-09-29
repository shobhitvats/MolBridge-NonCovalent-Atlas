"""Utility to load and inject CSS styles for Streamlit app."""
from pathlib import Path
import streamlit as st

STYLE_DIR = Path(__file__).resolve().parent.parent / "assets" / "styles"

CORE_ORDER = ["base.css", "custom_select.css"]
THEME_DARK = "theme_dark.css"
THEME_LIGHT = "theme_light.css"


def load_core_styles():
    for fname in CORE_ORDER:
        fpath = STYLE_DIR / fname
        if fpath.exists():
            st.markdown(f"<style>{fpath.read_text()}</style>", unsafe_allow_html=True)


def apply_theme(theme: str):
    """Apply theme-specific stylesheet. theme should be 'Dark' or 'Light'."""
    if theme == "Dark":
        f = STYLE_DIR / THEME_DARK
    else:
        f = STYLE_DIR / THEME_LIGHT
    if f.exists():
        st.markdown(f"<style>{f.read_text()}</style>", unsafe_allow_html=True)


def inject_all(theme: str):
    load_core_styles()
    apply_theme(theme)
