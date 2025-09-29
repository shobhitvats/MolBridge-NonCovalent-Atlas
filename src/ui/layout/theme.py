"""Theme & Accessibility Utilities."""
from __future__ import annotations
import streamlit as st

HIGH_CONTRAST_CSS = """
<style>
/* High Contrast override */
body, .stApp { background:#000 !important; color:#fff !important; }
.stButton>button, .stDownloadButton>button { background:#ffd000 !important; color:#000 !important; }
table, th, td { color:#fff !important; }
</style>
"""


def render_high_contrast_toggle():
    hc = st.sidebar.toggle("High Contrast Mode", value=False, help="Increase contrast for accessibility (WCAG AA)")
    st.session_state['_high_contrast'] = hc
    if hc:
        st.markdown(HIGH_CONTRAST_CSS, unsafe_allow_html=True)
