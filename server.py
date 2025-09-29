"""Protein Interaction Explorer - Main Streamlit Application
Production-ready server for analyzing noncovalent interactions in proteins.

Environment tweaks for resource-constrained dev systems:
- Force Streamlit file watcher to 'poll' (matches config.toml) to avoid inotify instance exhaustion.
"""

import sys
import os
os.environ.setdefault("STREAMLIT_SERVER_FILE_WATCHER_TYPE", "poll")
import streamlit as st
from pathlib import Path

# Add src to path for imports
sys.path.append(str(Path(__file__).parent / "src"))

from utils.config import AppConfig, load_config
try:
    from app.ui import style_loader  # new modular style loader
except Exception:
    style_loader = None  # fallback if refactor path not yet active
from utils.session import SessionManager
from utils.cache import CacheManager
from ui.main_interface import MainInterface

# Configure Streamlit page
st.set_page_config(
    page_title="MolBridge: NonCovalent Atlas",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://github.com/shobhitvats/Protein-Interaction-Analysis-Server',
        'Report a bug': 'https://github.com/shobhitvats/Protein-Interaction-Analysis-Server/issues',
    'About': """
    # MolBridge
    ### NonCovalent Atlas

    Charting the landscape of noncovalent interactions...    
    A comprehensive tool for analyzing noncovalent interactions in protein structures.
        
        **Features:**
        - Batch analysis of PDB structures
        - 11 types of noncovalent interactions
        - Interactive 3D visualization
        - Comprehensive reporting
        - REST API access
        
        **Version:** 1.0.0
        """
    }
)

def initialize_app():
    """Initialize application components."""
    if 'app_initialized' not in st.session_state:
        # Initialize configuration
        config = load_config()
        st.session_state.config = config
        
        # Initialize session manager
        st.session_state.session_manager = SessionManager()
        
        # Initialize cache manager
        st.session_state.cache_manager = CacheManager(config.cache_dir)
        
        # Initialize main interface
        st.session_state.main_interface = MainInterface(config)
        
        st.session_state.app_initialized = True

def main():
    """Main application entry point."""
    initialize_app()
    # Early style injection (core + theme) before heavy UI render if loader available
    if style_loader:
        # Default to dark theme until user picks; user theme radio later can re-apply
        style_loader.inject_all(theme=st.session_state.get('theme', 'Dark'))
    
    # Custom CSS with dropdown text visibility fixes (Dark theme defaults)
    st.markdown("""
    <style>
    .main > div { padding-top: 2rem; }
    .stAlert { margin-top: 1rem; }
    
    /* Labels: bold and white */
    .stSelectbox label,
    .stMultiSelect label {
        color: #F5F6FA !important;
        font-weight: 700 !important;
        font-size: 0.98rem !important;
        display: block !important;
    }

    /* Select container (BaseWeb) */
    .stSelectbox > div > div[data-baseweb="select"],
    .stMultiSelect > div > div[data-baseweb="select"] {
        background: rgba(36,37,38,0.85) !important;
        color: #F5F6FA !important;
        border: 1px solid rgba(0,180,216,0.35) !important;
        border-radius: 8px !important;
        min-height: 46px !important;
        display: flex !important;
        align-items: center !important;
        padding: 0 10px !important;
        position: relative !important; /* for absolute centering of value */
    }

    /* Make the inner value container fill available space and center its content */
    .stSelectbox > div > div[data-baseweb="select"] > div:first-child,
    .stMultiSelect > div > div[data-baseweb="select"] > div:first-child {
        flex: 1 1 auto !important;
        display: flex !important;
        align-items: center !important;
        min-height: 28px !important;
    }
    /* Keep the chevron compact */
    .stSelectbox > div > div[data-baseweb="select"] svg,
    .stMultiSelect > div > div[data-baseweb="select"] svg {
        margin-left: 6px !important;
        flex: 0 0 auto !important;
    }

    /* Text inside the select (selected value, input, spans) */
    .stSelectbox > div > div[data-baseweb="select"] > div,
    .stSelectbox > div > div[data-baseweb="select"] span,
    .stSelectbox > div > div[data-baseweb="select"] div,
    .stSelectbox [class*="singleValue"],
    .stSelectbox [class*="value"],
    .stSelectbox [class*="placeholder"],
    .stSelectbox input,
    .stMultiSelect > div > div[data-baseweb="select"] > div,
    .stMultiSelect > div > div[data-baseweb="select"] span,
    .stMultiSelect > div > div[data-baseweb="select"] div,
    .stMultiSelect input {
        color: #F5F6FA !important;
        -webkit-text-fill-color: #F5F6FA !important;
        background: transparent !important;
        font-weight: 600 !important;
        text-shadow: none !important;
        opacity: 1 !important;
    }

    /* Caret for text input */
    .stSelectbox input, .stMultiSelect input {
        caret-color: #F5F6FA !important;
    }

    /* Dropdown menu (portal) global fallback */
    ul[role="listbox"],
    .stSelectbox div[role="listbox"],
    .stMultiSelect div[role="listbox"] {
        background: rgba(36,37,38,0.95) !important;
        border: 1px solid rgba(0,180,216,0.35) !important;
        border-radius: 8px !important;
        box-shadow: 0 4px 18px rgba(0,0,0,0.35) !important;
    }
    ul[role="listbox"] li,
    div[role="option"],
    div[role="option"] * {
        color: #F5F6FA !important;
        -webkit-text-fill-color: #F5F6FA !important;
        background: rgba(36,37,38,0.85) !important;
    }
    div[role="option"]:hover, ul[role="listbox"] li:hover {
        background: rgba(0,180,216,0.3) !important;
        color: #FFFFFF !important;
    }

    /* Chevron icon color */
    .stSelectbox svg, .stMultiSelect svg {
        fill: #F5F6FA !important;
        color: #F5F6FA !important;
    }

    /* Focus ring */
    .stSelectbox > div > div[data-baseweb="select"]:focus-within,
    .stMultiSelect > div > div[data-baseweb="select"]:focus-within {
        border-color: #00B4D8 !important;
        box-shadow: 0 0 0 2px rgba(0,180,216,0.3) !important;
    }

    /* --- Additional hardening for cases where selected text appears invisible --- */
    /* Ensure the value container and its text are not scaled or faded by BaseWeb animations */
    .stSelectbox > div > div[data-baseweb="select"],
    .stMultiSelect > div > div[data-baseweb="select"] {
        /* Keep internal layout clipped so hidden sizers don't leak */
        overflow: hidden !important;
    }
    /* Hide BaseWeb measurement/sizer nodes which use a single '.' and aria-hidden="true" */
    .stSelectbox [aria-hidden="true"],
    .stMultiSelect [aria-hidden="true"] {
        display: none !important;
        visibility: hidden !important;
        width: 0 !important;
        height: 0 !important;
        overflow: hidden !important;
    }
    .stSelectbox [role="combobox"],
    .stSelectbox [role="combobox"] * {
        color: #F5F6FA !important;
        -webkit-text-fill-color: #F5F6FA !important;
        opacity: 1 !important;
        list-style: none !important;
    }
    .stSelectbox div[data-baseweb="select"] [class*="Value"],
    .stSelectbox div[data-baseweb="select"] [class*="SingleValue"],
    .stSelectbox div[data-baseweb="select"] [class*="Placeholder"],
    .stSelectbox div[data-baseweb="select"] [class*="value"] {
        color: #F5F6FA !important;
        -webkit-text-fill-color: #F5F6FA !important;
        opacity: 1 !important;
        /* allow proper centering transforms where needed */
        filter: none !important;
        mix-blend-mode: normal !important;
        /* position will be set more specifically below for value/placeholder */
        z-index: 2 !important;
        white-space: nowrap !important;
        text-overflow: ellipsis !important;
        overflow: hidden !important;
        font-size: 0.98rem !important;
        font-weight: 600 !important;
        visibility: visible !important;
    }

    /* Restore BaseWeb absolute centering for the visible value and the placeholder */
    .stSelectbox div[data-baseweb="select"] [class*="SingleValue"],
    .stSelectbox div[data-baseweb="select"] [class*="Placeholder"] {
        position: absolute !important;
        top: 52% !important;               /* push reference point slightly lower */
        transform: translateY(-52%) !important; /* net effect: slight upward shift */
        margin: 0 !important;
        left: 12px !important; /* typical padding; adjust if needed */
        right: 36px !important; /* leave room for chevron */
    }

    /* Eliminate bullets/dots from any list or pseudo-element markers inside select */
    .stSelectbox > div > div[data-baseweb="select"] *,
    .stMultiSelect > div > div[data-baseweb="select"] * {
        list-style: none !important;
    }
    .stSelectbox > div > div[data-baseweb="select"] *::marker,
    .stMultiSelect > div > div[data-baseweb="select"] *::marker {
        content: '' !important;
        display: none !important;
    }
    .stSelectbox > div > div[data-baseweb="select"] *::before,
    .stSelectbox > div > div[data-baseweb="select"] *::after,
    .stMultiSelect > div > div[data-baseweb="select"] *::before,
    .stMultiSelect > div > div[data-baseweb="select"] *::after {
        content: none !important;
    }

    /* Ensure the value container participates in layout properly and centers content */
    .stSelectbox > div > div[data-baseweb="select"] > div,
    .stMultiSelect > div > div[data-baseweb="select"] > div {
        min-height: 40px !important; /* match control height */
        line-height: 24px !important;
        display: flex !important;
        align-items: center !important;
        /* Slight top-light padding to push text visually up */
        padding: 4px 10px 6px 10px !important;
        box-sizing: border-box !important;
    }
    /* Explicitly center the combobox layer (where value/input lives) */
    .stSelectbox div[data-baseweb="select"] [role="combobox"],
    .stMultiSelect div[data-baseweb="select"] [role="combobox"] {
        display: flex !important;
        align-items: center !important;
        min-height: 28px !important;
        height: 100% !important;
        position: relative !important;
        top: 0 !important;
    }

    /* Based on observed DOM:
       select[data-baseweb="select"] > div(wrapper) > div:first-child(left content) > div:first-child(value)
       Center the left content area and the selected value element */
    .stSelectbox > div > div[data-baseweb="select"] > div:first-child { /* wrapper */
        display: flex !important;
        align-items: center !important;
        min-height: 46px !important;
    }
    .stSelectbox > div > div[data-baseweb="select"] > div:first-child > div:first-child { /* left content block */
        display: flex !important;
        align-items: center !important;
        min-height: 46px !important;
    }
    .stSelectbox > div > div[data-baseweb="select"] > div:first-child > div:first-child > div:first-child { /* selected value element */
        display: inline-flex !important;
        align-items: center !important;
        margin: 0 !important;
        padding: 0 !important;
        transform: none !important;
        line-height: 1.25 !important;
    }

    /* Robust targeting of the selected value element by attribute, covering both observed DOM variants */
    .stSelectbox div[data-baseweb="select"] > div > div[value],
    .stSelectbox div[data-baseweb="select"] > div > div:first-child > div[value] {
        display: inline-flex !important;
        align-items: center !important;
        max-width: calc(100% - 28px) !important;
        white-space: nowrap !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        margin: 0 !important;
        padding: 0 !important;
        line-height: 1.25 !important;
        transform: translateY(-2px) !important; /* slight upward nudge */
    }
    .stSelectbox div[data-baseweb="select"] [class*="SingleValue"],
    .stSelectbox div[data-baseweb="select"] [class*="Value"] {
        display: inline-flex !important;
        align-items: center !important;
        max-width: calc(100% - 28px) !important; /* leave room for chevron */
        line-height: 1.25 !important;
        margin: 0 !important;
        padding: 0 !important;
        align-self: center !important;
    }
    /* Remove manual upward nudge - rely on absolute centering above */
    .stSelectbox div[data-baseweb="select"] [class*="Value"] {
        transform: none !important;
    }
    .stSelectbox div[data-baseweb="select"] input,
    .stMultiSelect div[data-baseweb="select"] input {
        margin: 0 !important;
        padding: 0 !important;
        line-height: 1.25 !important;
        align-self: center !important;
        transform: none !important;
    }

    /* Keep input underneath the value layer to prevent overlay artifacts */
    .stSelectbox div[data-baseweb="select"] input {
        position: relative !important;
        z-index: 1 !important;
        background: transparent !important;
    }

    /* Avoid bullet-like artifacts from list styles */
    .stSelectbox div[data-baseweb="select"] [class*="Value"],
    .stSelectbox div[data-baseweb="select"] [class*="SingleValue"],
    .stSelectbox div[data-baseweb="select"] [class*="Placeholder"] {
        list-style: none !important;
        text-indent: 0 !important;
    }
    }
    </style>
    """, unsafe_allow_html=True)
    
    # App header
    # st.title("🧬 Protein Interaction Explorer")
    # st.markdown("*Comprehensive analysis of noncovalent interactions in protein structures*")
    
    # Render main interface
    st.session_state.main_interface.render()

if __name__ == "__main__":
    main()
