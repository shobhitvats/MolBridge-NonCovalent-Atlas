"""
Protein Interaction Explorer - Main Streamlit Application
Production-ready server for analyzing noncovalent interactions in proteins.
"""

import streamlit as st
import sys
import os
from pathlib import Path

# Add src to path for imports
sys.path.append(str(Path(__file__).parent / "src"))

from utils.config import AppConfig, load_config
from utils.session import SessionManager
from utils.cache import CacheManager
from ui.main_interface import MainInterface

# Configure Streamlit page
st.set_page_config(
    page_title="Protein Interaction Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://github.com/shobhitvats/Protein-Interaction-Analysis-Server',
        'Report a bug': 'https://github.com/shobhitvats/Protein-Interaction-Analysis-Server/issues',
        'About': """
        # Protein Interaction Explorer
        
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
    
    # Custom CSS
    st.markdown("""
    <style>
    .main > div {
        padding-top: 2rem;
    }
    .stAlert {
        margin-top: 1rem;
    }
    .interaction-chip {
        display: inline-block;
        background-color: #f0f2f6;
        border-radius: 16px;
        padding: 4px 12px;
        margin: 2px;
        font-size: 0.8rem;
        border: 1px solid #d0d0d0;
    }
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
        border-left: 4px solid #1f77b4;
    }
    </style>
    """, unsafe_allow_html=True)
    
    # App header
    st.title("🧬 Protein Interaction Explorer")
    st.markdown("*Comprehensive analysis of noncovalent interactions in protein structures*")
    
    # Render main interface
    st.session_state.main_interface.render()

if __name__ == "__main__":
    main()
