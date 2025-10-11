"""
Main user interface for Protein Interaction Explorer.
Enhanced with high-performance processing capabilities.
"""

import streamlit as st
import pandas as pd
from typing import List, Dict, Any, Optional
import time
import uuid
import io
import zipfile
import numpy as np
from loguru import logger
import os
import json
from pathlib import Path
try:
    from rapidfuzz import process as fuzz_process  # type: ignore
except Exception:  # optional dependency safeguard
    fuzz_process = None

from utils.config import AppConfig, get_interaction_types, get_interaction_display_names
from utils.pdb_handler import PDBHandler
from visualization.structure_viewer import StructureViewer
from visualization.plots import InteractionPlots
from reporting.report_generator import ReportGenerator
from analysis.batch_processor import HighPerformanceBatchProcessor
from performance.parallel_processor import HighPerformanceProcessor
from analysis.extensions import (
    compute_residue_profiles,
    compute_interface_analysis,
    compute_outliers,
    compute_provenance,
    compute_motifs,
    compute_secondary_structure,
    compute_sasa_bsa,
    compute_geometry_quality,
    compute_disulfides,
    compute_pockets,
    compute_conservation,
    compute_pi_pi_refinement,
    compute_hbond_subtypes,
)
from ui.parameters.registry import REGISTRY as PARAMETER_REGISTRY, get_parameters_for
from ui.components.parameter_editor import render_parameter_sliders
from ui.components.preset_manager import render_preset_manager
from ui.layout.theme import render_high_contrast_toggle
from ui.components.color_legend import render_interaction_legend


class MainInterface:
    """Main Streamlit UI for MolBridge."""

    def __init__(self, config: AppConfig):
        self.config = config
        # Schema / format versions
        self.LAYOUT_SNAPSHOT_VERSION = 1

        # Core components
        self.structure_viewer = StructureViewer(config)
        self.plots = InteractionPlots(config)
        self.report_generator = ReportGenerator(config)
        self.batch_processor = HighPerformanceBatchProcessor(config=config)
        self.pdb_handler = PDBHandler(config)
        self.performance_processor = HighPerformanceProcessor(config=config)

        # Session-state defaults
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = {}
        if 'current_pdb' not in st.session_state:
            st.session_state.current_pdb = None
        if 'selected_interactions' not in st.session_state:
            st.session_state.selected_interactions = {itype: False for itype in get_interaction_types()}
        if 'individual_strength_filters' not in st.session_state:
            st.session_state.individual_strength_filters = {itype: 'all' for itype in get_interaction_types()}
        if 'performance_metrics' not in st.session_state:
            st.session_state.performance_metrics = []
        if 'processing_mode' not in st.session_state:
            st.session_state.processing_mode = 'high_performance'
        if 'session_data' not in st.session_state:
            st.session_state.session_data = {'bookmarks': []}

    # ---------------- Utility / Helper Methods -----------------
    def _environment_dependency_check(self):
        """Warn about missing optional dependencies (soft-fail)."""
        if st.session_state.get('_deps_checked'):
            return
        optional = ['pandas', 'seaborn', 'plotly']
        missing = []
        for mod in optional:
            try:
                __import__(mod)
            except Exception:
                missing.append(mod)
        if missing:
            st.sidebar.info(
                f"Optional packages not found: {', '.join(missing)}. Some tables/plots will fallback to plain text.",
                icon="ℹ️"
            )
        st.session_state._deps_checked = True

    def _compute_aggregated_performance(self):
        """Aggregate detector timing & counts across all analyzed structures."""
        if not st.session_state.analysis_results:
            return None
        totals = {}
        counts = {}
        structures = 0
        for pdb_id, res in st.session_state.analysis_results.items():
            summary = res.get('summary', {})
            tmap = summary.get('detector_timings', {})
            cmap = summary.get('detector_counts', {})
            if tmap:
                structures += 1
            for k, v in tmap.items():
                totals[k] = totals.get(k, 0.0) + v
            for k, v in cmap.items():
                counts[k] = counts.get(k, 0) + v
        if not totals:
            return None
        grand_total = sum(totals.values())
        rows = []
        for det, t in sorted(totals.items(), key=lambda x: x[1], reverse=True):
            rows.append({
                'Detector': det,
                'Total Time (s)': round(t, 4),
                'Avg per Structure (s)': round(t/structures, 4) if structures else 0,
                'Total Count': counts.get(det, 0),
                'Overall Rate (int/s)': round(counts.get(det,0)/t,2) if t > 0 and counts.get(det,0) > 0 else 0,
                'Time %': (t / grand_total * 100) if grand_total > 0 else 0
            })
        return {
            'rows': rows,
            'structures': structures,
            'grand_total': grand_total
        }

    def _update_individual_strength_filter(self, itype: str):
        """Callback to update the central strength filter dictionary from the widget's state."""
        widget_key = f"strength_{itype}"
        if widget_key in st.session_state:
            st.session_state.individual_strength_filters[itype] = st.session_state[widget_key]

    def _update_interaction_filter(self, itype: str, widget_key: str):
        """Callback to update selection of interaction types from widget state."""
        if widget_key in st.session_state:
            st.session_state.selected_interactions[itype] = st.session_state[widget_key]

    def _render_interaction_strength_filters(self):
        """Render the integrated interaction and strength filter UI in the sidebar."""
        st.subheader("2️⃣ Interactions & Strength")
        interaction_types = get_interaction_types()
        display_names = get_interaction_display_names()
        # Externalized parameter registry
        self._interaction_parameter_registry = PARAMETER_REGISTRY
        if 'selected_interactions' not in st.session_state:
            st.session_state.selected_interactions = {itype: False for itype in interaction_types}
        strength_options = {'strong': 'Strong','strong_moderate': 'S+M','all': 'All'}
        strength_keys = list(strength_options.keys())
        col1, col2 = st.columns(2)
        with col1:
            if st.button("✅ All", help="Select all interaction types", key="btn_all_interactions"):
                # Set both model state and widget state so checkboxes visually update
                for itype in interaction_types:
                    st.session_state.selected_interactions[itype] = True
                    wkey = f"sidebar_interaction_{itype}"
                    st.session_state[wkey] = True
                st.rerun()
        with col2:
            if st.button("🧹 None", help="Deselect all interaction types", key="btn_clear_interactions"):
                for itype in interaction_types:
                    st.session_state.selected_interactions[itype] = False
                    wkey = f"sidebar_interaction_{itype}"
                    st.session_state[wkey] = False
                st.rerun()
        # One-click revert all to current preset
        col_rev_all, col_save_layout = st.columns(2)
        with col_rev_all:
            if st.button("Revert All to Preset", key="revert_all_preset"):
                applied = getattr(self.config, 'applied_preset', 'literature_default')
                preset_vals = self.config.presets.get(applied, {}) if hasattr(self.config, 'presets') else {}
                for field, val in preset_vals.items():
                    if hasattr(self.config.interactions, field):
                        setattr(self.config.interactions, field, val)
                st.toast("All parameters reverted – re-run to apply", icon="↩️")
        with col_save_layout:
            if st.button("Save Layout", key="save_layout_snapshot"):
                import datetime, json, pathlib
                snapshot = {
                    'version': self.LAYOUT_SNAPSHOT_VERSION,
                    'selected_interactions': st.session_state.selected_interactions,
                    'individual_strength_filters': st.session_state.individual_strength_filters
                }
                pathlib.Path('cache').mkdir(exist_ok=True, parents=True)
                (pathlib.Path('cache')/f"layout_{datetime.datetime.now().strftime('%H%M%S')}.json").write_text(json.dumps(snapshot, indent=2))
                st.toast("Layout saved", icon="💾")
        for itype in interaction_types:
            c1, c2 = st.columns([3, 2])
            with c1:
                widget_key = f"sidebar_interaction_{itype}"
                # Parameter modified indicator (●)
                modified = False
                params = get_parameters_for(itype)
                applied = getattr(self.config, 'applied_preset', 'literature_default')
                preset_vals = self.config.presets.get(applied, {}) if hasattr(self.config, 'presets') else {}
                for m in params:
                    field = m.get('field')
                    if field and field in preset_vals and hasattr(self.config.interactions, field):
                        curv = getattr(self.config.interactions, field)
                        target = preset_vals[field]
                        if isinstance(curv, (int,float)) and isinstance(target,(int,float)):
                            if abs(float(curv)-float(target))>1e-9:
                                modified = True; break
                        elif curv != target:
                            modified = True; break
                label_txt = ("● " if modified else "") + display_names[itype]
                st.checkbox(label_txt,
                            value=st.session_state.selected_interactions.get(itype, False),
                            key=widget_key,
                            on_change=self._update_interaction_filter,
                            args=(itype, widget_key))
            if st.session_state.selected_interactions.get(itype, False):
                with c2:
                    default_value = st.session_state.individual_strength_filters.get(itype, 'all')
                    current_index = strength_keys.index(default_value)
                    st.radio(
                        f"Strength for {display_names[itype]}",
                        options=strength_keys,
                        format_func=lambda x: strength_options[x],
                        index=current_index,
                        key=f"strength_{itype}",
                        on_change=self._update_individual_strength_filter,
                        args=(itype,),
                        label_visibility="collapsed",
                    )
                # Collapsible parameter section for selected interaction
                params = get_parameters_for(itype)
                if params:
                    with st.expander(f"⚙️ Parameters – {display_names[itype]}", expanded=False):
                        applied_preset = getattr(self.config, 'applied_preset', 'literature_default')
                        preset_values = self.config.presets.get(applied_preset, {}) if hasattr(self.config, 'presets') else {}
                        changed_any = render_parameter_sliders(itype, params, self.config, preset_values)
                        cols_btn = st.columns(3)
                        with cols_btn[0]:
                            if st.button(f"Apply", key=f"apply_{itype}"):
                                st.session_state['parameter_mismatch'] = False
                                st.toast("Parameters applied – re-run to reflect in results", icon="✅")
                        with cols_btn[1]:
                            if st.button("Save Preset", key=f"save_preset_{itype}"):
                                import datetime
                                snapshot = {m['field']: getattr(self.config.interactions, m['field']) for m in params if 'field' in m and hasattr(self.config.interactions, m['field'])}
                                if 'saved_parameter_presets' not in st.session_state:
                                    st.session_state.saved_parameter_presets = {}
                                pname = f"{itype}_{datetime.datetime.now().strftime('%H%M%S')}"
                                st.session_state.saved_parameter_presets[pname] = {'values': snapshot, 'created_at': datetime.datetime.now().isoformat(timespec='seconds')}
                                st.toast(f"Saved preset {pname}")
                        with cols_btn[2]:
                            if st.button("Revert", key=f"revert_{itype}"):
                                for m in params:
                                    if 'field' in m and m['field'] in preset_values:
                                        setattr(self.config.interactions, m['field'], preset_values[m['field']])
                                st.info("Reverted – re-run required")
                        st.caption("🔧 Modified" if changed_any else "✅ Matches preset")

    def render(self):
        """Render the main interface with a modular, modern dashboard layout."""
        # One-time optional dependency check
        self._environment_dependency_check()
        # Always-on dark theme (light mode removed per user request)
        with st.sidebar:
            st.session_state["theme"] = "Dark"
            if 'use_classic_region_colors' not in st.session_state:
                st.session_state.use_classic_region_colors = True
            st.session_state.use_classic_region_colors = st.checkbox(
                "Classic Ramachandran Colors",
                value=st.session_state.use_classic_region_colors,
                help="Toggle between classic (green/yellow/red) and unified palette (if implemented later)"
            )
            st.markdown(
                """
                <style>
                body, .stApp { background-color: #18191A !important; color: #F5F6FA !important; }
                .stButton>button, .stDownloadButton>button { background-color: #00B4D8 !important; color: #18191A !important; border-radius: 8px; font-weight: 600; }
                .stTabs [data-baseweb=\"tab\"] { background: #242526 !important; color: #F5F6FA !important; }
                .modular-card { background: #22242a; border-radius: 18px; box-shadow: 0 4px 24px 0 #0003; padding: 2.5rem 2rem 2rem 2rem; margin-bottom: 2rem; }
                .hero-section { background: linear-gradient(90deg, #00B4D8 0%, #18191A 100%); border-radius: 18px; box-shadow: 0 4px 24px 0 #0003; padding: 2.5rem 2rem 2rem 2rem; margin-bottom: 2rem; color: #fff; }
                </style>
                """,
                unsafe_allow_html=True
            )

        # Sidebar for navigation and settings
        self._render_sidebar()

        # Hero CSS
        st.markdown(
            """<style>
            body, .stApp { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, 'Noto Sans', 'Liberation Sans', sans-serif !important; background: #18191A !important; }
            .stApp::before { content: ''; position: fixed; inset: 0; z-index: 0; background: radial-gradient(1200px 600px at 10% 10%, rgba(0,180,216,0.10), transparent 60%), radial-gradient(1000px 500px at 80% 20%, rgba(72,202,228,0.08), transparent 60%), linear-gradient(180deg, rgba(24,25,26,0.95) 0%, rgba(24,25,26,0.95) 100%); pointer-events: none; }
            .glass-hero { position: relative; background: linear-gradient(135deg, rgba(24,25,26,0.85) 0%, rgba(0,180,216,0.10) 100%); border-radius: 2.5rem; box-shadow: 0 4px 16px 0 #00B4D850; backdrop-filter: blur(12px); -webkit-backdrop-filter: blur(12px); padding: 3.5rem 2.5rem 2.5rem 2.5rem; margin-bottom: 2.5rem; overflow: hidden; animation: floatHero 3.5s cubic-bezier(.23,1.01,.32,1) 0.2s both; }
            .glass-hero .wave-bg { position: absolute; left: 0; right: 0; top: 0; height: 120px; z-index: 0; pointer-events: none; animation: waveMove 8s linear infinite alternate; }
            .glass-hero-content { position: relative; z-index: 1; display: flex; align-items: center; gap: 2.5rem; }
            .glass-title { font-size: clamp(2.2rem, 4.5vw, 3.2rem); font-weight: 800; line-height: 1.1; letter-spacing: 0.02em; color: #F5F6FA; text-shadow: 0 3px 14px rgba(0,180,216,0.45); }
            .glass-subtitle { margin-top: 0.35rem; font-size: clamp(1.0rem, 1.6vw, 1.25rem); font-weight: 500; color: #E2F3F9; opacity: 0.95; }
            .glass-subtitle .subtitle-prefix { font-weight: 800; color: #FFFFFF; display: inline; }
            .glass-subtitle .subtitle-main { display: inline; }
            .glass-subtitle .subtitle-line2 { display: block; margin-top: 0.18rem; opacity: 0.95; font-size: 1.01em; font-weight: 400; letter-spacing: 0.01em; margin-left: calc(1.1em * 2 + 13ch); text-indent: 0; }
            .glass-logo { width: 90px; height: 90px; border-radius: 50%; background: linear-gradient(135deg, #00B4D8 0%, #48CAE4 100%); box-shadow: 0 2px 8px #00B4D850; display: flex; align-items: center; justify-content: center; backdrop-filter: blur(4px); animation: pulseGlow 2.5s infinite alternate; }
            @keyframes floatHero { 0% { transform: translateY(40px) scale(0.98); opacity: 0; } 100% { transform: translateY(0) scale(1); opacity: 1; } }
            @keyframes waveMove { 0% { transform: translateY(0); } 100% { transform: translateY(10px) scaleX(1.03); } }
            @keyframes pulseGlow { 0% { box-shadow: 0 2px 8px #00B4D850; } 100% { box-shadow: 0 4px 16px #00B4D8aa; } }
            </style>""",
            unsafe_allow_html=True,
        )

        # Hero HTML
        st.markdown(
            """
            <div class=\"glass-hero\">
                <div class=\"wave-bg\">
                    <svg width=\"100%\" height=\"120\" viewBox=\"0 0 1440 120\" fill=\"none\" xmlns=\"http://www.w3.org/2000/svg\">
                        <path d=\"M0,80 C360,160 1080,0 1440,80 L1440,120 L0,120 Z\" fill=\"#00B4D8\" fill-opacity=\"0.18\"/>
                        <path d=\"M0,100 C400,0 1040,200 1440,100 L1440,120 L0,120 Z\" fill=\"#00B4D8\" fill-opacity=\"0.12\"/>
                    </svg>
                </div>
                <div class=\"glass-hero-content\">
                    <div class=\"glass-logo\">
                        <img src=\"https://img.icons8.com/color/96/000000/dna-helix.png\" width=\"60\" height=\"60\" alt=\"Protein Logo\"/>
                    </div>
                    <div>
                        <div class=\"glass-title\">MolBridge</div>
                        <div class=\"glass-subtitle\">
                            <span class=\"subtitle-prefix\">NonCovalent Atlas — </span><span class=\"subtitle-main\">Charting the landscape of noncovalent interactions...</span><br>
                            <span class=\"subtitle-line2\">&nbsp;&nbsp;&nbsp;Advanced protein structure & interaction analysis</span>
                        </div>
                    </div>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
        # Quick README / Docs links just below hero
        st.markdown(
            """
            <div style='margin-top:-1.2rem; margin-bottom:2.2rem; text-align:center; display:flex; gap:0.85rem; justify-content:center; flex-wrap:wrap;'>
                <a href=\"https://github.com/shobhitvats/Protein-Interaction-Analysis-Server?tab=readme-ov-file#molbridge-noncovalent-atlas\" target=\"_blank\" style=\"text-decoration:none;\" aria-label=\"Open README on GitHub\">
                    <span style=\"background:linear-gradient(90deg,#00B4D8,#48CAE4); color:#18191A; padding:0.55rem 1.15rem; border-radius:999px; font-weight:700; font-size:0.94rem; box-shadow:0 2px 10px #00B4D880; display:inline-block; white-space:nowrap;\">
                        📘 README
                    </span>
                </a>
                <a href=\"https://github.com/shobhitvats/Protein-Interaction-Analysis-Server/blob/main/docs/SCIENTIFIC_DOCUMENTATION.md#molbridge-scientific-documentation\" target=\"_blank\" style=\"text-decoration:none;\" aria-label=\"Open Scientific Documentation\">
                    <span style=\"background:linear-gradient(90deg,#48CAE4,#90E0EF); color:#18191A; padding:0.55rem 1.15rem; border-radius:999px; font-weight:700; font-size:0.94rem; box-shadow:0 2px 10px #48CAE480; display:inline-block; white-space:nowrap;\">
                        🧪 Scientific Criteria
                    </span>
                </a>
                <a href=\"https://github.com/shobhitvats/Protein-Interaction-Analysis-Server/blob/main/docs/TECHNICAL_DOCUMENTATION.md#molbridge-technical-documentation\" target=\"_blank\" style=\"text-decoration:none;\" aria-label=\"Open Technical Documentation\">
                    <span style=\"background:linear-gradient(90deg,#0096C7,#00B4D8); color:#18191A; padding:0.55rem 1.15rem; border-radius:999px; font-weight:700; font-size:0.94rem; box-shadow:0 2px 10px #0096C770; display:inline-block; white-space:nowrap;\">
                        🛠️ Technical Docs
                    </span>
                </a>
            </div>
            """,
            unsafe_allow_html=True,
        )

        # Cards and controls CSS
        st.markdown("""
        <style>
        .glass-card { background: rgba(36,37,38,0.82); border-radius: 1.5rem; box-shadow: 0 8px 32px 0 #00B4D880; backdrop-filter: blur(14px); -webkit-backdrop-filter: blur(14px); padding: 2.5rem 2rem 2rem 2rem; margin-bottom: 2.5rem; border: 1.5px solid rgba(0,180,216,0.18); }
        .gradient-btn button, .stButton>button, .stDownloadButton>button, .stFormSubmitButton>button { background: linear-gradient(90deg, #00B4D8 0%, #48CAE4 100%) !important; color: #18191A !important; border-radius: 12px !important; font-weight: 800 !important; font-size: 1.13rem !important; box-shadow: 0 2px 12px #00B4D880; border: none !important; padding: 0.7rem 1.7rem !important; margin: 0.2rem; letter-spacing: 0.02em; transition: background 0.3s, box-shadow 0.3s, transform 0.15s; }
        .gradient-btn button:hover, .stButton>button:hover, .stDownloadButton>button:hover, .stFormSubmitButton>button:hover { background: linear-gradient(90deg, #48CAE4 0%, #00B4D8 100%) !important; box-shadow: 0 6px 24px #00B4D8cc; transform: translateY(-2px) scale(1.04); }
        .stSelectbox>div>div, .stMultiSelect>div>div, .stTextInput>div>input, .stTextArea>div>textarea { background: rgba(36,37,38,0.55) !important; border-radius: 10px !important; color: #F5F6FA !important; border: 1.5px solid rgba(0,180,216,0.18) !important; font-size: 0.96rem !important; box-shadow: 0 1px 6px #00B4D830; padding: 0.35rem 0.75rem !important; min-height: 40px !important; line-height: 24px !important; }
        .stSelectbox>div>div[data-baseweb='select'] *::marker { content: none !important; }
        .stSelectbox>div>div[data-baseweb='select'] *::before, .stSelectbox>div>div[data-baseweb='select'] *::after { content: none !important; }
        .stSelectbox div[data-baseweb='select'] [class*='SingleValue'], .stSelectbox div[data-baseweb='select'] [class*='Value'] { max-width: calc(100% - 28px) !important; white-space: nowrap !important; overflow: hidden !important; text-overflow: ellipsis !important; }
        .stSelectbox>div>div:focus, .stMultiSelect>div>div:focus, .stTextInput>div>input:focus, .stTextArea>div>textarea:focus { border: 1.5px solid #00B4D8 !important; box-shadow: 0 0 0 2px #00B4D880; }
        .stSelectbox div[data-baseweb='select'] > div, .stSelectbox div[data-baseweb='select'] span, .stSelectbox div[data-baseweb='select'] div > div, .stSelectbox [role='combobox'], .stSelectbox [role='option'] { color: #F5F6FA !important; background: rgba(36,37,38,0.85) !important; }
        .stSelectbox div[data-baseweb='popover'] { background: rgba(36,37,38,0.95) !important; border: 1px solid #00B4D8 !important; border-radius: 8px !important; }
        .stSelectbox ul[role='listbox'] li { color: #F5F6FA !important; background: rgba(36,37,38,0.85) !important; }
        .stSelectbox ul[role='listbox'] li:hover { background: rgba(0,180,216,0.3) !important; color: #FFFFFF !important; }
        .stCheckbox>label>div:first-child { border-radius: 6px !important; border: 1.5px solid #00B4D8 !important; background: rgba(0,180,216,0.12) !important; }
        .stRadio>div>label { border-radius: 8px !important; background: rgba(0,180,216,0.10) !important; color: #F5F6FA !important; padding: 0.3rem 0.8rem !important; margin: 0.1rem 0.2rem !important; font-weight: 600 !important; }
        .stRadio>div>label[data-selected='true'] { background: linear-gradient(90deg, #00B4D8 0%, #48CAE4 100%) !important; color: #18191A !important; }
        .section-divider { height: 2px; background: linear-gradient(90deg, #00B4D8 0%, #18191A 100%); border: none; margin: 2.5rem 0 2rem 0; border-radius: 2px; }
        .section-header { font-size: 1.5rem; font-weight: 700; color: #00B4D8; margin-bottom: 1.2rem; display: flex; align-items: center; gap: 0.7rem; }
        </style>
        """, unsafe_allow_html=True)

        # Tabs
        tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs([
            "🔍 Analysis","📊 Visualization","📋 Results","📄 Reports",
            "🧪 Structure Quality","🔄 Batch Comparison","⚙️ Settings","ℹ️ Info"
        ]) 
        with tab1:
            self._render_analysis_tab()
        with tab2:
            self._render_visualization_tab()
        with tab3:
            self._render_results_tab()
        with tab4:
            self._render_reports_tab()
        with tab5:
            self._render_structure_quality_tab()
        with tab6:
            self._render_batch_comparison_tab()
        with tab7:
            self._render_settings_tab()
        with tab8:
            self._render_info_tab()

    def _render_structure_quality_tab(self):
        """Render structural quality & annotation extensions in organized expanders."""
        st.header("🧪 Structural Quality & Annotations")
        if not st.session_state.analysis_results:
            st.info("Run an analysis first to see structural quality metrics.")
            return
        if not st.session_state.current_pdb:
            st.info("Select or load a structure to view its structural annotations.")
            return
        pdb_id = st.session_state.current_pdb
        result = st.session_state.analysis_results.get(pdb_id)
        if not result:
            st.warning("No result for current structure.")
            return
        # Ensure extensions (lazy compute)
        self._ensure_extensions(pdb_id, result)
        exts = result.get('extensions', {})
        if not exts:
            st.info("No structural extensions computed (check settings toggles).")
            return
        # Secondary Structure
        if 'secondary_structure' in exts:
            with st.expander("Secondary Structure & Torsions", expanded=True):
                ss = exts['secondary_structure']
                counts = ss.get('counts', {})
                fractions = ss.get('fractions', {})
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Helix Residues", counts.get('H',0))
                    st.metric("Sheet Residues", counts.get('E',0))
                    st.metric("Coil Residues", counts.get('C',0))
                with col2:
                    st.metric("Helix %", f"{fractions.get('H_frac',0)*100:.1f}%")
                    st.metric("Sheet %", f"{fractions.get('E_frac',0)*100:.1f}%")
                    st.metric("Coil %", f"{fractions.get('C_frac',0)*100:.1f}%")
                torsions = ss.get('torsions', [])
                if torsions:
                    import pandas as pd
                    df_t = pd.DataFrame(torsions)[:500]
                    st.caption("First 500 torsion rows (phi/psi/omega, assigned SS)")
                    st.dataframe(df_t, width='stretch')
        # Geometry Quality
        if 'geometry_quality' in exts:
            with st.expander("Geometry Quality (Ramachandran, Clashes)", expanded=True):
                gq = exts['geometry_quality']
                rama = gq.get('ramachandran', {})
                clashes = gq.get('clashes', [])
                col1, col2, col3 = st.columns(3)
                col1.metric("Favored", rama.get('favored',0))
                col2.metric("Allowed", rama.get('allowed',0))
                col3.metric("Outliers", rama.get('outliers',0))
                torsion_outliers = gq.get('torsion_outliers', [])
                if torsion_outliers:
                    st.subheader("Torsion Outliers")
                    import pandas as pd
                    st.dataframe(pd.DataFrame(torsion_outliers), width='stretch')
                if clashes:
                    st.subheader(f"Steric Clashes ({len(clashes)})")
                    import pandas as pd
                    st.dataframe(pd.DataFrame(clashes)[:300], width='stretch')
        # Surface & Interface
        if 'sasa_bsa' in exts or 'disulfides' in exts:
            with st.expander("Surface / Interface (SASA, BSA, Disulfides)", expanded=False):
                if 'sasa_bsa' in exts:
                    sas = exts['sasa_bsa']
                    chain_sasa = sas.get('chain_sasa', {})
                    if chain_sasa:
                        import pandas as pd
                        st.subheader("Chain SASA")
                        st.dataframe(pd.DataFrame([
                            {'Chain': c, 'SASA': v.get('sasa',0.0)} for c,v in chain_sasa.items()
                        ]), width='stretch')
                    bsa = sas.get('buried_surface', {})
                    if bsa:
                        st.subheader("Approx. Buried Surface (Pairwise)")
                        import pandas as pd
                        rows = []
                        for pair, val in bsa.items():
                            rows.append({'Pair': pair, 'Buried_SASA': val})
                        st.dataframe(pd.DataFrame(rows), width='stretch')
                if 'disulfides' in exts:
                    ds = exts['disulfides'].get('bonds', [])
                    st.subheader(f"Disulfide Bonds ({len(ds)})")
                    if ds:
                        import pandas as pd
                        st.dataframe(pd.DataFrame(ds), width='stretch')
        # Pockets
        if 'pockets' in exts:
            with st.expander("Predicted Pockets", expanded=False):
                pockets = exts['pockets'].get('pockets', [])
                st.caption("Heuristic grid clustering; interpret volumes comparatively")
                if pockets:
                    import pandas as pd
                    dfp = pd.DataFrame(pockets)
                    # Harmonize volume column naming (extension uses 'approx_volume')
                    if 'volume' not in dfp.columns and 'approx_volume' in dfp.columns:
                        dfp = dfp.rename(columns={'approx_volume': 'volume'})
                    # Only sort if column present
                    if 'volume' in dfp.columns:
                        dfp = dfp.sort_values('volume', ascending=False)
                    st.dataframe(dfp[:300], width='stretch')
        # Conservation (stub)
        if 'conservation' in exts:
            with st.expander("Residue Conservation (Approximate)", expanded=False):
                cons = exts['conservation'].get('residues', {})
                if cons:
                    import pandas as pd
                    rows = []
                    for res, rec in list(cons.items())[:500]:
                        rows.append({'Residue': res, 'Score': rec.get('score'), 'Rank': rec.get('rank')})
                    st.dataframe(pd.DataFrame(rows), width='stretch')
        # Pi-Pi refinement
        if 'pi_pi_refinement' in exts:
            with st.expander("π-π Stacking Refinement", expanded=False):
                ref = exts['pi_pi_refinement']
                subtypes = ref.get('subtype_counts') or ref.get('subtypes', {})
                total_pi = ref.get('count', sum(subtypes.values()) if isinstance(subtypes, dict) else 0)
                st.caption(f"Total π-π interactions: {total_pi}")
                col1, col2, col3 = st.columns(3)
                col1.metric("Parallel", subtypes.get('parallel',0))
                col2.metric("T-shaped", subtypes.get('t_shaped',0))
                col3.metric("Offset Parallel", subtypes.get('offset_parallel', subtypes.get('intermediate',0)))
                details = ref.get('details', [])
                if details:
                    import pandas as pd
                    st.dataframe(pd.DataFrame(details)[:300], width='stretch')

        # A11y enhancement script
        st.markdown("""
        <script>(function(){function s(t){return(t||'').toLowerCase().replace(/[^a-z0-9]+/g,'_').replace(/^_|_$/g,'')}function f(e){var l=e.closest('[data-testid]')?e.closest('[data-testid]').querySelector('label'):null;if(!l){var p=e.parentElement;while(p&&!l){if(p.querySelector)l=p.querySelector('label');p=p.parentElement}}return l}function patch(e){var label=e.getAttribute('aria-label')||'';if(!label){var L=f(e);if(L){label=(L.textContent||'').trim()}}if(label&&!e.getAttribute('aria-label'))e.setAttribute('aria-label',label);if(!e.id){e.id='fld_'+s(label||e.placeholder||'field')+'_'+Math.random().toString(36).slice(2,7)}if(!e.name){e.name=s(label||e.placeholder||e.id)}if(e.autocomplete==='')e.setAttribute('autocomplete','off');if(e.type==='password')e.setAttribute('autocomplete','new-password');if(/pdb id/i.test(label))e.setAttribute('autocomplete','off');var ln=f(e);if(ln&&!ln.getAttribute('for'))ln.setAttribute('for',e.id)}function run(r){r.querySelectorAll('input, select, textarea').forEach(patch)}function init(){run(document);var mo=new MutationObserver(function(m){m.forEach(function(rec){rec.addedNodes&&rec.addedNodes.forEach(function(n){if(n.nodeType===1)run(n)})})});mo.observe(document.body,{subtree:true,childList:true})}if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',init);else init();})();</script>
        """, unsafe_allow_html=True)

        # Final select centering CSS
        st.markdown("""
        <style>
        .stSelectbox > div > div[data-baseweb='select'], .stMultiSelect > div > div[data-baseweb='select'] { position: relative !important; min-height: 46px !important; display: flex !important; align-items: center !important; padding: 0 10px !important; }
        .stSelectbox > div > div[data-baseweb='select'] > div:first-child, .stSelectbox > div > div[data-baseweb='select'] > div:first-child > div:first-child { display: flex !important; align-items: center !important; min-height: 46px !important; }
        .stSelectbox div[data-baseweb='select'] > div > div[value], .stSelectbox div[data-baseweb='select'] > div:first-child > div[value], .stSelectbox div[data-baseweb='select'] > div:first-child > div:first-child > div[value] { display: inline-flex !important; align-items: center !important; max-width: calc(100% - 28px) !important; white-space: nowrap !important; overflow: hidden !important; text-overflow: ellipsis !important; margin: 0 !important; padding: 0 !important; line-height: 1.25 !important; transform: none !important; }
        </style>
        """, unsafe_allow_html=True)
        
        

    def _render_batch_comparison_tab(self):
        """Render batch comparison and advanced visualizations."""
        st.header("🔄 Batch Comparative Analysis")
        if len(st.session_state.analysis_results) < 2:
            st.info("Run a batch analysis to compare multiple structures.")
            return
        # Use the comparative analysis from plots
        self.plots.render_comparative_analysis(st.session_state.analysis_results)
        # Advanced: show interaction network for each structure
        st.write("---")
        st.subheader("🕸️ Interaction Networks (per structure)")
        for pdb_id, result in st.session_state.analysis_results.items():
            with st.expander(f"Network: {pdb_id}", expanded=False):
                self.plots.render_interaction_network(result)


    def _render_info_tab(self):
        """Render the info tab with interaction criteria via modular panel."""
        try:
            from app.ui.panels import info_panel
        except Exception:
            # Fallback if path differs or during transitional refactor
            try:
                from ui.panels import info_panel  # type: ignore
            except Exception:
                st.error("Info panel module not found.")
                return
        info_panel.render()
    
    def _render_sidebar(self):
        """Render the sidebar with input options and quick actions."""
        with st.sidebar:
            # Small screen / mobile adaptation toggle (heuristic)
            if 'is_small_screen' not in st.session_state:
                st.session_state.is_small_screen = False
            st.checkbox("Small Screen Mode", value=st.session_state.is_small_screen, key='is_small_screen', help="Condense UI for narrow viewports (manual toggle due to Streamlit sandbox).")
            if st.session_state.is_small_screen:
                st.caption("Compact mode active: heavy plotting expanders start collapsed.")
            render_high_contrast_toggle()
            st.header("1️⃣ Input")
            
            # Quick PDB Entry Section
            st.subheader("PDB Entry")
            with st.form("sidebar_pdb_form"):
                pdb_id = st.text_input(
                    "PDB ID:",
                    placeholder="e.g., 6SPO",
                    help="Enter a 4-character PDB identifier"
                )
                submitted = st.form_submit_button("🚀 Analyze", type="primary")
            
            if submitted and pdb_id:
                # Check if any interactions are selected
                selected_interactions = [k for k, v in st.session_state.get('selected_interactions', {}).items() if v]
                if selected_interactions:
                    self._run_single_analysis(pdb_id.upper().strip())
                else:
                    st.warning("⚠️ Please select at least one interaction type below.")
            
            st.divider()
            self._render_interaction_strength_filters()
            # Preset manager panel
            render_preset_manager(
                apply_callback=lambda name, vals: [setattr(self.config.interactions, k, v) for k, v in vals.items()],
                get_current_params_callback=lambda: {k: getattr(self.config.interactions, k) for k in dir(self.config.interactions) if not k.startswith('_')}
            )
            if st.session_state.get('parameter_mismatch'):
                st.markdown('<div style="background:#ffd000;color:#000;padding:4px 8px;border-radius:6px;font-size:0.75rem;font-weight:600;display:inline-block;">Modified params – re-run needed</div>', unsafe_allow_html=True)
            st.divider()
            
            st.header("3️⃣ Advanced Input")
            
            # Input mode selection
            input_mode = st.radio(
                "Choose input mode:",
                ["Single PDB ID", "Batch Upload", "Load Session"],
                help="Select how you want to provide protein structures for analysis"
            )
            
            if input_mode == "Single PDB ID":
                st.info("💡 Use quick entry above or Analysis tab for full options")
            elif input_mode == "Batch Upload":
                self._render_batch_input()
            else:
                self._render_session_loader()
            
            st.divider()
            st.subheader("4️⃣ Run & Settings")
            
            # Structure display options
            show_assembly = st.selectbox(
                "Structure Display:",
                ["Biological Assembly", "Asymmetric Unit"],
                index=0 if self.config.default_assembly == "biological" else 1
            )
            
            include_ligands = st.checkbox("Include Ligands", value=True)
            exclude_waters = st.checkbox("Exclude Waters", value=True)
            
            # Interaction preset
            preset = st.selectbox(
                "Analysis Preset:",
                ["Conservative", "Literature Default", "Exploratory"],
                index=1,
                help="Predefined parameter sets for different analysis stringency"
            )
            
            # Store settings in session state
            st.session_state.structure_settings = {
                "assembly": "biological" if "Biological" in show_assembly else "asymmetric",
                "include_ligands": include_ligands,
                "exclude_waters": exclude_waters,
                "preset": preset.lower().replace(" ", "_")
            }
            
            st.divider()
            st.subheader("5️⃣ Session")
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Save Session", width="stretch"):
                    session_id = st.session_state.session_manager.save_session()
                    if session_id:
                        st.success("Session saved!")
            
            with col2:
                if st.button("New Session", width="stretch"):
                    st.session_state.session_data = st.session_state.session_manager._create_empty_session()
                    st.session_state.analysis_results = {}
                    st.rerun()
            
            st.divider()
            st.subheader("6️⃣ Performance")
            if st.button("Dry-Run Preview", help="Geometry-only preview using dry_run; estimates candidate funnel without full materialization"):
                if not st.session_state.get('current_pdb'):
                    st.info("Run an analysis first to load structure context.")
                else:
                    structure = st.session_state.get('rama_source_structure')
                    if structure is None:
                        st.warning("Structure object not cached; re-run analysis once first.")
                    else:
                        from analysis.registry import DETECTOR_REGISTRY
                        det_classes = [cls for cls,_ in DETECTOR_REGISTRY.values()]
                        preview_proc = HighPerformanceProcessor(config=self.config)
                        progress_events = []
                        def _cb(ev):
                            progress_events.append(ev)
                        preview_proc.process_interactions_parallel(structure, det_classes, config=self.config, dry_run=True, progress_callback=_cb)
                        import pandas as pd
                        rows = []
                        for ev in progress_events:
                            f = ev.get('funnel') or {}
                            rows.append({
                                'detector': ev['detector'],
                                'raw': f.get('raw_pairs'),
                                'candidates': f.get('candidate_pairs'),
                                'accepted': f.get('accepted_pairs'),
                                'accept_ratio': f.get('acceptance_ratio'),
                                'duration': ev.get('duration')
                            })
                        if rows:
                            st.dataframe(pd.DataFrame(rows), use_container_width=True)
                        else:
                            st.info("No preview instrumentation available.")
            
            # Processing mode selector
            processing_mode = st.selectbox(
                "Processing Mode:",
                ["high_performance", "standard", "conservative"],
                index=0,
                format_func=lambda x: {
                    "high_performance": "🚀 High Performance (Max Speed)",
                    "standard": "⚖️ Standard (Balanced)",
                    "conservative": "🛡️ Conservative (Safe)"
                }[x],
                help="Select processing mode for analysis speed vs. safety"
            )
            st.session_state.processing_mode = processing_mode
            
            # Show current performance metrics
            if st.session_state.performance_metrics:
                recent_metrics = st.session_state.performance_metrics[-5:]  # Last 5 operations
                
                # Calculate average performance
                if recent_metrics:
                    avg_time = sum(m.get('analysis_time', m.get('batch_time', 0)) for m in recent_metrics) / len(recent_metrics)
                    total_proteins = sum(m.get('total_proteins', 1) for m in recent_metrics)
                    
                    with st.container():
                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("Avg Time", f"{avg_time:.1f}s")
                        with col2:
                            st.metric("Total Processed", total_proteins)
                
                # Performance trend indicator
                if len(recent_metrics) >= 2:
                    current_perf = recent_metrics[-1].get('analysis_time', recent_metrics[-1].get('batch_time', 0))
                    previous_perf = recent_metrics[-2].get('analysis_time', recent_metrics[-2].get('batch_time', 0))
                    
                    if current_perf < previous_perf:
                        st.success("📈 Performance improving!")
                    elif current_perf > previous_perf * 1.2:
                        st.warning("📉 Performance declining")
                
                # Global processor stats
                global_perf = self.performance_processor.get_performance_summary()
                if 'cache_size' in global_perf:
                    st.caption(f"🗄️ Cache: {global_perf['cache_size']} entries")
                    
                    if st.button("🧹 Clear Cache", help="Clear performance cache to free memory"):
                        self.performance_processor.clear_cache()
                        st.success("Cache cleared!")
                        time.sleep(1)
                        st.rerun()

            # Layout snapshot restore UI (previously only save existed in interaction panel)
            st.divider()
            st.subheader("🗂️ Layout Snapshots")
            try:
                layout_files = sorted(Path('cache').glob('layout_*.json'))
            except Exception:
                layout_files = []
            if layout_files:
                names = [f.name for f in layout_files]
                chosen = st.selectbox("Choose snapshot", names, key='layout_restore_select')
                col_ls1, col_ls2 = st.columns(2)
                with col_ls1:
                    if st.button("Restore Layout", key='restore_layout_btn'):
                        try:
                            payload = json.loads((Path('cache')/chosen).read_text())
                            sel = payload.get('selected_interactions'); strengths = payload.get('individual_strength_filters')
                            if isinstance(sel, dict):
                                for k,v in sel.items():
                                    if k in st.session_state.selected_interactions:
                                        st.session_state.selected_interactions[k] = v
                            if isinstance(strengths, dict):
                                for k,v in strengths.items():
                                    st.session_state.individual_strength_filters[k] = v
                            st.success("Layout restored – rerun to apply if needed.")
                        except Exception as e:
                            st.error(f"Restore failed: {e}")
                with col_ls2:
                    if st.button("Delete Snapshot", key='delete_layout_btn'):
                        try:
                            (Path('cache')/chosen).unlink()
                            st.experimental_rerun()
                        except Exception as e:
                            st.error(f"Delete failed: {e}")
            else:
                st.caption("No saved layouts yet (use 'Save Layout' in interactions panel).")

            # Scenario Profiles
            st.divider()
            st.subheader("🧪 Scenario Profiles")
            profiles = getattr(self.config, 'scenario_profiles', {}) or {}
            if profiles:
                prof_keys = list(profiles.keys())
                choice = st.selectbox("Profile", ["(none)"] + prof_keys, key='scenario_profile_select', format_func=lambda k: k.replace('_',' ').title())
                if choice != "(none)":
                    prof = profiles.get(choice, {})
                    csp1, csp2 = st.columns(2)
                    with csp1:
                        if st.button("Apply Profile", key='apply_profile_btn'):
                            ints = prof.get('interactions') or []
                            for it in st.session_state.selected_interactions.keys():
                                st.session_state.selected_interactions[it] = it in ints
                            preset = prof.get('preset')
                            if preset and preset in self.config.presets:
                                for f,v in self.config.presets[preset].items():
                                    if hasattr(self.config.interactions, f):
                                        setattr(self.config.interactions, f, v)
                                self.config.applied_preset = preset
                            st.success("Scenario profile applied – run analysis.")
                    with csp2:
                        if st.button("Explain", key='explain_profile_btn'):
                            st.info(prof.get('description','No description provided'))
                    if prof.get('notes'):
                        with st.expander("Profile Notes"):
                            st.markdown(prof['notes'])
            else:
                st.caption("No scenario profiles loaded (add templates/scenario_profiles.yaml)")

            st.caption("Press Ctrl/Cmd+K or use button below for Command Palette")

        # Outside sidebar: command palette modal / container
        self._render_command_palette()

    # ---------------- Command Palette & Actions -----------------
    def _command_items(self) -> List[Dict[str,str]]:
        items: List[Dict[str,str]] = []
        # Tabs navigation (labels only; actual navigation limited by Streamlit tab API)
        tab_items = [
            ("Go: Analysis", "nav:tab:analysis"),
            ("Go: Visualization", "nav:tab:visualization"),
            ("Go: Results", "nav:tab:results"),
            ("Go: Reports", "nav:tab:reports"),
            ("Go: Structure Quality", "nav:tab:quality"),
            ("Go: Batch Comparison", "nav:tab:batch"),
            ("Go: Settings", "nav:tab:settings"),
            ("Go: Info", "nav:tab:info")
        ]
        for label, cmd in tab_items:
            items.append({'label': label, 'command': cmd, 'group': 'Navigation'})
        # Presets
        for pname in self.config.presets.keys():
            items.append({'label': f"Apply Preset: {pname}", 'command': f"preset:{pname}", 'group': 'Presets'})
        # Scenario profiles
        for sname in getattr(self.config, 'scenario_profiles', {}).keys():
            items.append({'label': f"Apply Profile: {sname}", 'command': f"profile:{sname}", 'group': 'Profiles'})
        # Utility
        items.extend([
            {'label': 'Toggle Small Screen Mode', 'command': 'ui:toggle_small', 'group': 'UI'},
            {'label': 'Save Layout Snapshot', 'command': 'layout:save', 'group': 'Layout'},
            {'label': 'Restore Last Layout', 'command': 'layout:restore_last', 'group': 'Layout'},
        ])
        return items

    def _render_command_palette(self):
        if 'show_command_palette' not in st.session_state:
            st.session_state.show_command_palette = False
        # Inject JS to capture Ctrl/Cmd+K and click hidden button
        st.markdown(
            """
            <script>
            document.addEventListener('keydown', function(e){
              if((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k'){
                e.preventDefault();
                const streamlitDoc = window.parent.document;
                const btn = streamlitDoc.querySelector('button[data-command-palette-open="1"]');
                if(btn){ btn.click(); }
              }
            });
            </script>
            """,
            unsafe_allow_html=True
        )
        # Hidden open button (styled normally) with attribute for JS
        open_btn = st.button("⌨️ Command Palette", key='open_command_palette_btn', help="Ctrl/Cmd+K", kwargs={'data-command-palette-open': '1'})
        if open_btn:
            st.session_state.show_command_palette = True
        if not st.session_state.show_command_palette:
            return
        items = self._command_items()
        st.markdown("### 🔎 Command Palette")
        query = st.text_input("Type a command or search…", key='cmd_palette_query')
        matches = items
        if query:
            if fuzz_process:
                labels = [it['label'] for it in items]
                try:
                    scored = fuzz_process.extract(query, labels, limit=15)
                    lbls = {s[0] for s in scored}
                    matches = [it for it in items if it['label'] in lbls]
                except Exception:
                    ql = query.lower(); matches = [it for it in items if ql in it['label'].lower()]
            else:
                ql = query.lower(); matches = [it for it in items if ql in it['label'].lower()]
        if matches:
            for it in matches[:15]:
                if st.button(it['label'], key=f"cmd_{it['command']}"):
                    self._execute_command_palette_action(it['command'])
                    st.session_state.show_command_palette = False
                    st.experimental_rerun()
        else:
            st.caption("No matches")
        if st.button("Close", key='close_cmd_palette'):
            st.session_state.show_command_palette = False

    def _execute_command_palette_action(self, command: str):
        # Preset
        if command.startswith('preset:'):
            pname = command.split(':',1)[1]
            preset = self.config.presets.get(pname)
            if preset:
                for f,v in preset.items():
                    if hasattr(self.config.interactions, f):
                        setattr(self.config.interactions, f, v)
                self.config.applied_preset = pname
            return
        # Scenario profile
        if command.startswith('profile:'):
            prof_name = command.split(':',1)[1]
            prof = self.config.get_scenario_profile(prof_name)
            if prof:
                ints = prof.get('interactions') or []
                for it in st.session_state.selected_interactions.keys():
                    st.session_state.selected_interactions[it] = it in ints
                preset = prof.get('preset')
                if preset and preset in self.config.presets:
                    for f,v in self.config.presets[preset].items():
                        if hasattr(self.config.interactions, f):
                            setattr(self.config.interactions, f, v)
                    self.config.applied_preset = preset
            return
        # UI toggle
        if command == 'ui:toggle_small':
            st.session_state.is_small_screen = not st.session_state.is_small_screen
            return
        # Layout save
        if command == 'layout:save':
            Path('cache').mkdir(exist_ok=True)
            snap = {
                'version': self.LAYOUT_SNAPSHOT_VERSION,
                'selected_interactions': st.session_state.selected_interactions,
                'individual_strength_filters': st.session_state.individual_strength_filters
            }
            fname = Path('cache')/f"layout_{int(time.time())}.json"
            fname.write_text(json.dumps(snap, indent=2))
            return
        # Layout restore last
        if command == 'layout:restore_last':
            snaps = sorted(Path('cache').glob('layout_*.json'))
            if snaps:
                try:
                    payload = json.loads(snaps[-1].read_text())
                    # Version guard (future migrations possible)
                    version = payload.get('version', 0)
                    if version > self.LAYOUT_SNAPSHOT_VERSION:
                        st.warning(f"Snapshot version {version} is newer than supported ({self.LAYOUT_SNAPSHOT_VERSION}); skipping restore.")
                        return
                    sel = payload.get('selected_interactions'); strengths = payload.get('individual_strength_filters')
                    if isinstance(sel, dict):
                        for k,v in sel.items():
                            if k in st.session_state.selected_interactions:
                                st.session_state.selected_interactions[k] = v
                    if isinstance(strengths, dict):
                        for k,v in strengths.items():
                            st.session_state.individual_strength_filters[k] = v
                except Exception:
                    pass
            return
        # Navigation placeholder (Streamlit cannot programmatically switch tabs easily)
        if command.startswith('nav:tab:'):
            st.session_state.desired_tab = command.split(':')[-1]
            return
    
    def _render_single_input(self):
        """Render single PDB ID input interface."""
        
        # Create a form for PDB input and submission
        with st.form(key="pdb_analysis_form", clear_on_submit=False):
            pdb_id = st.text_input(
                "PDB ID:",
                placeholder="e.g., 1A2B",
                help="Enter a 4-character PDB identifier and press Enter to analyze"
            ).upper()
            
            # Submit button (Enter key will trigger this)
            submitted = st.form_submit_button(
                "🚀 Analyze Structure", 
                type="primary", 
                width="stretch",
                disabled=(sum(1 for selected in st.session_state.selected_interactions.values() if selected) == 0)
            )
            
            if submitted:
                if pdb_id and len(pdb_id) == 4:
                    selected_count = sum(1 for selected in st.session_state.selected_interactions.values() if selected)
                    if selected_count > 0:
                        # Filter selected interactions
                        selected_interaction_types = [
                            itype for itype, selected in st.session_state.selected_interactions.items() 
                            if selected
                        ]
                        self._run_single_analysis(pdb_id, selected_interaction_types)
                    else:
                        st.error("Please select at least one interaction type to analyze")
                else:
                    st.error("Please enter a valid 4-character PDB ID")
    
    def _render_batch_input(self):
        """Render batch upload interface."""
        st.write("Upload a CSV/Excel file with PDB IDs:")
        
        uploaded_file = st.file_uploader(
            "Choose file",
            type=['csv', 'xlsx', 'xls'],
            help="File should contain PDB IDs in the first column"
        )
        
        if uploaded_file:
            try:
                # Read the file
                if uploaded_file.name.endswith('.csv'):
                    df = pd.read_csv(uploaded_file)
                else:
                    df = pd.read_excel(uploaded_file)
                
                # Extract PDB IDs from first column
                pdb_ids = df.iloc[:, 0].astype(str).str.upper().str.strip().tolist()
                pdb_ids = [pid for pid in pdb_ids if len(pid) == 4]
                
                st.write(f"Found {len(pdb_ids)} valid PDB IDs")
                
                if len(pdb_ids) > self.config.processing.max_batch_size:
                    st.warning(f"Batch size limited to {self.config.processing.max_batch_size}. Only the first {self.config.processing.max_batch_size} structures will be processed.")
                    pdb_ids = pdb_ids[:self.config.processing.max_batch_size]
                
                if pdb_ids:
                    st.dataframe(pd.DataFrame({"PDB IDs": pdb_ids[:10]}))  # Show first 10
                    if len(pdb_ids) > 10:
                        st.info(f"... and {len(pdb_ids) - 10} more")
                    
                    # Interaction selection for batch (reuse the same session state)
                    st.subheader("🔬 Select Interactions for Batch Analysis")
                    
                    # Ensure selected_interactions exists
                    if 'selected_interactions' not in st.session_state:
                        interaction_types = get_interaction_types()
                        st.session_state.selected_interactions = {itype: True for itype in interaction_types}
                    
                    # Show current selection
                    selected_count = sum(1 for selected in st.session_state.selected_interactions.values() if selected)
                    display_names = get_interaction_display_names()
                    selected_names = [
                        display_names[itype] for itype, selected in st.session_state.selected_interactions.items()
                        if selected
                    ]
                    
                    if selected_count > 0:
                        st.info(f"✅ Will analyze {selected_count} interaction types: {', '.join(selected_names)}")
                        st.caption("💡 Modify interaction selection in the sidebar")
                    else:
                        st.warning("⚠️ No interactions selected! Please select interactions in the sidebar first.")
                    
                    if st.button("🚀 Analyze Batch", type="primary", width="stretch", disabled=(selected_count == 0)):
                        # Filter selected interactions
                        selected_interaction_types = [
                            itype for itype, selected in st.session_state.selected_interactions.items() 
                            if selected
                        ]
                        self._run_batch_analysis(pdb_ids, selected_interaction_types)
                
            except Exception as e:
                st.error(f"Error reading file: {e}")
    
    def _render_session_loader(self):
        """Render session loading interface."""
        sessions = st.session_state.session_manager.list_sessions()
        
        if sessions:
            session_options = [
                f"{s['name']} ({s['session_id'][:8]}...)" 
                for s in sessions
            ]
            
            selected_idx = st.selectbox(
                "Select session:",
                range(len(session_options)),
                format_func=lambda x: session_options[x]
            )
            
            if st.button("Load Session", width="stretch"):
                session_id = sessions[selected_idx]['session_id']
                if st.session_state.session_manager.load_session(session_id):
                    st.success("Session loaded!")
                    st.rerun()
        else:
            st.info("No saved sessions found")
    
    def _run_single_analysis(self, pdb_id: str, selected_interactions: List[str] = None):
        """Run analysis for a single PDB structure with selected interaction types."""
        # If no interactions are specified, get the currently selected ones from session state
        if selected_interactions is None:
            selected_interactions = [
                itype for itype, selected in st.session_state.selected_interactions.items() 
                if selected
            ]
            
        # If still no interactions selected, show error
        if not selected_interactions:
            st.error("No interaction types selected for analysis")
            return
        progress_bar = st.progress(0)
        status_text = st.empty()
        live_status = st.empty()

        try:
            # Step 1: Download/load structure
            status_text.text(f"Loading structure {pdb_id}...")
            progress_bar.progress(0.2)
            
            structure_settings = st.session_state.structure_settings
            structure = self.pdb_handler.fetch_structure_variant(
                pdb_id,
                assembly=structure_settings.get("assembly", "biological"),
                include_ligands=structure_settings.get("include_ligands", True),
                exclude_waters=structure_settings.get("exclude_waters", True)
            )
            # Cache structure for downstream Ramachandran computation
            st.session_state.rama_source_structure = structure
            
            if not structure:
                st.error(f"❌ Failed to load structure {pdb_id}")
                st.warning("""
                **Possible reasons:**
                - The PDB server is experiencing high load or timeouts
                - The PDB ID may be incorrect
                - Network connectivity issues
                
                **Try:**
                - Wait a moment and try again
                - Verify the PDB ID at [RCSB PDB](https://www.rcsb.org/)
                - Check your internet connection
                """)
                return
            
            # Step 2: Run high-performance analysis with selected interactions
            selected_display_names = [get_interaction_display_names()[itype] for itype in selected_interactions]
            status_text.text(f"🚀 Analyzing {len(selected_interactions)} interaction types: {', '.join(selected_display_names[:3])}{'...' if len(selected_interactions) > 3 else ''}")
            progress_bar.progress(0.5)
            
            # Show performance mode and selected interactions
            perf_info = st.info(f"⚡ Processing Mode: {st.session_state.processing_mode.replace('_', ' ').title()}")
            interaction_info = st.info(f"🔬 Analyzing: {', '.join(selected_display_names)}")
            
            # Track analysis start time
            analysis_start = time.time()
            
            # Use enhanced processor for single protein with selected interactions.
            # Provide structure object directly if processor supports it; fallback to pdb_id.
            # Initialize progress stream container
            # Initialize or append to persistent detector progress log (multi-run)
            if 'detector_progress_panel' not in st.session_state:
                st.session_state.detector_progress_panel = []
            run_id = str(uuid.uuid4())
            run_start_ts = time.time()
            def _progress_cb(ev):
                ev.setdefault('timestamp', time.time())
                ev.setdefault('run_id', run_id)
                ev.setdefault('pdb_id', pdb_id)
                ev.setdefault('run_start', run_start_ts)
                st.session_state.detector_progress_panel.append(ev)
                if 'detector' in ev:
                    live_status.markdown(f"<div style='background:#00B4D8;padding:4px 10px;border-radius:12px;display:inline-block;color:#000;font-size:0.7rem;font-weight:600;'>Running: {ev['detector']} · {ev.get('count','?')}</div>", unsafe_allow_html=True)
            try:
                results = self.batch_processor.process_single_protein(pdb_id, interaction_filters=selected_interactions, structure=structure, progress_callback=_progress_cb)
            except TypeError:
                # Older signature without structure parameter (no progress streaming)
                results = self.batch_processor.process_single_protein(pdb_id, interaction_filters=selected_interactions)
            
            analysis_time = time.time() - analysis_start
            
            progress_bar.progress(0.8)
            
            # Store results with performance metrics
            st.session_state.analysis_results[pdb_id] = results
            st.session_state.current_pdb = pdb_id
            
            # Track performance metrics
            performance_data = {
                'pdb_id': pdb_id,
                'analysis_time': analysis_time,
                'total_interactions': results.get('summary', {}).get('total_interactions', 0),
                'timestamp': time.time(),
                'performance_metrics': results.get('summary', {}).get('performance_metrics', {})
            }
            st.session_state.performance_metrics.append(performance_data)
            
            # Update session metadata
            st.session_state.session_manager.update_metadata(
                total_structures=len(st.session_state.analysis_results),
                total_interactions=sum(
                    r.get('summary', {}).get('total_interactions', 0)
                    for r in st.session_state.analysis_results.values()
                )
            )
            
            progress_bar.progress(1.0)
            status_text.text("✅ High-performance analysis complete!")
            
            # Show performance summary
            if results.get('success', False):
                total_interactions = results.get('summary', {}).get('total_interactions', 0)
                perf_metrics = results.get('summary', {}).get('performance_metrics', {})
                
                success_msg = f"🎉 Successfully analyzed {pdb_id}!\n"
                success_msg += f"📊 Found {total_interactions} interactions in {analysis_time:.2f}s"
                
                if 'average_throughput' in perf_metrics:
                    success_msg += f" ({perf_metrics['average_throughput']} tasks/sec)"
                
                st.success(success_msg)

                # Display applied structural variant settings
                with st.expander("🧬 Structure Variant Settings Applied"):
                    assembly_txt = "Biological Assembly" if structure_settings.get("assembly") == "biological" else "Asymmetric Unit"
                    include_ligands_txt = "Yes" if structure_settings.get("include_ligands") else "No (ligands removed)"
                    exclude_waters_txt = "Yes (waters removed)" if structure_settings.get("exclude_waters") else "No (waters retained)"
                    st.markdown(
                        f"Assembly: **{assembly_txt}**  \n"
                        f"Include Ligands: **{include_ligands_txt}**  \n"
                        f"Exclude Waters: **{exclude_waters_txt}**"
                    )
                
                # Show performance details in expander
                with st.expander("⚡ Performance Details"):
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Analysis Time", f"{analysis_time:.2f}s")
                        
                    with col2:
                        st.metric("Total Interactions", total_interactions)
                        
                    with col3:
                        throughput = total_interactions / analysis_time if analysis_time > 0 else 0
                        st.metric("Interactions/sec", f"{throughput:.1f}")
                    
                    if perf_metrics:
                        st.markdown("**Performance Metrics:**")
                        # Render key metrics as simple bullet points (plain text style)
                        for k, v in perf_metrics.items():
                            st.markdown(f"- {k.replace('_',' ').title()}: {v}")
            else:
                st.error(f"❌ Analysis failed for {pdb_id}: {results.get('error', 'Unknown error')}")
            
            perf_info.empty()  # Clear the performance info
            
        except Exception as e:
            st.error(f"Analysis failed: {e}")
        finally:
            progress_bar.empty(); status_text.empty(); live_status.empty()
            # Autosave minimal session meta
            try:
                import json, pathlib
                pathlib.Path('cache').mkdir(parents=True, exist_ok=True)
                meta = {
                    'selected_interactions': st.session_state.selected_interactions,
                    'analysis_result_ids': list(st.session_state.analysis_results.keys())
                }
                (pathlib.Path('cache')/ 'session_state.json').write_text(json.dumps(meta, indent=2))
            except Exception:
                pass
    
    def _run_batch_analysis(self, pdb_ids: List[str], selected_interactions: List[str] = None):
        """Run high-performance batch analysis for multiple PDB structures with selected interactions."""
        if selected_interactions is None:
            selected_interactions = get_interaction_types()
            
        batch_id = str(uuid.uuid4())
        
        # Create progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
        # Skeleton loader container per selected interaction (visual placeholder)
        skel_container = st.empty()
        try:
            with skel_container.container():
                if selected_interactions:
                    cols = st.columns(min(4, len(selected_interactions)))
                    for idx, itype in enumerate(selected_interactions[:8]):
                        with cols[idx % len(cols)]:
                            st.markdown(f"<div style='padding:10px;border:1px solid #444;border-radius:8px;background:linear-gradient(90deg,#2a2b2d,#333);animation:pulse 1.2s ease-in-out infinite;'>⏳ {itype}</div>", unsafe_allow_html=True)
                st.markdown("<style>@keyframes pulse{0%{opacity:.55}50%{opacity:1}100%{opacity:.55}}</style>", unsafe_allow_html=True)
        except Exception:
            pass
        results_container = st.empty()
        
        # Show batch processing info with interaction details
        selected_display_names = [get_interaction_display_names()[itype] for itype in selected_interactions]
        batch_info = st.info(f"🚀 High-Performance Batch Mode: Processing {len(pdb_ids)} proteins with {self.batch_processor.processor.max_workers} parallel workers")
        interaction_info = st.info(f"🔬 Analyzing {len(selected_interactions)} interaction types: {', '.join(selected_display_names)}")
        
        try:
            # Track batch analysis start time
            batch_start = time.time()
            
            # Run high-performance batch processing with selected interactions
            status_text.text(f"🚀 Starting parallel batch analysis of {len(pdb_ids)} proteins for {len(selected_interactions)} interaction types...")
            
            results = self.batch_processor.process_multiple_proteins(pdb_ids, interaction_filters=selected_interactions)
            
            batch_time = time.time() - batch_start
            
            # Process and store results
            successful_results = [r for r in results if r.get('success', False)]
            failed_results = [r for r in results if not r.get('success', False)]
            
            # Store results in session
            for result in successful_results:
                pdb_id = result['pdb_id']
                st.session_state.analysis_results[pdb_id] = result
            
            # Track batch performance
            batch_performance = {
                'batch_id': batch_id,
                'total_proteins': len(pdb_ids),
                'successful': len(successful_results),
                'failed': len(failed_results),
                'batch_time': batch_time,
                'throughput': len(successful_results) / batch_time if batch_time > 0 else 0,
                'timestamp': time.time(),
                'parallel_workers': self.batch_processor.processor.max_workers,
                'analyzed_interactions': selected_interactions
            }
            st.session_state.performance_metrics.append(batch_performance)
            
            # Update session metadata
            st.session_state.session_manager.update_metadata(
                total_structures=len(st.session_state.analysis_results),
                total_interactions=sum(
                    r.get('summary', {}).get('total_interactions', 0)
                    for r in st.session_state.analysis_results.values()
                )
            )
            
            progress_bar.progress(1.0)
            status_text.text("✅ High-performance batch analysis complete!")
            
            # Show comprehensive results
            if successful_results:
                total_interactions = sum(r.get('summary', {}).get('total_interactions', 0) for r in successful_results)
                
                success_msg = f"🎉 Batch analysis complete!\n"
                success_msg += f"✅ Successfully processed {len(successful_results)}/{len(pdb_ids)} proteins\n"
                success_msg += f"📊 Found {total_interactions} total interactions in {batch_time:.2f}s\n"
                success_msg += f"⚡ Throughput: {batch_performance['throughput']:.2f} proteins/sec"
                
                st.success(success_msg)
                
                # Show detailed batch performance
                with st.expander("🚀 Batch Performance Analysis"):
                    col1, col2, col3, col4 = st.columns(4)
                    
                    with col1:
                        st.metric("Batch Time", f"{batch_time:.2f}s")
                        
                    with col2:
                        st.metric("Success Rate", f"{len(successful_results)/len(pdb_ids)*100:.1f}%")
                        
                    with col3:
                        st.metric("Parallel Workers", batch_performance['parallel_workers'])
                        
                    with col4:
                        st.metric("Avg per Protein", f"{batch_time/len(pdb_ids):.2f}s")
                    
                    # Show per-protein breakdown
                    if st.checkbox("Show per-protein timing"):
                        perf_df = pd.DataFrame([
                            {
                                'PDB ID': r['pdb_id'],
                                'Success': '✅' if r.get('success', False) else '❌',
                                'Processing Time': f"{r.get('processing_time', 0):.2f}s",
                                'Total Interactions': r.get('summary', {}).get('total_interactions', 0)
                            }
                            for r in results
                        ])
                        st.dataframe(perf_df, width="stretch")
                
                # Show failed analyses if any
                if failed_results:
                    with st.expander(f"❌ Failed Analyses ({len(failed_results)})"):
                        for result in failed_results:
                            st.error(f"{result['pdb_id']}: {result.get('error', 'Unknown error')}")
            else:
                st.error("❌ All batch analyses failed!")
            
            batch_info.empty()  # Clear batch info
            
            # Offer download of results
            if st.button("📥 Download Batch Results"):
                zip_data = self._create_batch_download(results)
                st.download_button(
                    label="Download ZIP",
                    data=zip_data,
                    file_name=f"batch_results_{batch_id[:8]}.zip",
                    mime="application/zip"
                )
                
        except Exception as e:
            st.error(f"Batch analysis failed: {e}")
        finally:
            progress_bar.empty()
            status_text.empty()
            skel_container.empty()
    
    def _update_batch_progress(self, progress_bar, status_text, results_container, 
                             current: int, total: int, pdb_id: str):
        """Update batch processing progress."""
        progress = current / total
        progress_bar.progress(progress)
        status_text.text(f"Analyzing {pdb_id}... ({current}/{total})")
        
        # Show intermediate results
        with results_container.container():
            if current > 0:
                st.metric("Structures Processed", f"{current}/{total}")
    
    def _create_batch_download(self, results: Dict[str, Any]) -> bytes:
        """Create downloadable ZIP file with batch results."""
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for pdb_id, result in results.items():
                # Add CSV report for each structure
                csv_data = self.report_generator.generate_csv_report(pdb_id, result)
                zip_file.writestr(f"{pdb_id}_interactions.csv", csv_data)
                
                # Add summary JSON
                import json
                json_data = json.dumps(result, indent=2, default=str)
                zip_file.writestr(f"{pdb_id}_summary.json", json_data)
        
        return zip_buffer.getvalue()
    
    def _render_analysis_tab(self):
        """Render the analysis configuration tab with comprehensive controls."""
        st.header("🔍 Interaction Analysis Configuration")
        # Parameter diff badge (active preset vs current config values)
        try:
            applied = getattr(self.config, 'applied_preset', 'literature_default')
            preset_vals = self.config.presets.get(applied, {})
            diffs = []
            for k, v in preset_vals.items():
                cur = getattr(self.config.interactions, k, None)
                if cur is None:
                    continue
                if isinstance(cur, (int, float)) and isinstance(v, (int, float)):
                    if abs(float(cur) - float(v)) > 1e-9:
                        diffs.append((k, v, cur))
                else:
                    if cur != v:
                        diffs.append((k, v, cur))
            badge_color = "orange" if diffs else "green"
            st.markdown(f"<div style='padding:6px 10px;border-radius:6px;background:rgba(0,180,216,0.08);border:1px solid #00B4D8;display:inline-block;margin-bottom:0.5rem;'>Preset: <strong>{applied}</strong> · Status: <span style='color:{'var(--warning)' if diffs else '#3adb76'}'>{'Custom' if diffs else 'In Sync'}</span></div>", unsafe_allow_html=True)
            if diffs:
                with st.expander(f"View parameter differences ({len(diffs)})", expanded=False):
                    for k, base, cur in diffs[:120]:
                        st.write(f"{k}: {base} → {cur}")
                    if len(diffs) > 120:
                        st.caption("(Truncated list)")
        except Exception:
            pass
        
        if not st.session_state.analysis_results:
            st.info("👆 Use the sidebar to input structures for analysis")
            self._render_single_input()
            return
        
        # Ensure a PDB is selected if results exist, to prevent KeyErrors
        if not st.session_state.current_pdb or st.session_state.current_pdb not in st.session_state.analysis_results:
            st.session_state.current_pdb = list(st.session_state.analysis_results.keys())[0]
            st.rerun()

        # Display current results summary
        self._render_analysis_summary()

        # (Removed) Duplicate parameter sliders moved exclusively to sidebar to reduce redundancy.
        # If needed later for large-screen workflows, we can re-enable with a feature flag.

        # Detailed point-wise analysis summary
        st.write("---")
        self._render_pointwise_analysis_summary()

        # Performance profiling panel (if timings available for current structure)
        cur_result = st.session_state.analysis_results.get(st.session_state.current_pdb, {}) if st.session_state.current_pdb else {}
        timings = cur_result.get('summary', {}).get('detector_timings', {})
        counts = cur_result.get('summary', {}).get('detector_counts', {})
        if timings:
            with st.expander("⏱️ Detector Performance", expanded=False):
                try:
                    import pandas as pd  # type: ignore
                except Exception:
                    pd = None  # fallback
                rows = []
                for name, t in sorted(timings.items(), key=lambda x: x[1]):
                    rows.append({
                        'Detector': name,
                        'Time (s)': round(t, 4),
                        'Count': counts.get(name, 0),
                        'Rate (int/s)': round(counts.get(name,0)/t,2) if t>0 and counts.get(name,0)>0 else 0
                    })
                if rows and pd is not None:
                    dfp = pd.DataFrame(rows)
                    st.dataframe(dfp, use_container_width=True, height=min(400, 35*len(rows)+60))
                else:
                    for r in rows:
                        st.write(f"{r['Detector']}: {r['Time (s)']} s | {r['Count']} interactions | {r['Rate (int/s)']} int/s")
                total_time = sum(timings.values())
                st.progress(0.0)
                # Mini horizontal bar summary
                for name, t in sorted(timings.items(), key=lambda x: x[1], reverse=True):
                    pct = (t/total_time)*100 if total_time>0 else 0
                    st.write(f"{name} — {round(t,4)} s ({pct:.1f}%)")
                st.caption("Timings gathered via lightweight re-run per detector (approximate; excludes parallel scheduling overhead).")
                # Aggregate across session
                agg = self._compute_aggregated_performance()
                if agg:
                    st.markdown("**Session Aggregate**")
                    try:
                        import pandas as pd  # type: ignore
                    except Exception:
                        pd = None
                    if pd:
                        agg_df = pd.DataFrame(agg['rows'])
                        st.dataframe(agg_df, use_container_width=True, height=min(400, 40*len(agg['rows'])+60))
                    else:
                        for r in agg['rows']:
                            st.write(f"{r['Detector']}: {r['Total Time (s)']} s total | {r['Avg per Structure (s)']} s/structure | {r['Total Count']} interactions | {r['Time %']:.1f}%")
                    st.caption(f"Structures with timing data: {agg['structures']} | Cumulative detector time: {round(agg['grand_total'],3)} s")

        # Global parameter lock
        if 'parameters_locked' not in st.session_state:
            st.session_state.parameters_locked = False
        lock_col1, lock_col2 = st.columns([1,5])
        with lock_col1:
            locked = st.toggle("🔒 Lock", value=st.session_state.parameters_locked, help="Prevent accidental slider changes")
            st.session_state.parameters_locked = locked
        if locked:
            st.caption("Parameters locked. Unlock to edit in sidebar or settings tab.")

        
        # Parameter signature + re-analysis controls
        st.write("---")

        # Build a stable parameter signature hash
        import hashlib, json, datetime
        param_source = st.session_state.get('custom_parameters', {})
        # Include selected interactions & strength filters so those changes also trigger enabled state
        signature_payload = {
            'params': {k: param_source.get(k) for k in sorted(param_source.keys())},
            'selected_interactions': {k: v for k, v in sorted(st.session_state.selected_interactions.items())},
            'strength_filters': st.session_state.get('individual_strength_filters', {}),
        }
        signature_str = json.dumps(signature_payload, sort_keys=True, default=str)
        signature_hash = hashlib.sha256(signature_str.encode()).hexdigest()[:16]

        if 'last_analysis_signature' not in st.session_state:
            st.session_state.last_analysis_signature = None
        if 'last_analysis_timestamp' not in st.session_state:
            st.session_state.last_analysis_timestamp = None

        changed = signature_hash != st.session_state.last_analysis_signature

        # UI layout
        c1, c2, c3 = st.columns([2,2,3])
        with c1:
            st.markdown(f"**Parameter Signature:** `{signature_hash}`")
        with c2:
            if st.session_state.last_analysis_timestamp:
                ts = st.session_state.last_analysis_timestamp
                st.markdown(f"🕒 Last analyzed: `{ts}`")
            else:
                st.markdown("🕒 Last analyzed: _never_")
        with c3:
            if st.session_state.last_analysis_signature is None:
                # First run: neutral guidance, don't imply change
                st.caption("Adjust parameters then run analysis when ready.")
            else:
                if changed:
                    st.success("Changes detected — re-analysis needed to update results.")
                else:
                    # Avoid noisy message; show subtle caption instead of info box
                    st.caption("Parameters unchanged since last analysis.")

        st.write("")
        btn_disabled = not changed
        if st.button("🔄 Re-analyze with Current Settings", type="primary", disabled=btn_disabled, width='stretch'):
            if st.session_state.current_pdb:
                selected_interactions = [itype for itype, selected in st.session_state.selected_interactions.items() if selected]
                if selected_interactions:
                    with st.spinner("Re-analyzing with new parameters..."):
                        self._run_single_analysis(st.session_state.current_pdb, selected_interactions)
                    # Persist new signature + timestamp
                    st.session_state.last_analysis_signature = signature_hash
                    st.session_state.last_analysis_timestamp = datetime.datetime.now().strftime('%H:%M:%S')
                    st.success("✅ Re-analysis complete!")
                else:
                    st.error("❌ Please select at least one interaction type")
            else:
                st.error("❌ No structure selected")

        # Quick per-interaction re-analysis buttons (selected interactions only)
        selected_now = [itype for itype, sel in st.session_state.selected_interactions.items() if sel]
        if selected_now:
            with st.expander("🔁 Re-run Specific Interaction Types", expanded=False):
                cols_int = st.columns(min(4, len(selected_now)))
                for idx, itype in enumerate(selected_now):
                    col = cols_int[idx % len(cols_int)]
                    with col:
                        if st.button(f"Re-run {itype}", key=f"rerun_{itype}"):
                            # Run single-detector style analysis by re-calling processor with one filter
                            if st.session_state.current_pdb:
                                try:
                                    res_partial = self.batch_processor.process_single_protein(st.session_state.current_pdb, interaction_filters=[itype])
                                    # Merge updated interactions for that type only
                                    existing = st.session_state.analysis_results.get(st.session_state.current_pdb, {})
                                    if existing and 'interactions' in existing and 'interactions' in res_partial:
                                        existing['interactions'][itype] = res_partial['interactions'].get(itype, [])
                                        # Update counts & timings selectively
                                        if 'summary' in existing and 'summary' in res_partial:
                                            existing['summary']['interaction_counts'][itype] = len(existing['interactions'][itype])
                                            # Recompute total interactions
                                            existing['summary']['total_interactions'] = sum(existing['summary']['interaction_counts'].values())
                                            # Merge timing/counts for detector
                                            pt = res_partial.get('summary', {}).get('detector_timings', {})
                                            pc = res_partial.get('summary', {}).get('detector_counts', {})
                                            if pt:
                                                existing['summary'].setdefault('detector_timings', {}).update(pt)
                                            if pc:
                                                existing['summary'].setdefault('detector_counts', {}).update(pc)
                                    st.success(f"Re-ran {itype}")
                                except Exception as e:
                                    st.error(f"Failed to re-run {itype}: {e}")

        # Session snapshot export / import
        st.write("---")
        with st.expander("💾 Session Snapshot Export / Import", expanded=False):
            import json, datetime
            # Export snapshot
            if st.button("Export Snapshot", key="btn_export_snapshot"):
                try:
                    # Extract interaction parameters dynamically from config.interactions (public attributes not starting with _)
                    param_obj = getattr(self.config, 'interactions', None)
                    params = {}
                    if param_obj:
                        for attr in dir(param_obj):
                            if attr.startswith('_'):
                                continue
                            val = getattr(param_obj, attr)
                            if isinstance(val, (int, float, str, bool)):
                                params[attr] = val
                    snapshot = {
                        'timestamp': datetime.datetime.utcnow().isoformat()+'Z',
                        'applied_preset': getattr(self.config, 'applied_preset', None),
                        'interaction_parameters': params,
                        'selected_interactions': st.session_state.get('selected_interactions'),
                        'strength_filters': st.session_state.get('individual_strength_filters'),
                        'structures_analyzed': list(st.session_state.analysis_results.keys()),
                        'interaction_counts_summary': {pid: res.get('summary', {}).get('interaction_counts', {}) for pid, res in st.session_state.analysis_results.items()},
                        'detector_timings_summary': {pid: res.get('summary', {}).get('detector_timings', {}) for pid, res in st.session_state.analysis_results.items()},
                    }
                    st.download_button(
                        label="Download Snapshot JSON",
                        data=json.dumps(snapshot, indent=2),
                        file_name="interaction_session_snapshot.json",
                        mime="application/json",
                        key="download_snapshot_btn"
                    )
                except Exception as e:
                    st.error(f"Failed to build snapshot: {e}")
            uploaded = st.file_uploader("Import Snapshot JSON", type=['json'], key="snapshot_import_uploader")
            if uploaded is not None:
                try:
                    raw = uploaded.read().decode('utf-8')
                    data = json.loads(raw)
                    # Restore interaction parameters
                    params = data.get('interaction_parameters', {})
                    param_obj = getattr(self.config, 'interactions', None)
                    if param_obj and isinstance(params, dict):
                        for k, v in params.items():
                            if hasattr(param_obj, k):
                                try:
                                    setattr(param_obj, k, v)
                                except Exception:
                                    pass
                    # Restore selections & filters
                    sel = data.get('selected_interactions')
                    if isinstance(sel, dict):
                        st.session_state.selected_interactions.update(sel)
                    sfilters = data.get('strength_filters')
                    if isinstance(sfilters, dict):
                        st.session_state.individual_strength_filters.update(sfilters)
                    if data.get('applied_preset'):
                        self.config.applied_preset = data['applied_preset']
                    st.success("Snapshot imported. Re-run analysis for detailed interactions.")
                except Exception as e:
                    st.error(f"Failed to import snapshot: {e}")
    
    def _render_analysis_summary(self):
        """Render summary of current analysis results."""
        total_structures = len(st.session_state.analysis_results)
        total_interactions = sum(
            len(result.get('interactions', {})) 
            for result in st.session_state.analysis_results.values()
        )
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Structures Analyzed", total_structures)
        
        with col2:
            st.metric("Total Interactions", total_interactions)
        
        with col3:
            avg_interactions = total_interactions / total_structures if total_structures > 0 else 0
            st.metric("Avg. Interactions/Structure", f"{avg_interactions:.1f}")
        
        # Structure selector
        if total_structures > 1:
            structure_options = list(st.session_state.analysis_results.keys())
            selected_structure = st.selectbox(
                "Select structure to view:",
                structure_options,
                index=structure_options.index(st.session_state.current_pdb) 
                if st.session_state.current_pdb in structure_options else 0
            )
            st.session_state.current_pdb = selected_structure
    
    def _render_parameter_controls(self):
        """Render comprehensive parameter control sliders for all interactions."""
        st.write("🎛️ **Adjust detection parameters for precise interaction identification**")
        # Auto-bind quick preset (from sidebar) once per change
        quick_preset = st.session_state.structure_settings.get("preset") if 'structure_settings' in st.session_state else None
        preset_map = {
            'conservative': {
                'strong_threshold': 0.9, 'moderate_threshold': 0.75,
                'hbond_distance': 3.2, 'halogen_distance': 3.8,
                'pi_pi_distance': 5.0, 'ionic_distance': 5.5
            },
            'literature_default': {},  # uses central defaults
            'exploratory': {
                'strong_threshold': 0.7, 'moderate_threshold': 0.5,
                'hbond_distance': 4.0, 'halogen_distance': 4.5,
                'pi_pi_distance': 6.0, 'ionic_distance': 7.0
            }
        }
        if quick_preset and 'last_quick_preset' not in st.session_state:
            st.session_state.last_quick_preset = None
        if quick_preset and quick_preset != st.session_state.get('last_quick_preset'):
            # Apply mapped values WITHOUT wiping other parameters
            for k, v in preset_map.get(quick_preset, {}).items():
                st.session_state[k] = v
            st.session_state.last_quick_preset = quick_preset
            st.session_state.last_preset_applied = quick_preset.replace('_',' ').title()
            st.caption(f"⚡ Quick Preset Applied: {st.session_state.last_preset_applied}")
        # Central defaults so we can reliably reset / apply presets (must be declared first)
        parameter_defaults = {
            'hbond_distance': 3.5, 'hbond_angle': 120,
            'ionic_distance': 6.0,
            'hydrophobic_distance': 5.0,
            'dispersion_min': 3.5, 'dispersion_max': 5.0,
            'halogen_distance': 4.0, 'halogen_angle': 140,
            'chalcogen_distance': 4.0, 'chalcogen_angle': 140,
            'pnictogen_distance': 4.0, 'pnictogen_angle': 140,
            'tetrel_distance_min': 2.5, 'tetrel_distance_max': 3.6, 'tetrel_angle': 160,
            'pi_pi_distance': 5.5, 'pi_pi_angle': 30,
            'ch_pi_distance': 4.5, 'ch_pi_angle': 90,
            'anion_pi_distance': 5.0,
            'n_pi_star_distance': 3.5, 'n_pi_star_angle': 120,
            'strong_threshold': 0.85, 'moderate_threshold': 0.65,
            'exclude_intraresidue': False, 'min_sequence_separation': 1,
        }

        # Ensure container for custom parameters exists
        if 'custom_parameters' not in st.session_state:
            st.session_state.custom_parameters = {}

        # Initialize any missing defaults BEFORE widgets are instantiated
        for k, v in parameter_defaults.items():
            if k not in st.session_state:
                st.session_state[k] = v

        # Helper used by button callbacks (executed before sliders when placed first)
        def _apply_preset(values: dict, reset: bool = False, label: str | None = None):
            if reset:
                for dk, dv in parameter_defaults.items():
                    st.session_state[dk] = dv
            for pk, pv in values.items():
                st.session_state[pk] = pv
            if label:
                st.session_state['last_preset_applied'] = label

        # Preset controls FIRST so their callbacks run before slider widgets are created
        preset_col1, preset_col2, preset_col3 = st.columns(3)
        with preset_col1:
            st.button("🔄 Reset Defaults", key="btn_reset_defaults", help="Reset all parameters to default values", on_click=_apply_preset, kwargs={'values': {}, 'reset': True, 'label': 'Defaults'})
        with preset_col2:
            st.button("📊 Conservative", key="btn_conservative", help="High-confidence stricter cutoffs", on_click=_apply_preset, kwargs={'values': {
                'strong_threshold': 0.9, 'moderate_threshold': 0.75,
                'hbond_distance': 3.2, 'halogen_distance': 3.8,
                'pi_pi_distance': 5.0, 'ionic_distance': 5.5
            }, 'reset': False, 'label': 'Conservative'})
        with preset_col3:
            st.button("🔍 Exploratory", key="btn_exploratory", help="Liberal discovery-oriented cutoffs", on_click=_apply_preset, kwargs={'values': {
                'strong_threshold': 0.7, 'moderate_threshold': 0.5,
                'hbond_distance': 4.0, 'halogen_distance': 4.5,
                'pi_pi_distance': 6.0, 'ionic_distance': 7.0
            }, 'reset': False, 'label': 'Exploratory'})

        if st.session_state.get('last_preset_applied'):
            st.caption(f"Preset applied: {st.session_state['last_preset_applied']}. Remember to click 'Re-analyze' to recompute interactions.")

        # Customization badge – compare current params to defaults or last applied preset values
        tracked_keys = list(parameter_defaults.keys())
        current_values = {k: st.session_state.get(k) for k in tracked_keys}
        baseline = parameter_defaults.copy()
        # If a quick preset or button preset applied, fold in its modifications to baseline for comparison
        lp = st.session_state.get('last_preset_applied')
        reverse_lookup = {
            'Conservative': 'conservative',
            'Exploratory': 'exploratory',
            'Literature Default': 'literature_default'
        }
        applied_key = reverse_lookup.get(lp, None)
        if applied_key and applied_key in preset_map:
            for k, v in preset_map[applied_key].items():
                baseline[k] = v
        # Determine differences
        diffs = {k: (baseline[k], current_values[k]) for k in tracked_keys if baseline[k] != current_values[k]}
        if diffs:
            with st.expander(f"✨ Customized Parameters ({len(diffs)})", expanded=False):
                for k, (b, c) in diffs.items():
                    st.write(f"{k}: {b} → {c}")
        else:
            st.caption("✅ Parameters match applied preset/defaults (no custom changes)")

        # --- User-defined preset management -------------------------------------------------
        import json, datetime, copy
        if 'saved_parameter_presets' not in st.session_state:
            st.session_state.saved_parameter_presets = {}
        if 'active_parameter_preset' not in st.session_state:
            st.session_state.active_parameter_preset = None

        with st.expander("💾 Save / Load Custom Presets", expanded=False):
            csp1, csp2, csp3, csp4 = st.columns([2,2,2,2])
            with csp1:
                preset_name = st.text_input("Preset Name", key="preset_name_input", placeholder="e.g. cavity_scan")
            with csp2:
                if st.button("💾 Save Preset", disabled=not preset_name.strip()):
                    snapshot = {k: st.session_state.get(k) for k in parameter_defaults.keys()}
                    st.session_state.saved_parameter_presets[preset_name.strip()] = {
                        'values': snapshot,
                        'created_at': datetime.datetime.now().isoformat(timespec='seconds')
                    }
                    st.session_state.active_parameter_preset = preset_name.strip()
                    st.success(f"Preset '{preset_name.strip()}' saved")
                    st.rerun()
            with csp3:
                existing = list(st.session_state.saved_parameter_presets.keys())
                selected_preset = st.selectbox("Load Preset", options=["(none)"] + existing, index=0 if not existing else (existing.index(st.session_state.active_parameter_preset)+1 if st.session_state.active_parameter_preset in existing else 0), key="preset_loader")
            with csp4:
                if st.button("📂 Apply Loaded", disabled=(selected_preset == "(none)")):
                    if selected_preset != "(none)":
                        data = st.session_state.saved_parameter_presets[selected_preset]['values']
                        for k, v in data.items():
                            st.session_state[k] = v
                        st.session_state.active_parameter_preset = selected_preset
                        st.session_state.last_preset_applied = selected_preset
                        st.info(f"Preset '{selected_preset}' applied. Re-analyze to update results.")
                        st.rerun()

            # Manage presets list
            if st.session_state.saved_parameter_presets:
                with st.expander("Manage Saved Presets", expanded=False):
                    for pname, meta in list(st.session_state.saved_parameter_presets.items()):
                        colsP = st.columns([3,2,1,1])
                        with colsP[0]:
                            st.markdown(f"**{pname}**")
                        with colsP[1]:
                            st.caption(meta.get('created_at',''))
                        with colsP[2]:
                            if st.button("⬇️", key=f"download_preset_{pname}", help="Download preset as JSON"):
                                json_bytes = json.dumps(meta['values'], indent=2).encode()
                                st.download_button(label=f"Download {pname}.json", data=json_bytes, file_name=f"preset_{pname}.json", mime='application/json', key=f"dlbtn_{pname}")
                        with colsP[3]:
                            if st.button("🗑️", key=f"delete_preset_{pname}", help="Delete preset"):
                                del st.session_state.saved_parameter_presets[pname]
                                if st.session_state.active_parameter_preset == pname:
                                    st.session_state.active_parameter_preset = None
                                st.warning(f"Deleted preset '{pname}'")
                                st.rerun()

            # Upload external preset JSON
            uploaded = st.file_uploader("Import Preset JSON", type=['json'], key="preset_upload")
            if uploaded is not None:
                try:
                    payload = json.loads(uploaded.read().decode())
                    if isinstance(payload, dict):
                        imp_name = f"import_{len(st.session_state.saved_parameter_presets)+1}"
                        # Filter only recognized keys
                        filtered = {k: v for k, v in payload.items() if k in parameter_defaults}
                        st.session_state.saved_parameter_presets[imp_name] = {'values': filtered, 'created_at': datetime.datetime.now().isoformat(timespec='seconds')}
                        st.success(f"Imported preset as '{imp_name}'")
                        st.rerun()
                    else:
                        st.error("Invalid preset file structure.")
                except Exception as e:
                    st.error(f"Failed to import preset: {e}")
        
        # Organize parameters by interaction type
        # Compact mode toggle (collapses rarely used groups into expanders instead of tabs)
        if 'parameter_compact_mode' not in st.session_state:
            st.session_state.parameter_compact_mode = False
        compact_bar = st.toggle("Compact Mode (hide advanced groups)", value=st.session_state.parameter_compact_mode, key="compact_mode_toggle")
        st.session_state.parameter_compact_mode = compact_bar

        if not st.session_state.parameter_compact_mode:
            parameter_tabs = st.tabs([
                "🔹 Basic Interactions", 
                "🔸 σ-hole Interactions", 
                "🔺 π-System Interactions",
                "🔻 Other Interactions"
            ])
        else:
            # In compact mode, show only first two groups as tabs; others as inline expanders later
            parameter_tabs = st.tabs([
                "🔹 Basic", 
                "🔸 σ-hole"
            ])
        
        if not st.session_state.parameter_compact_mode:
            with parameter_tabs[2]:
                st.write("**π-System Related Interactions**")
                col1, col2 = st.columns(2)
                
                with col1:
                    # π-π Stacking (render only if selected)
                    if st.session_state.get('selected_interactions', {}).get('pipi', True):
                        st.write("**π-π Stacking**")
                        st.session_state.custom_parameters["pi_pi_distance"] = st.slider(
                            "Distance Cutoff (Å)",
                            min_value=3.0, max_value=8.0, value=st.session_state.get('pi_pi_distance', 5.5), step=0.1,
                            key="pi_pi_distance", help="Maximum centroid-centroid distance"
                        )
                        st.session_state.custom_parameters["pi_pi_angle"] = st.slider(
                            "Angle Cutoff (°)",
                            min_value=0, max_value=45, value=st.session_state.get('pi_pi_angle', 30), step=5,
                            key="pi_pi_angle", help="Maximum dihedral angle between rings"
                        )

                    # C-H···π Interactions (render only if selected)
                    if st.session_state.get('selected_interactions', {}).get('chpi', True):
                        st.write("**C-H···π Interactions**")
                        st.session_state.custom_parameters["ch_pi_distance"] = st.slider(
                            "Distance Cutoff (Å)",
                            min_value=3.0, max_value=6.0, value=st.session_state.get('ch_pi_distance', 4.5), step=0.1,
                            key="ch_pi_distance", help="Maximum C···centroid distance"
                        )
                        st.session_state.custom_parameters["ch_pi_angle"] = st.slider(
                            "Angle Cutoff (°)",
                            min_value=60, max_value=120, value=st.session_state.get('ch_pi_angle', 90), step=5,
                            key="ch_pi_angle", help="Maximum C-H···centroid angle"
                        )
                
                with col2:
                    # Anion-π Interactions (render only if selected)
                    if st.session_state.get('selected_interactions', {}).get('anionpi', True):
                        st.write("**Anion-π Interactions**")
                        st.session_state.custom_parameters["anion_pi_distance"] = st.slider(
                            "Distance Cutoff (Å)",
                            min_value=3.0, max_value=7.0, value=st.session_state.get('anion_pi_distance', 5.0), step=0.1,
                            key="anion_pi_distance", help="Maximum anion-centroid distance"
                        )

                    # n→π* Interactions (render only if selected)
                    if st.session_state.get('selected_interactions', {}).get('npistar', True):
                        st.write("**n→π* Interactions**")
                        st.session_state.custom_parameters["n_pi_star_distance"] = st.slider(
                            "Distance Cutoff (Å)",
                            min_value=2.5, max_value=4.5, value=st.session_state.get('n_pi_star_distance', 3.5), step=0.1,
                            key="n_pi_star_distance", help="Maximum n···C=O distance"
                        )
                        st.session_state.custom_parameters["n_pi_star_angle"] = st.slider(
                            "Angle Cutoff (°)",
                            min_value=90, max_value=150, value=st.session_state.get('n_pi_star_angle', 120), step=5,
                            key="n_pi_star_angle", help="Minimum n···C=O angle"
                        )
        else:
            with st.expander("🔺 π-System Interactions", expanded=False):
                col1, col2 = st.columns(2)
                with col1:
                    if st.session_state.get('selected_interactions', {}).get('pipi', True):
                        st.write("**π-π Stacking**")
                        st.session_state.custom_parameters["pi_pi_distance"] = st.slider(
                            "Distance Cutoff (Å)", min_value=3.0, max_value=8.0, value=st.session_state.get('pi_pi_distance', 5.5), step=0.1, key="pi_pi_distance")
                        st.session_state.custom_parameters["pi_pi_angle"] = st.slider(
                            "Angle Cutoff (°)", min_value=0, max_value=45, value=st.session_state.get('pi_pi_angle', 30), step=5, key="pi_pi_angle")
                    if st.session_state.get('selected_interactions', {}).get('chpi', True):
                        st.write("**C-H···π**")
                        st.session_state.custom_parameters["ch_pi_distance"] = st.slider(
                            "Distance (Å)", min_value=3.0, max_value=6.0, value=st.session_state.get('ch_pi_distance', 4.5), step=0.1, key="ch_pi_distance")
                        st.session_state.custom_parameters["ch_pi_angle"] = st.slider(
                            "Angle (°)", min_value=60, max_value=120, value=st.session_state.get('ch_pi_angle', 90), step=5, key="ch_pi_angle")
                with col2:
                    if st.session_state.get('selected_interactions', {}).get('anionpi', True):
                        st.write("**Anion-π**")
                        st.session_state.custom_parameters["anion_pi_distance"] = st.slider(
                            "Distance (Å)", min_value=3.0, max_value=7.0, value=st.session_state.get('anion_pi_distance', 5.0), step=0.1, key="anion_pi_distance")
                    if st.session_state.get('selected_interactions', {}).get('npistar', True):
                        st.write("**n→π***")
                        st.session_state.custom_parameters["n_pi_star_distance"] = st.slider(
                            "Distance (Å)", min_value=2.5, max_value=4.5, value=st.session_state.get('n_pi_star_distance', 3.5), step=0.1, key="n_pi_star_distance")
                        st.session_state.custom_parameters["n_pi_star_angle"] = st.slider(
                            "Angle (°)", min_value=90, max_value=150, value=st.session_state.get('n_pi_star_angle', 120), step=5, key="n_pi_star_angle")
            st.write("**Basic Non-Covalent Interactions**")
        if not st.session_state.parameter_compact_mode:
            with parameter_tabs[3]:
                st.write("**Specialized Parameters**")
                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Quality Thresholds**")
                    # (single authoritative sliders defined here only)
                    st.session_state.custom_parameters["strong_threshold"] = st.slider(
                        "Strong Interaction Threshold", min_value=0.7, max_value=0.95,
                        value=st.session_state.get('strong_threshold', 0.85), step=0.05,
                        key="strong_threshold", help="Score threshold for strong classification"
                    )
                    st.session_state.custom_parameters["moderate_threshold"] = st.slider(
                        "Moderate Interaction Threshold", min_value=0.5, max_value=0.8,
                        value=st.session_state.get('moderate_threshold', 0.65), step=0.05,
                        key="moderate_threshold", help="Score threshold for moderate classification"
                    )
                with col2:
                    st.write("**Environmental Factors**")
                    st.session_state.custom_parameters["exclude_intraresidue"] = st.checkbox(
                        "Exclude Intraresidue Interactions", value=st.session_state.get('exclude_intraresidue', False),
                        key="exclude_intraresidue", help="Skip interactions within the same residue"
                    )
                    st.session_state.custom_parameters["min_sequence_separation"] = st.slider(
                        "Min Sequence Separation", min_value=0, max_value=5,
                        value=st.session_state.get('min_sequence_separation', 1), step=1,
                        key="min_sequence_separation", help="Minimum residue separation in sequence"
                    )
        else:
            with st.expander("🔻 Other / Specialized", expanded=False):
                col1, col2 = st.columns(2)
                with col1:
                    st.write("**Thresholds (read-only here)**")
                    st.caption("Threshold sliders shown only in full mode to prevent duplicate widgets.")
                    st.metric("Strong Threshold", st.session_state.get('strong_threshold', 0.85))
                    st.metric("Moderate Threshold", st.session_state.get('moderate_threshold', 0.65))
                with col2:
                    st.write("**Env**")
                    st.session_state.custom_parameters["exclude_intraresidue"] = st.checkbox(
                        "Exclude Intraresidue", value=st.session_state.get('exclude_intraresidue', False), key="exclude_intraresidue")
                    st.session_state.custom_parameters["min_sequence_separation"] = st.slider(
                        "Min Seq Sep", min_value=0, max_value=5, value=st.session_state.get('min_sequence_separation', 1), step=1, key="min_sequence_separation")
            
            with col1:
                # Hydrogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('hydrogenbond', True):
                    st.write("**Hydrogen Bonds**")
                    st.session_state.custom_parameters["hbond_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.0, max_value=5.0, value=st.session_state.get('hbond_distance', 3.5), step=0.1,
                        key="hbond_distance", help="Maximum donor-acceptor distance"
                    )
                    st.session_state.custom_parameters["hbond_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=90, max_value=180, value=st.session_state.get('hbond_angle', 120), step=5,
                        key="hbond_angle", help="Minimum donor-H-acceptor angle"
                    )
                
                # Ionic Interactions
                st.write("**Ionic Interactions**")
                st.session_state.custom_parameters["ionic_distance"] = st.slider(
                    "Distance Cutoff (Å)",
                    min_value=3.0, max_value=8.0, value=st.session_state.get('ionic_distance', 6.0), step=0.1,
                    key="ionic_distance", help="Maximum charge-charge distance"
                )
            
            with col2:
                # Hydrophobic Contacts (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('hydrophobiccontact', True):
                    st.write("**Hydrophobic Contacts**")
                    st.session_state.custom_parameters["hydrophobic_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=6.0, value=st.session_state.get('hydrophobic_distance', 5.0), step=0.1,
                        key="hydrophobic_distance", help="Maximum carbon-carbon distance"
                    )
                
                # London Dispersion (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('dispersion', True):
                    st.write("**London Dispersion**")
                    st.session_state.custom_parameters["dispersion_min"] = st.slider(
                        "Min Distance (Å)",
                        min_value=3.0, max_value=4.0, value=st.session_state.get('dispersion_min', 3.5), step=0.1,
                        key="dispersion_min", help="Minimum contact distance"
                    )
                    st.session_state.custom_parameters["dispersion_max"] = st.slider(
                        "Max Distance (Å)",
                        min_value=4.0, max_value=6.0, value=st.session_state.get('dispersion_max', 5.0), step=0.1,
                        key="dispersion_max", help="Maximum contact distance"
                    )
        
        with parameter_tabs[1]:
            st.write("**σ-hole Mediated Interactions**")
            col1, col2 = st.columns(2)
            
            with col1:
                # Halogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('halogenbond', True):
                    st.write("**Halogen Bonds**")
                    st.session_state.custom_parameters["halogen_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=5.0, value=st.session_state.get('halogen_distance', 4.0), step=0.1,
                        key="halogen_distance", help="Maximum X···A distance"
                    )
                    st.session_state.custom_parameters["halogen_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=120, max_value=180, value=st.session_state.get('halogen_angle', 140), step=5,
                        key="halogen_angle", help="Minimum C-X···A angle"
                    )
                
                # Chalcogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('chalcogenbond', True):
                    st.write("**Chalcogen Bonds**")
                    st.session_state.custom_parameters["chalcogen_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=5.0, value=st.session_state.get('chalcogen_distance', 4.0), step=0.1,
                        key="chalcogen_distance", help="Maximum Ch···A distance"
                    )
                    st.session_state.custom_parameters["chalcogen_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=120, max_value=180, value=st.session_state.get('chalcogen_angle', 140), step=5,
                        key="chalcogen_angle", help="Minimum C-Ch···A angle"
                    )
            
            with col2:
                # Pnictogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('pnictogenbond', True):
                    st.write("**Pnictogen Bonds**")
                    st.session_state.custom_parameters["pnictogen_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=5.0, value=st.session_state.get('pnictogen_distance', 4.0), step=0.1,
                        key="pnictogen_distance", help="Maximum Pn···A distance"
                    )
                    st.session_state.custom_parameters["pnictogen_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=120, max_value=180, value=st.session_state.get('pnictogen_angle', 140), step=5,
                        key="pnictogen_angle", help="Minimum C-Pn···A angle"
                    )
                
                # Tetrel Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('tetrelbond', True):
                    st.write("**Tetrel Bonds**")
                    st.session_state.custom_parameters["tetrel_distance_min"] = st.slider(
                        "Min Distance (Å)",
                        min_value=2.0, max_value=3.0, value=st.session_state.get('tetrel_distance_min', 2.5), step=0.1,
                        key="tetrel_distance_min", help="Minimum T···A distance"
                    )
                    st.session_state.custom_parameters["tetrel_distance_max"] = st.slider(
                        "Max Distance (Å)",
                        min_value=3.0, max_value=4.5, value=st.session_state.get('tetrel_distance_max', 3.6), step=0.1,
                        key="tetrel_distance_max", help="Maximum T···A distance"
                    )
                    st.session_state.custom_parameters["tetrel_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=140, max_value=180, value=st.session_state.get('tetrel_angle', 160), step=5,
                        key="tetrel_angle", help="Minimum C-T···A angle"
                    )
        
        with parameter_tabs[2]:
            st.write("**π-System Related Interactions**")
            col1, col2 = st.columns(2)
            
            with col1:
                # π-π Stacking (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('pipi', True):
                    st.write("**π-π Stacking**")
                    st.session_state.custom_parameters["pi_pi_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=8.0, value=st.session_state.get('pi_pi_distance', 5.5), step=0.1,
                        key="pi_pi_distance", help="Maximum centroid-centroid distance"
                    )
                    st.session_state.custom_parameters["pi_pi_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=0, max_value=45, value=st.session_state.get('pi_pi_angle', 30), step=5,
                        key="pi_pi_angle", help="Maximum dihedral angle between rings"
                    )

                # C-H···π Interactions (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('chpi', True):
                    st.write("**C-H···π Interactions**")
                    st.session_state.custom_parameters["ch_pi_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=6.0, value=st.session_state.get('ch_pi_distance', 4.5), step=0.1,
                        key="ch_pi_distance", help="Maximum C···centroid distance"
                    )
                    st.session_state.custom_parameters["ch_pi_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=60, max_value=120, value=st.session_state.get('ch_pi_angle', 90), step=5,
                        key="ch_pi_angle", help="Maximum C-H···centroid angle"
                    )
            
            with col2:
                # Anion-π Interactions (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('anionpi', True):
                    st.write("**Anion-π Interactions**")
                    st.session_state.custom_parameters["anion_pi_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=7.0, value=st.session_state.get('anion_pi_distance', 5.0), step=0.1,
                        key="anion_pi_distance", help="Maximum anion-centroid distance"
                    )

                # n→π* Interactions (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('npistar', True):
                    st.write("**n→π* Interactions**")
                    st.session_state.custom_parameters["n_pi_star_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=4.5, value=st.session_state.get('n_pi_star_distance', 3.5), step=0.1,
                        key="n_pi_star_distance", help="Maximum n···C=O distance"
                    )
                    st.session_state.custom_parameters["n_pi_star_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=90, max_value=150, value=st.session_state.get('n_pi_star_angle', 120), step=5,
                        key="n_pi_star_angle", help="Minimum n···C=O angle"
                    )
        
        # Mirror final values into custom_parameters mapping (ensures any non-rendered defaults still captured)
        st.session_state.custom_parameters.update({k: st.session_state[k] for k in parameter_defaults.keys() if k in st.session_state})

    def _identify_hotspots(self, all_interactions: Dict[str, List[Any]], selected_interactions: List[str]) -> List[Dict[str, Any]]:
        """Identify interaction hotspots from a given set of interactions."""
        residue_counts = {}
        
        # Only consider selected interaction types that are also present in the results
        for itype, interactions in all_interactions.items():
            resolved_key = self._resolve_filter_key(itype)
            if resolved_key not in selected_interactions:
                continue
            
            for interaction in interactions:
                res1_id = f"{self._get_interaction_property(interaction, 'residue1', '')}{self._get_interaction_property(interaction, 'chain1', '')}"
                res2_id = f"{self._get_interaction_property(interaction, 'residue2', '')}{self._get_interaction_property(interaction, 'chain2', '')}"
                
                if res1_id:
                    residue_counts[res1_id] = residue_counts.get(res1_id, 0) + 1
                if res2_id:
                    residue_counts[res2_id] = residue_counts.get(res2_id, 0) + 1
        
        if not residue_counts:
            return []

        # Sort by interaction count
        sorted_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)
        
        hotspots = []
        for residue, count in sorted_residues:
            # A residue is a hotspot if it has a significant number of interactions
            if count >= 4:  # Hotspot threshold
                hotspots.append({
                    'residue': residue,
                    'interaction_count': count,
                })
                
        return hotspots

    def _get_interaction_inference(self, itype: str, count: int, strengths: Dict[str, int]) -> str:
        """Generate an insightful inference based on interaction data."""
        inferences = {
            'hydrogenbond': "A high number of hydrogen bonds, especially strong ones, is a primary contributor to the protein's structural stability and defines its secondary structure elements.",
            'ionicinteraction': "The presence of salt bridges on the protein surface can be crucial for protein-protein interactions, ligand binding, and maintaining stability across different pH environments.",
            'hydrophobiccontact': "Extensive hydrophobic contacts suggest a well-packed hydrophobic core, which is a major driving force for protein folding and overall stability.",
            'pipi': "π-π stacking interactions between aromatic residues are key for stabilizing protein structure, particularly in β-sheets, and are often found in ligand-binding pockets.",
            'halogenbond': "Halogen bonds are highly directional and can be critical for high-affinity ligand binding, acting as specific 'molecular velcro'.",
            'chalcogenbond': "Chalcogen bonds, involving sulfur or selenium, contribute to structural organization and can play a role in catalysis or redox sensing.",
            'pnictogenbond': "These interactions, involving elements like phosphorus, are less common but can provide unique geometric constraints important for specific molecular recognition events.",
            'tetrelbond': "Tetrel bonds involving carbon or silicon are subtle but can influence conformation and packing within the protein structure.",
            'anionpi': "Anion-π interactions are important for recognizing and binding negatively charged substrates or cofactors to aromatic rings.",
            'chpi': "C-H···π interactions are weak but numerous, collectively contributing significantly to the stability of protein structures and complexes.",
            'npistar': "n→π* interactions are short-range stereoelectronic interactions that can influence peptide bond conformation and protein secondary structure.",
            'dispersion': "London dispersion forces, while individually weak, are ubiquitous and collectively provide a significant stabilizing contribution to the overall protein structure."
        }
        
        resolved_key = self._resolve_filter_key(itype)
        base_inference = inferences.get(resolved_key, "")
        
        if not base_inference:
            return ""
            
        # Add more specific inferences based on counts and strengths
        if count > 50 and resolved_key == 'hydrogenbond':
            base_inference = f"With **{count}** hydrogen bonds detected, this interaction is a dominant force in maintaining the protein's architecture. " + base_inference
        
        strong_count = strengths.get('Strong', 0)
        if count > 0 and strong_count / count > 0.4 and resolved_key == 'hydrogenbond':
            base_inference += " The prevalence of strong bonds suggests a particularly rigid and stable structure."
        
        if count > 10 and resolved_key == 'ionicinteraction':
            base_inference = f"The presence of **{count}** salt bridges suggests they play a significant role in this protein's function, possibly in defining protein-protein interfaces or ensuring stability across a range of pH levels. " + base_inference

        return base_inference

    def _generate_analysis_summary(self, result: Dict[str, Any], selected_interactions: List[str]) -> str:
        """Generates a point-wise analysis summary for the given results."""
        pdb_id = result.get('pdb_id', 'Unknown')
        summary_lines = []

        # The result passed here should already be filtered by strength
        filtered_interactions = result.get('interactions', {})
        
        total_interactions = sum(len(interactions) for itype, interactions in filtered_interactions.items() if self._resolve_filter_key(itype) in selected_interactions)
        
        if total_interactions == 0:
            summary_lines.append("No interactions found with the current filter settings. Try adjusting the interaction types or strength filters in the sidebar or in the configuration expander below.")
            return "\n".join(summary_lines)
            
        summary_lines.append(f"**Overall, {total_interactions} interactions of the selected types were detected with the current filter settings.**\n")

        # Per-interaction analysis
        for itype in get_interaction_types(): # Iterate in a consistent order
            resolved_key = self._resolve_filter_key(itype)
            if resolved_key not in selected_interactions:
                continue

            # The key in filtered_interactions might have an underscore
            interactions = filtered_interactions.get(itype, [])
            if not interactions:
                # Try with the resolved key without underscore
                interactions = filtered_interactions.get(resolved_key, [])

            if not interactions:
                continue

            display_name = get_interaction_display_names().get(resolved_key, resolved_key.title())
            summary_lines.append(f"#### {display_name} ({len(interactions)} found)\n")

            # Stats
            distances = [self._get_interaction_property(i, 'distance', 0) for i in interactions if self._get_interaction_property(i, 'distance') is not None]
            if distances:
                avg_dist = np.mean(distances)
                min_dist = min(distances)
                max_dist = max(distances)
                summary_lines.append(f"- **Statistics**: Average distance of **{avg_dist:.2f} Å** (range: {min_dist:.2f} - {max_dist:.2f} Å).")

            strengths = [self._get_interaction_property(i, 'strength', 'unknown').capitalize() for i in interactions]
            if strengths:
                strength_counts = pd.Series(strengths).value_counts().to_dict()
                strength_str = ", ".join([f"**{count}** {name}" for name, count in strength_counts.items()])
                summary_lines.append(f"- **Strength Distribution**: {strength_str}.")

            # Inferences
            inference = self._get_interaction_inference(resolved_key, len(interactions), strength_counts)
            if inference:
                summary_lines.append(f"- **Inference**: {inference}\n")

        # Hotspot analysis
        hotspots = self._identify_hotspots(filtered_interactions, selected_interactions)
        if hotspots:
            summary_lines.append("### Key Findings & Hotspots\n")
            summary_lines.append("Certain residues act as critical hubs for interactions, potentially indicating functional or structural importance.\n")
            for i, hotspot in enumerate(hotspots[:3]): # Top 3
                res_name = hotspot['residue']
                count = hotspot['interaction_count']
                summary_lines.append(f"- **Hotspot {i+1}**: Residue **{res_name}** is highly active, participating in **{count} interactions**. This suggests it may be a key residue for stability or function.")
        
        return "\n".join(summary_lines)

    def _render_pointwise_analysis_summary(self):
        """Render the point-wise analysis summary for the current protein."""
        st.subheader(f"🔬 Point-wise Analysis for {st.session_state.current_pdb}")
        
        result = st.session_state.analysis_results[st.session_state.current_pdb]
        
        # Get currently selected interactions for the summary
        selected_interactions = [
            itype for itype, selected in st.session_state.get('selected_interactions', {}).items() if selected
        ]
        
        if selected_interactions:
            # Apply individual strength filters for the summary
            filtered_interactions = self._apply_individual_strength_filters(result.get('interactions', {}))
            
            # Create a temporary result object with filtered data for summary generation
            summary_result = result.copy()
            summary_result['interactions'] = filtered_interactions
            
            analysis_summary = self._generate_analysis_summary(summary_result, selected_interactions)
            st.markdown(analysis_summary, unsafe_allow_html=True)
        else:
            st.warning("No interaction types selected. Please select interactions in the sidebar to see an analysis.")


    def _render_visualization_tab(self):
        """Render the visualization tab with advanced controls and accessibility."""
        if not st.session_state.current_pdb:
            st.info("Select a structure from the Analysis tab to visualize")
            return

        st.header(f"📊 Visualization - {st.session_state.current_pdb}")
        # Open a scoped container to limit CSS/JS fixes to the entire Visualization tab
        st.markdown('<div id="viz-scope">', unsafe_allow_html=True)

        # Accessibility: Keep headings semantic; avoid extra empty containers that can create layout artifacts

        # Show active individual filters summary
        active_filters = {itype: filter_val for itype, filter_val in st.session_state.individual_strength_filters.items() if filter_val != 'all'}
        if active_filters:
            st.info("🔍 **Individual strength filters active:**")
            filter_groups = {}
            for itype, filter_val in active_filters.items():
                if filter_val not in filter_groups:
                    filter_groups[filter_val] = []
                filter_groups[filter_val].append(get_interaction_display_names()[itype])
            strength_names = {
                'strong': 'Strong only',
                'strong_moderate': 'Strong & Moderate',
                'all': 'All (Strong, Moderate, Weak)'
            }
            for filter_val, types in filter_groups.items():
                st.caption(f"**{strength_names.get(filter_val)}**: {', '.join(types)}")

        # 3D Structure viewer with extra controls
        result = st.session_state.analysis_results[st.session_state.current_pdb]
        if 'interactions' in result and result['interactions']:
            filtered_interactions = self._apply_individual_strength_filters(result['interactions'])
            filtered_result = result.copy()
            filtered_result['interactions'] = filtered_interactions
        else:
            filtered_result = result

        # Render 3D viewer with filtered interactions and extra info/links
        self.structure_viewer.render_structure(
            st.session_state.current_pdb,
            filtered_result,
            [itype for itype, selected in st.session_state.selected_interactions.items() if selected]
        )

        # Enhanced CSS specifically for dropdown text visibility (scoped to #viz-scope)
        st.markdown(
            """
            <style>
            /* General labels inside Visualization tab: bold, white, slightly larger */
            #viz-scope label, #viz-scope .stCheckbox label {
                color: #F5F6FA !important;
                font-weight: 700 !important;
                font-size: 0.98rem !important;
            }

            /* Comprehensive dropdown text visibility fix */
            #viz-scope [data-testid="stSelectbox"] label {
                color: #F5F6FA !important;
                font-weight: 700 !important;
                font-size: 0.98rem !important;
            }
            
            /* Target all selectbox input containers */
            #viz-scope [data-testid="stSelectbox"] > div > div[data-baseweb="select"] {
                background: rgba(36,37,38,0.85) !important;
                color: #F5F6FA !important;
                border: 1.5px solid rgba(0,180,216,0.3) !important;
                border-radius: 8px !important;
            }
            
            /* Target the text content specifically */
            #viz-scope [data-testid="stSelectbox"] > div > div[data-baseweb="select"] > div {
                color: #F5F6FA !important;
                background: transparent !important;
                font-size: 0.98rem !important;
                font-weight: 600 !important;
            }

            /* Ensure input/value text is visible (covers WebKit fill) */
            #viz-scope [data-testid="stSelectbox"] div[data-baseweb="select"] input {
                color: #F5F6FA !important;
                -webkit-text-fill-color: #F5F6FA !important;
                caret-color: #F5F6FA !important;
                text-shadow: none !important;
                opacity: 1 !important;
                font-size: 0.98rem !important;
                font-weight: 600 !important;
            }

            /* Selected value and placeholder across dynamic classnames */
            #viz-scope [data-testid="stSelectbox"] [class*="singleValue"],
            #viz-scope [data-testid="stSelectbox"] [class*="SingleValue"],
            #viz-scope [data-testid="stSelectbox"] [class*="value"],
            #viz-scope [data-testid="stSelectbox"] [class*="placeholder"],
            #viz-scope [data-testid="stSelectbox"] [class*="Placeholder"] {
                color: #F5F6FA !important;
                -webkit-text-fill-color: #F5F6FA !important;
                opacity: 1 !important;
                text-shadow: none !important;
                font-size: 0.98rem !important;
                font-weight: 600 !important;
                /* prevent ellipsis-only visibility */
                overflow: visible !important;
                text-overflow: clip !important;
                white-space: nowrap !important;
                max-width: 100% !important;
            }
            
            /* Target nested spans and divs that contain the actual text */
            #viz-scope [data-testid="stSelectbox"] > div > div[data-baseweb="select"] span,
            #viz-scope [data-testid="stSelectbox"] > div > div[data-baseweb="select"] div,
            #viz-scope [data-testid="stSelectbox"] > div > div[data-baseweb="select"] > div > div {
                color: #F5F6FA !important;
                background: transparent !important;
            }
            
            /* Dropdown menu styling (when popover is rendered inside scope) */
            #viz-scope [data-testid="stSelectbox"] div[role="listbox"] {
                background: rgba(36,37,38,0.95) !important;
                border: 1px solid #00B4D8 !important;
                border-radius: 8px !important;
            }
            
            /* Individual dropdown options */
            #viz-scope [data-testid="stSelectbox"] div[role="listbox"] > div,
            #viz-scope [data-testid="stSelectbox"] div[role="option"] {
                color: #F5F6FA !important;
                background: rgba(36,37,38,0.85) !important;
                padding: 8px 12px !important;
            }

            /* Option text, including nested spans */
            #viz-scope [data-testid="stSelectbox"] [role="option"],
            #viz-scope [data-testid="stSelectbox"] [role="option"] * {
                color: #F5F6FA !important;
                -webkit-text-fill-color: #F5F6FA !important;
                text-shadow: none !important;
            }
            
            /* Hover effect for options */
            #viz-scope [data-testid="stSelectbox"] div[role="option"]:hover {
                background: rgba(0,180,216,0.3) !important;
                color: #FFFFFF !important;
            }
            
            /* Arrow/chevron icon */
            #viz-scope [data-testid="stSelectbox"] svg {
                fill: #F5F6FA !important;
                color: #F5F6FA !important;
            }
            
            /* Focus state */
            #viz-scope [data-testid="stSelectbox"] > div > div[data-baseweb="select"]:focus-within {
                border-color: #00B4D8 !important;
                box-shadow: 0 0 0 2px rgba(0,180,216,0.3) !important;
            }
            
            /* Ensure all text elements are visible */
            #viz-scope [data-testid="stSelectbox"] * {
                color: #F5F6FA !important;
            }
            
            /* Override any conflicting styles */
            #viz-scope [data-testid="stSelectbox"] [class*="singleValue"],
            #viz-scope [data-testid="stSelectbox"] [class*="placeholder"] {
                color: #F5F6FA !important;
            }

            /* Global fallback: BaseWeb portal renders dropdown outside #viz-scope */
            ul[role="listbox"] {
                background: rgba(36,37,38,0.95) !important;
                border: 1px solid rgba(0,180,216,0.35) !important;
                border-radius: 8px !important;
            }
            ul[role="listbox"] li,
            div[role="option"] {
                color: #F5F6FA !important;
                -webkit-text-fill-color: #F5F6FA !important;
                background: rgba(36,37,38,0.85) !important;
            }
            ul[role="listbox"] li:hover,
            div[role="option"]:hover {
                background: rgba(0,180,216,0.3) !important;
                color: #FFFFFF !important;
            }
            </style>
            """,
            unsafe_allow_html=True
        )


        # External database links
        pdb_id = st.session_state.current_pdb
        st.markdown(f"[🔗 View {pdb_id} on RCSB PDB](https://www.rcsb.org/structure/{pdb_id}) | [UniProt Search](https://www.uniprot.org/uniprot/?query={pdb_id}) | [Literature](https://pubmed.ncbi.nlm.nih.gov/?term={pdb_id})")

        # Interaction plots
        st.subheader("📈 Interaction Analysis Plots")
        tab1, tab2, tab3, tab4 = st.tabs(["🔥 Chain Heatmap", "📊 Distribution", "🔬 Ramachandran", "🕸️ Network"])
        with tab1:
            st.subheader("Chain Interaction Heatmap")
            st.caption("Shows interactions between different protein chains")
            self.plots.render_chain_heatmap(filtered_result)
        with tab2:
            st.subheader("Interaction Distribution")
            st.caption("Distribution of different interaction types")
            self.plots.render_interaction_distribution(filtered_result)
        with tab3:
            st.subheader("Ramachandran Plot")
            st.caption("Protein backbone conformation analysis")
            self.plots.render_ramachandran_plot(result)
        with tab4:
            st.subheader("Interaction Network")
            self.plots.render_interaction_network(filtered_result)

        # Close scoped container
        st.markdown('</div>', unsafe_allow_html=True)

    

    def _render_results_tab(self):
        """Render the results summary tab with annotation, legend, progress, and diff features."""
        if not st.session_state.analysis_results:
            st.info("No analysis results available")
            return

        st.header("📋 Analysis Results")
        # Unified interaction legend (always visible for orientation)
        render_interaction_legend()

        # Show active individual filters summary for context
        active_filters = {itype: filter_val for itype, filter_val in st.session_state.individual_strength_filters.items() if filter_val != 'all'}
        if active_filters:
            st.info("🔍 **Note:** Individual strength filters are active. The table below reflects these settings.")
            filter_groups = {}
            for itype, filter_val in active_filters.items():
                if filter_val not in filter_groups:
                    filter_groups[filter_val] = []
                filter_groups[filter_val].append(get_interaction_display_names()[itype])
            for filter_val, types in filter_groups.items():
                strength_names = {
                    'strong': 'Strong only',
                    'strong_moderate': 'Strong & Moderate',
                    'all': 'All (Strong, Moderate, Weak)'
                }
                st.caption(f"**{strength_names.get(filter_val)}**: {', '.join(types)}")
        else:
            st.info("🔍 Showing all interaction strengths for all types. Use the sidebar to configure filters.")

        # Results for current structure
        if st.session_state.current_pdb:
            result = st.session_state.analysis_results[st.session_state.current_pdb]
            # Compute extensions lazily (non-destructive to core analysis) based on feature toggles
            self._ensure_extensions(st.session_state.current_pdb, result)
            st.subheader(f"Interactions in {st.session_state.current_pdb}")
            interactions_df = self._create_interactions_dataframe(result)
            
            # Add interaction type filter (inter-chain vs intra-chain)
            interaction_category = st.radio(
                "Show:",
                ["All Interactions", "Inter-chain Only", "Intra-chain Only"],
                horizontal=True,
                help="Filter by whether interactions occur between different chains or within the same chain"
            )
            
            # Interaction table filters (chain, strength, residue substring)
            if not interactions_df.empty:
                cols_f = st.columns(3)
                with cols_f[0]:
                    chain_filter = st.selectbox("Chain Filter", options=["All"] + sorted({c for c in interactions_df.get('chain', []) if isinstance(c, str)}), index=0)
                with cols_f[1]:
                    strength_filter = st.selectbox("Strength", options=["All","Strong","Moderate","Weak"], index=0)
                with cols_f[2]:
                    residue_search = st.text_input("Residue Contains", placeholder="e.g. ASP")
                
                # Apply interaction category filter
                if 'Chain 1' in interactions_df.columns and 'Chain 2' in interactions_df.columns:
                    if interaction_category == "Inter-chain Only":
                        interactions_df = interactions_df[interactions_df['Chain 1'] != interactions_df['Chain 2']]
                    elif interaction_category == "Intra-chain Only":
                        interactions_df = interactions_df[interactions_df['Chain 1'] == interactions_df['Chain 2']]
                
                if chain_filter != "All" and 'chain' in interactions_df.columns:
                    interactions_df = interactions_df[interactions_df['chain'] == chain_filter]
                if strength_filter != "All" and 'strength' in interactions_df.columns:
                    interactions_df = interactions_df[interactions_df['strength'].str.capitalize() == strength_filter]
                if residue_search and 'residue' in interactions_df.columns:
                    interactions_df = interactions_df[interactions_df['residue'].str.contains(residue_search, case=False, na=False)]
            if not interactions_df.empty:
                total_interactions = sum(len(interactions) for interactions in result.get('interactions', {}).values())
                filtered_count = len(interactions_df)
                
                # Calculate inter-chain and intra-chain counts from the original unfiltered dataframe
                original_df = self._create_interactions_dataframe(result)
                if 'Chain 1' in original_df.columns and 'Chain 2' in original_df.columns:
                    inter_chain_count = len(original_df[original_df['Chain 1'] != original_df['Chain 2']])
                    intra_chain_count = len(original_df[original_df['Chain 1'] == original_df['Chain 2']])
                    
                    # Display metrics
                    col_m1, col_m2, col_m3 = st.columns(3)
                    with col_m1:
                        st.metric("Total Interactions", total_interactions)
                    with col_m2:
                        st.metric("Inter-chain", inter_chain_count, help="Interactions between different chains")
                    with col_m3:
                        st.metric("Intra-chain", intra_chain_count, help="Interactions within the same chain")
                
                if active_filters:
                    st.caption(f"Showing {filtered_count} of {total_interactions} total interactions after applying filters")
                
                # Create display dataframe without the Chain 1 and Chain 2 columns
                display_df = interactions_df.drop(columns=['Chain 1', 'Chain 2'], errors='ignore')
                
                # Virtualized table for large datasets
                if len(display_df) > 1200:
                    try:
                        try:
                            from st_aggrid import AgGrid, GridOptionsBuilder  # type: ignore
                            gob = GridOptionsBuilder.from_dataframe(display_df)
                            gob.configure_default_column(resizable=True, sortable=True, filter=True)
                            AgGrid(display_df, gridOptions=gob.build(), height=420, theme='streamlit')
                        except ModuleNotFoundError:
                            st.info("Install streamlit-aggrid to enable virtualized large table rendering")
                            st.dataframe(display_df, width="stretch")
                    except Exception:
                        st.dataframe(display_df, width="stretch")
                else:
                    st.dataframe(display_df, width="stretch")
                # Export, copy, and share options
                col1, col2, col3, col4, col5 = st.columns(5)
                with col1:
                    csv_data = interactions_df.to_csv(index=False)
                    st.download_button(
                        "📄 Download CSV",
                        csv_data,
                        f"{st.session_state.current_pdb}_interactions.csv",
                        "text/csv"
                    )
                with col2:
                    excel_data = io.BytesIO()
                    interactions_df.to_excel(excel_data, index=False)
                    st.download_button(
                        "📊 Download Excel",
                        excel_data.getvalue(),
                        f"{st.session_state.current_pdb}_interactions.xlsx",
                        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                with col3:
                    st.code(csv_data, language="csv")
                    st.caption("Copy CSV to clipboard for sharing.")
                with col4:
                    st.markdown(f"[🔗 Share this PDB on RCSB](https://www.rcsb.org/structure/{st.session_state.current_pdb})")
                with col5:
                    # Bundle export (zip JSON + CSV)
                    import zipfile, json, tempfile, base64, numpy as _np, datetime as _dt
                    def _safe(obj):
                        # Handle numpy scalar types and arrays
                        if isinstance(obj, (_np.integer,)):
                            return int(obj)
                        if isinstance(obj, (_np.floating,)):
                            return float(obj)
                        if isinstance(obj, (_np.ndarray,)):
                            return obj.tolist()
                        # Dataclass-like objects with __dict__
                        if hasattr(obj, '__dict__'):
                            return {k: v for k,v in obj.__dict__.items() if not k.startswith('_')}
                        # Fallback to string
                        return str(obj)
                    try:
                        data_json = json.dumps(result, indent=2, default=_safe)
                    except Exception as _serr:
                        data_json = json.dumps({'serialization_error': str(_serr), 'timestamp': _dt.datetime.utcnow().isoformat(), 'fallback_str': str(result)[:10000]}, indent=2)
                    with tempfile.TemporaryDirectory() as td:
                        import pathlib
                        p_json = pathlib.Path(td)/f"{st.session_state.current_pdb}_result.json"
                        p_csv = pathlib.Path(td)/f"{st.session_state.current_pdb}_interactions.csv"
                        p_json.write_text(data_json)
                        p_csv.write_text(csv_data)
                        zip_path = pathlib.Path(td)/f"{st.session_state.current_pdb}_bundle.zip"
                        with zipfile.ZipFile(zip_path, 'w') as zf:
                            zf.write(p_json, p_json.name)
                            zf.write(p_csv, p_csv.name)
                        st.download_button("📦 Bundle", data=zip_path.read_bytes(), file_name=zip_path.name, mime="application/zip", help="Download JSON + CSV bundle")
            else:
                st.warning("No interactions found with current filters")

        # Bookmarks and notes
        self._render_bookmarks_section()

        # Detector Progress Panel (if streaming events captured)
        prog_events = st.session_state.get('detector_progress_panel')
        # Periodic auto-persist every 5 minutes if events present
        try:
            import time as _t, json as _json, pathlib as _pl
            if prog_events:
                now = _t.time()
                last = st.session_state.get('_last_progress_autosave', 0)
                interval = 300  # 5 minutes
                if now - last > interval:
                    outp = _pl.Path('cache')/ 'progress_events.json'
                    outp.parent.mkdir(parents=True, exist_ok=True)
                    outp.write_text(_json.dumps(prog_events[-5000:], indent=2))  # cap persisted size
                    st.session_state._last_progress_autosave = now
        except Exception:
            pass
        if prog_events:
            with st.expander("⚙️ Detector Progress", expanded=False):
                try:
                    df_prog = pd.DataFrame(prog_events)
                    if 'timestamp' in df_prog.columns:
                        df_prog['ts_readable'] = df_prog['timestamp'].apply(lambda x: _dt.datetime.fromtimestamp(x).strftime('%H:%M:%S'))
                    cols_pref = [c for c in ['ts_readable','detector','count','duration','pdb_id','dry_run','run_id'] if c in df_prog.columns]
                    st.dataframe(df_prog[cols_pref], use_container_width=True)
                    # Aggregations
                    agg = df_prog.groupby('detector').agg(total_events=('detector','count'), avg_duration=('duration','mean'), last_count=('count','last'))
                    st.markdown("**Per-Detector Summary**")
                    st.dataframe(agg, use_container_width=True)
                    if st.checkbox("Group by Run", value=False):
                        if 'run_id' in df_prog.columns:
                            run_agg = df_prog.groupby('run_id').agg(events=('detector','count'), detectors=('detector','nunique'), total_duration=('duration','sum'))
                            st.dataframe(run_agg, use_container_width=True)
                    colp1, colp2, colp3 = st.columns(3)
                    with colp1:
                        if st.button("Reset Progress Log", help="Clear accumulated detector progress events"):
                            st.session_state.detector_progress_panel = []
                            st.experimental_rerun()
                    with colp2:
                        if st.button("Persist Progress Log", help="Write events to disk under cache/progress_events.json"):
                            try:
                                import json, pathlib
                                outp = pathlib.Path('cache')/ 'progress_events.json'
                                outp.parent.mkdir(parents=True, exist_ok=True)
                                outp.write_text(json.dumps(df_prog.to_dict(orient='records'), indent=2))
                                st.success(f"Saved {len(df_prog)} events to {outp}")
                            except Exception as _e:
                                st.error(f"Persist failed: {_e}")
                    with colp3:
                        if st.button("Rotate & Cap", help="Trim oldest events keeping latest N (default 2000)"):
                            cap = st.number_input("Max Events", min_value=100, max_value=20000, value=2000, step=100, key="progress_cap")
                            if len(st.session_state.detector_progress_panel) > cap:
                                st.session_state.detector_progress_panel = st.session_state.detector_progress_panel[-cap:]
                                st.success(f"Trimmed log to {cap} most recent events")
                                st.experimental_rerun()
                    with colp3:
                        if st.button("Load Persisted Log", help="Append events from disk if available", key="load_persisted_log"):
                            try:
                                import json, pathlib
                                inp = pathlib.Path('cache')/ 'progress_events.json'
                                if inp.exists():
                                    data = json.loads(inp.read_text())
                                    st.session_state.detector_progress_panel.extend(data)
                                    st.success(f"Loaded {len(data)} events from disk")
                                    st.experimental_rerun()
                                else:
                                    st.info("No persisted file present")
                            except Exception as _e:
                                st.error(f"Load failed: {_e}")
                except Exception as e:
                    st.warning(f"Unable to render progress panel: {e}")

        # Comparative Diff View (structure-level gain/loss)
        if len(st.session_state.analysis_results) >= 2:
            with st.expander("🔄 Comparative Diff (Two Structures)", expanded=False):
                diff_filter = st.text_input("Filter signatures", key="diff_filter_sig", placeholder="Substring e.g. ALA|LYS")
                structure_ids = list(st.session_state.analysis_results.keys())
                diff_c1, diff_c2 = st.columns(2)
                with diff_c1:
                    struct_a = st.selectbox("Structure A", structure_ids, key="diff_struct_a")
                with diff_c2:
                    struct_b = st.selectbox("Structure B", structure_ids, key="diff_struct_b")
                if struct_a != struct_b:
                    r_a = st.session_state.analysis_results.get(struct_a, {})
                    r_b = st.session_state.analysis_results.get(struct_b, {})
                    ints_a = r_a.get('interactions', {})
                    ints_b = r_b.get('interactions', {})
                    added_rows, lost_rows = [], []
                    def _sig(rec):
                        """Produce a coordinate-aware stable signature for interaction record.
                        Falls back to residue/chain composite if no coordinates present.
                        Expected dict fields (if available): residue1, residue2, chain1, chain2, atom1, atom2, x1,y1,z1,x2,y2,z2.
                        """
                        import hashlib
                        if isinstance(rec, dict):
                            # Prefer coordinate-based hash if coordinate pairs present
                            coord_fields = ['x1','y1','z1','x2','y2','z2']
                            if all(f in rec for f in coord_fields):
                                coord_part = ','.join(f"{rec[f]:.2f}" if isinstance(rec[f], (int,float)) else str(rec[f]) for f in coord_fields)
                                base_meta = '|'.join(str(rec.get(k,'')) for k in ['residue1','residue2','chain1','chain2','atom1','atom2'])
                                h = hashlib.sha1(f"{base_meta}|{coord_part}".encode()).hexdigest()[:10]
                                return f"{base_meta}|{coord_part}|{h}"
                            # Fallback textual composite
                            parts = []
                            for key in ['residue1','residue2','residue','resname','chain1','chain2','chain','atom1','atom2']:
                                if key in rec and rec[key] is not None:
                                    parts.append(str(rec[key]))
                            if parts:
                                return '|'.join(parts)
                            return str(rec)
                        return str(rec)
                    all_keys = set(ints_a.keys()) | set(ints_b.keys())
                    for k in all_keys:
                        set_a = {_sig(x) for x in ints_a.get(k, [])}
                        set_b = {_sig(x) for x in ints_b.get(k, [])}
                        for sig in (set_b - set_a):
                            added_rows.append({'interaction': k, 'signature': sig})
                        for sig in (set_a - set_b):
                            lost_rows.append({'interaction': k, 'signature': sig})
                    if added_rows or lost_rows:
                        col_add, col_lost = st.columns(2)
                        import pandas as _pd, io as _io
                        df_added = _pd.DataFrame(added_rows)
                        df_lost = _pd.DataFrame(lost_rows)
                        if diff_filter:
                            if not df_added.empty:
                                df_added = df_added[df_added['signature'].str.contains(diff_filter, case=False, na=False)]
                            if not df_lost.empty:
                                df_lost = df_lost[df_lost['signature'].str.contains(diff_filter, case=False, na=False)]
                        # Mini summary counts
                        from collections import Counter
                        add_cnt = Counter([r['interaction'] for r in added_rows])
                        lost_cnt = Counter([r['interaction'] for r in lost_rows])
                        summary_df = _pd.DataFrame([{ 'interaction':k, 'added':add_cnt.get(k,0), 'lost':lost_cnt.get(k,0)} for k in set(add_cnt)|set(lost_cnt)]).sort_values('added', ascending=False)
                        st.markdown("**Interaction Type Change Summary**")
                        st.dataframe(summary_df, use_container_width=True)
                        with col_add:
                            st.markdown(f"**Added (B vs A)** – {len(added_rows)}")
                            if not df_added.empty:
                                st.dataframe(df_added[:300], use_container_width=True)
                            else:
                                st.caption("No new interactions")
                        with col_lost:
                            st.markdown(f"**Lost (A vs B)** – {len(lost_rows)}")
                            if not df_lost.empty:
                                st.dataframe(df_lost[:300], use_container_width=True)
                            else:
                                st.caption("No lost interactions")
                        # Export combined diff
                        combined = []
                        for row in added_rows:
                            r = dict(row)
                            r['change'] = 'added'
                            combined.append(r)
                        for row in lost_rows:
                            r = dict(row)
                            r['change'] = 'lost'
                            combined.append(r)
                        if combined:
                            df_combined = _pd.DataFrame(combined)
                            csv_bytes = df_combined.to_csv(index=False).encode()
                            st.download_button("📄 Download Diff CSV", data=csv_bytes, file_name=f"diff_{struct_a}_vs_{struct_b}.csv", mime="text/csv")
                    else:
                        st.success("Structures share identical interaction signatures (based on simplified diff heuristic).")
                else:
                    st.info("Select two different structures to compute diff")

        # Advanced extension analytics section
        if st.session_state.current_pdb:
            ext_result = st.session_state.analysis_results[st.session_state.current_pdb].get('extensions', {})
            if ext_result:
                st.write("---")
                st.header("🧪 Advanced Protein Interaction Extensions")
                # Residue Profiles
                if 'residue_profiles' in ext_result:
                    with st.expander("Residue Interaction Profiles", expanded=False):
                        profiles = ext_result['residue_profiles'].get('profiles', {})
                        if profiles:
                            import pandas as pd
                            rows = []
                            for residue, rec in profiles.items():
                                base = {
                                    'Residue': residue,
                                    'Total': rec.get('total', 0)
                                }
                                # Flatten top interaction types (limit to 4 for width)
                                by_type = rec.get('by_type', {})
                                top_items = sorted(by_type.items(), key=lambda x: x[1], reverse=True)[:4]
                                for i, (itype, cnt) in enumerate(top_items, 1):
                                    base[f'Type{i}'] = f"{itype}:{cnt}"
                                rows.append(base)
                            df_profiles = pd.DataFrame(rows).sort_values('Total', ascending=False)
                            st.dataframe(df_profiles, width='stretch')
                            csv_bytes = df_profiles.to_csv(index=False).encode()
                            st.download_button("Download Profiles CSV", csv_bytes, f"{st.session_state.current_pdb}_residue_profiles.csv", "text/csv")
                        else:
                            st.info("No residue profile data available")
                # Interface Analysis
                if 'interface_analysis' in ext_result:
                    with st.expander("Chain-Chain Interface Analysis", expanded=False):
                        interfaces = ext_result['interface_analysis'].get('interfaces', [])
                        if interfaces:
                            import pandas as pd
                            df_if = pd.DataFrame(interfaces)
                            st.dataframe(df_if, width='stretch')
                            st.download_button("Download Interfaces CSV", df_if.to_csv(index=False).encode(), f"{st.session_state.current_pdb}_interfaces.csv", "text/csv")
                        else:
                            st.info("No inter-chain interfaces detected")
                # Outliers
                if 'outliers' in ext_result:
                    with st.expander("Threshold-Adjacent / Outlier Interactions", expanded=False):
                        outliers = ext_result['outliers'].get('records', [])
                        if outliers:
                            import pandas as pd
                            df_out = pd.DataFrame(outliers)
                            st.dataframe(df_out, width='stretch')
                            st.download_button("Download Outliers CSV", df_out.to_csv(index=False).encode(), f"{st.session_state.current_pdb}_outliers.csv", "text/csv")
                        else:
                            st.info("No near-threshold interactions flagged")
                # Provenance
                if 'provenance' in ext_result:
                    with st.expander("Analysis Provenance", expanded=False):
                        prov = ext_result['provenance']
                        st.json(prov)
                # Hydrogen bond subtypes
                if 'hbond_subtypes' in ext_result:
                    with st.expander("Hydrogen Bond Subtypes", expanded=False):
                        hbext = ext_result['hbond_subtypes']
                        counts = hbext.get('counts', {})
                        total_hb = hbext.get('total_hbonds', 0)
                        if total_hb > 0:
                            import pandas as pd
                            rows = []
                            fr = hbext.get('fractions', {})
                            for subtype, count in counts.items():
                                rows.append({
                                    'Subtype': subtype,
                                    'Count': count,
                                    'Fraction': f"{fr.get(subtype,0):.2f}"
                                })
                            df_sub = pd.DataFrame(rows).sort_values('Count', ascending=False)
                            st.metric("Total H-Bonds Classified", total_hb)
                            st.dataframe(df_sub, width='stretch')
                            st.download_button(
                                "Download H-Bond Subtypes CSV",
                                df_sub.to_csv(index=False).encode(),
                                f"{st.session_state.current_pdb}_hbond_subtypes.csv",
                                "text/csv"
                            )
                            if st.checkbox("Show Annotated H-Bonds", value=False, key="show_hbond_annotated"):
                                annotated = hbext.get('annotated', [])
                                if annotated:
                                    df_ann = pd.DataFrame(annotated)
                                    st.dataframe(df_ann, use_container_width=True)
                        else:
                            st.info("No hydrogen bonds available to classify")
                # Motifs (experimental)
                if 'motifs' in ext_result:
                    with st.expander("Structural Motifs (Experimental)", expanded=False):
                        motifs = ext_result['motifs'].get('motifs', [])
                        if motifs:
                            import pandas as pd
                            df_m = pd.DataFrame(motifs)
                            st.dataframe(df_m, width='stretch')
                            st.download_button("Download Motifs CSV", df_m.to_csv(index=False).encode(), f"{st.session_state.current_pdb}_motifs.csv", "text/csv")
                        else:
                            st.info("No motifs detected with current heuristic")

    def _ensure_extensions(self, pdb_id: str, result: Dict[str, Any]):
        """Compute extension analytics once per structure based on feature toggles.

        Adds results under result['extensions'] without modifying core 'interactions'.
        Safe to call repeatedly (idempotent)."""
        try:
            ext_store = result.setdefault('extensions', {})
            # --------------------------------------------------------------
            # Lazy backfill: if structural coordinate data absent, attempt
            # to reconstruct a minimal internal representation from the
            # currently cached / loaded Bio.PDB structure using PDBHandler.
            # This enables secondary structure, geometry, SASA, pockets, etc.
            # --------------------------------------------------------------
            if (not result.get('structures')) and hasattr(self, 'pdb_handler'):
                try:
                    biostruct = self.pdb_handler.fetch_structure_variant(pdb_id,
                                                                         assembly=self.config.default_assembly,
                                                                         include_ligands=self.config.include_ligands,
                                                                         exclude_waters=self.config.exclude_waters)
                    if biostruct is not None:
                        internal = self.pdb_handler.to_internal_representation(biostruct)
                        if internal:
                            result['structures'] = internal
                except Exception as bf_exc:
                    logger.warning(f"Backfill of structure coords failed for {pdb_id}: {bf_exc}")
            # Residue Profiles
            if self.config.enable_residue_profiles and 'residue_profiles' not in ext_store:
                ext_store['residue_profiles'] = compute_residue_profiles(result, self.config)
            # Interface Analysis
            if self.config.enable_interface_analysis and 'interface_analysis' not in ext_store:
                ext_store['interface_analysis'] = compute_interface_analysis(result, self.config)
            # Outlier Detection
            if self.config.enable_outlier_detection and 'outliers' not in ext_store:
                ext_store['outliers'] = compute_outliers(result, self.config)
            # Provenance Panel
            if self.config.enable_provenance_panel and 'provenance' not in ext_store:
                ext_store['provenance'] = compute_provenance(result, self.config)
            # Motif Detection (experimental)
            if self.config.enable_motif_detection and 'motifs' not in ext_store:
                try:
                    ext_store['motifs'] = compute_motifs(result, self.config)
                except Exception as me:
                    ext_store['motifs'] = {'error': str(me), 'motif_count': 0, 'motifs': []}
            # Structural quality / annotation modules
            if self.config.enable_secondary_structure and 'secondary_structure' not in ext_store:
                ext_store['secondary_structure'] = compute_secondary_structure(result, self.config)
            # Repair / augment fractions if pre-patch cached result lacks suffixed keys
            elif 'secondary_structure' in ext_store:
                try:
                    ss_ext = ext_store['secondary_structure']
                    counts = ss_ext.get('counts', {}) or {}
                    fr = ss_ext.get('fractions', {}) or {}
                    total_ct = sum(counts.values())
                    needs_patch = False
                    # Missing any *_frac key OR all zero while counts non-zero
                    if any(f"{k}_frac" not in fr for k in ['H','E','C']):
                        needs_patch = True
                    else:
                        if total_ct > 0 and all(fr.get(f"{k}_frac",0) == 0 for k in ['H','E','C']) and any(counts.get(k,0) > 0 for k in ['H','E','C']):
                            needs_patch = True
                    if needs_patch and total_ct > 0:
                        for k in ['H','E','C']:
                            fr[f"{k}_frac"] = counts.get(k,0)/total_ct
                        ss_ext['fractions'] = fr
                except Exception as rep_exc:
                    logger.warning(f"Secondary structure fraction repair failed for {pdb_id}: {rep_exc}")
            if self.config.enable_sasa_bsa and 'sasa_bsa' not in ext_store:
                ext_store['sasa_bsa'] = compute_sasa_bsa(result, self.config)
            if self.config.enable_geometry_quality and 'geometry_quality' not in ext_store:
                # ensure secondary structure torsions present first
                if 'secondary_structure' not in ext_store and self.config.enable_secondary_structure:
                    ext_store['secondary_structure'] = compute_secondary_structure(result, self.config)
                ext_store['geometry_quality'] = compute_geometry_quality(result, self.config)
            if self.config.enable_disulfide_analysis and 'disulfides' not in ext_store:
                ext_store['disulfides'] = compute_disulfides(result, self.config)
            if self.config.enable_pocket_detection and 'pockets' not in ext_store:
                ext_store['pockets'] = compute_pockets(result, self.config)
            if self.config.enable_conservation and 'conservation' not in ext_store:
                try:
                    ext_store['conservation'] = compute_conservation(result, self.config)
                except Exception as ce:
                    ext_store['conservation'] = {'error': str(ce), 'residues': {}, 'note': 'conservation stub error'}
            if self.config.enable_pi_pi_refinement and 'pi_pi_refinement' not in ext_store:
                try:
                    ext_store['pi_pi_refinement'] = compute_pi_pi_refinement(result, self.config)
                except Exception as pe:
                    ext_store['pi_pi_refinement'] = {'error': str(pe), 'note': 'pi-pi refinement failed'}
            # Hydrogen bond subtype classification
            if self.config.enable_hbond_subtypes and 'hbond_subtypes' not in ext_store:
                try:
                    ext_store['hbond_subtypes'] = compute_hbond_subtypes(result, self.config)
                except Exception as he:
                    ext_store['hbond_subtypes'] = {'error': str(he), 'counts': {}, 'total_hbonds': 0, 'fractions': {}, 'annotated': []}
        except Exception as e:
            logger.error(f"Extension computation failed for {pdb_id}: {e}")
    
    def _render_reports_tab(self):
        """Render the reports generation tab."""
        if not st.session_state.analysis_results:
            st.info("No analysis results available for reporting")
            return
        
        st.header("📄 Report Generation")
        
        # Report options
        report_type = st.selectbox(
            "Report Type:",
            ["PDF Report", "PowerPoint Presentation", "LaTeX Export", "Complete Package"]
        )
        
        # Structure selection for reporting
        available_structures = list(st.session_state.analysis_results.keys())
        
        if len(available_structures) == 1:
            selected_structures = available_structures
            st.info(f"Generating report for: {available_structures[0]}")
        else:
            selected_structures = st.multiselect(
                "Select structures to include:",
                available_structures,
                default=available_structures[:5]  # Limit default selection
            )
        
        if not selected_structures:
            st.warning("Please select at least one structure")
            return
        
        # Report customization
        with st.expander("Report Settings", expanded=True):
            include_metadata = st.checkbox("Include Metadata", value=True)
            include_methodology = st.checkbox("Include Methodology", value=True)
            include_3d_views = st.checkbox("Include 3D Structure Views", value=True)
            include_plots = st.checkbox("Include Analysis Plots", value=True)
            
            # User notes
            report_notes = st.text_area(
                "Additional Notes:",
                placeholder="Add any additional information for the report..."
            )
        
        # Generate report
        if st.button("🚀 Generate Report", type="primary"):
            with st.spinner("Generating report..."):
                try:
                    if report_type == "PDF Report":
                        report_data = self.report_generator.generate_pdf_report(
                            selected_structures,
                            st.session_state.analysis_results,
                            {
                                "include_metadata": include_metadata,
                                "include_methodology": include_methodology,
                                "include_3d_views": include_3d_views,
                                "include_plots": include_plots,
                                "notes": report_notes
                            }
                        )
                        
                        st.download_button(
                            "📥 Download PDF Report",
                            report_data,
                            f"protein_interaction_report_{len(selected_structures)}_structures.pdf",
                            "application/pdf"
                        )
                    
                    elif report_type == "PowerPoint Presentation":
                        pptx_data = self.report_generator.generate_powerpoint_report(
                            selected_structures,
                            st.session_state.analysis_results
                        )
                        
                        st.download_button(
                            "📥 Download PowerPoint",
                            pptx_data,
                            f"protein_interaction_presentation_{len(selected_structures)}_structures.pptx",
                            "application/vnd.openxmlformats-officedocument.presentationml.presentation"
                        )

                    elif report_type == "LaTeX Export":
                        latex_text = self.report_generator.generate_latex_export(
                            selected_structures,
                            st.session_state.analysis_results
                        )
                        st.code(latex_text, language="latex")
                        st.download_button(
                            "📥 Download LaTeX (.tex)",
                            latex_text,
                            f"protein_interaction_report_{len(selected_structures)}_structures.tex",
                            "text/x-tex"
                        )

                    elif report_type == "Complete Package":
                        # Build a comprehensive ZIP containing multiple formats
                        zip_bytes = self.report_generator.generate_complete_package(
                            selected_structures,
                            st.session_state.analysis_results,
                            options={
                                "include_metadata": include_metadata,
                                "include_methodology": include_methodology,
                                "include_3d_views": include_3d_views,
                                "include_plots": include_plots,
                                "notes": report_notes
                            }
                        )
                        st.download_button(
                            "📦 Download Complete Package (ZIP)",
                            zip_bytes,
                            f"protein_interaction_complete_{len(selected_structures)}_structures.zip",
                            "application/zip"
                        )
                    
                    st.success("Report generated successfully!")
                    
                except Exception as e:
                    st.error(f"Report generation failed: {e}")
    
    def _render_settings_tab(self):
        """Render the settings and configuration tab."""
        st.header("⚙️ Settings & Configuration")
        # (Deployment & environment configuration sections removed per request)
        st.caption("Deployment and environment variable settings have been removed to streamline this tab.")
        
        # API Configuration
        st.subheader("🔗 API Integration")
        
        col1, col2 = st.columns(2)
        
        with col1:
            openai_key = st.text_input(
                "OpenAI API Key (optional):",
                type="password",
                value=st.session_state.get("openai_api_key", ""),
                help="For AI-enhanced analysis insights"
            )
            if openai_key:
                st.session_state.openai_api_key = openai_key
        
        with col2:
            anthropic_key = st.text_input(
                "Anthropic API Key (optional):",
                type="password", 
                value=st.session_state.get("anthropic_api_key", ""),
                help="Alternative AI provider for analysis insights"
            )
            if anthropic_key:
                st.session_state.anthropic_api_key = anthropic_key
        
        # Cache Management
        st.subheader("💾 Cache Management")
        
        cache_stats = st.session_state.cache_manager.get_cache_stats()
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Cache Size (MB)", f"{cache_stats.get('cache_size_mb', 0):.1f}")
        with col2:
            st.metric("Cached Items", cache_stats.get('total_items', 0))
        with col3:
            st.metric("PDB Files", cache_stats.get('pdb_files_cached', 0))
        
        if st.button("🗑️ Clear Cache"):
            cleared = st.session_state.cache_manager.clear_cache()
            st.success(f"Cleared {cleared} cache entries")
        
        # Performance Settings
        st.subheader("⚡ Performance")
        st.caption("Select a performance profile or switch to manual to fine-tune individual acceleration features.")
        perf_profile = st.selectbox(
            "Performance Profile",
            options=["auto", "full", "minimal", "manual"],
            index=["auto","full","minimal","manual"].index(st.session_state.get('perf_profile','auto')),
            help="auto: heuristic tuning (recommended); full: enable all; minimal: safe fast-path only; manual: use explicit toggles below"
        )
        st.session_state.perf_profile = perf_profile
        # Persist as env for processing layer to read
        import os as _os
        _os.environ['MOLBRIDGE_PERF_PROFILE'] = perf_profile
        manual = perf_profile == 'manual'
        cols_perf = st.columns(4)
        with cols_perf[0]:
            vec_on = st.checkbox("Vector Geometry", value=True, disabled=not manual, help="Batched NumPy geometry kernels")
        with cols_perf[1]:
            proc_on = st.checkbox("Process Pool", value=False, disabled=not manual, help="Multiprocess detectors + shared memory")
        with cols_perf[2]:
            shm_on = st.checkbox("Shared Memory", value=False, disabled=not manual or not proc_on, help="Share coords & features across processes")
        with cols_perf[3]:
            numba_on = st.checkbox("Numba JIT", value=False, disabled=not manual, help="JIT accelerated distance kernels (fallback if Rust unavailable)")
        cols_perf2 = st.columns(3)
        with cols_perf2[0]:
            task_graph_on = st.checkbox("Task Graph", value=True, disabled=not manual, help="Precompute reusable features once per structure")
        with cols_perf2[1]:
            rust_pref = st.checkbox("Prefer Rust", value=True, disabled=not manual, help="Use PyO3 geometry extension if built")
        with cols_perf2[2]:
            adapt_info = st.checkbox("Emit Adaptive Instrumentation", value=True, disabled=not manual, help="Include adaptive threshold stats in metrics stream")
        if manual:
            if st.button("Apply Manual Performance Flags"):
                _os.environ['MOLBRIDGE_ENABLE_VECTOR_GEOM'] = '1' if vec_on else '0'
                _os.environ['MOLBRIDGE_USE_PROCESS_POOL'] = '1' if proc_on else '0'
                _os.environ['MOLBRIDGE_USE_SHM'] = '1' if shm_on else '0'
                _os.environ['MOLBRIDGE_USE_NUMBA'] = '1' if numba_on else '0'
                _os.environ['MOLBRIDGE_TASK_GRAPH'] = '1' if task_graph_on else '0'
                _os.environ['MOLBRIDGE_USE_RUST'] = '1' if rust_pref else '0'
                st.success("Manual performance flags applied. Re-run analysis to take effect.")
        else:
            st.info(f"Profile '{perf_profile}' will be auto-applied at run time.")
        
        max_workers = st.slider(
            "Max Parallel Workers:",
            min_value=1,
            max_value=8,
            value=self.config.processing.max_workers,
            help="Number of parallel processes for batch analysis"
        )
        
        memory_limit = st.slider(
            "Memory Limit (GB):",
            min_value=1.0,
            max_value=16.0,
            value=self.config.processing.memory_limit_gb,
            step=0.5,
            help="Maximum memory usage for analysis"
        )
        
        # Save settings
        if st.button("💾 Save Settings"):
            # Update config
            self.config.processing.max_workers = max_workers
            self.config.processing.memory_limit_gb = memory_limit
            st.success("Settings saved!")

        # Data output feature toggles (new)
        st.subheader("📦 Data Output Options")
        st.caption("Control normalization, provenance, and experimental columnar storage modes.")
        import os as _dos

        def _env_truth(name: str, default: str = '0') -> bool:
            return _dos.getenv(name, default).lower() in {'1', 'true', 'yes'}

        n_default = _env_truth('MOLBRIDGE_ENABLE_NORMALIZATION')
        p_default = _env_truth('MOLBRIDGE_ENABLE_PROVENANCE')
        c_default = _env_truth('MOLBRIDGE_ENABLE_COLUMNAR')
        dc_default = _env_truth('MOLBRIDGE_DIRECT_COLUMNAR_JSON')

        norm_on = st.checkbox(
            "Normalized Records",
            value=n_default,
            help="Emit unified normalized interaction schema alongside canonical records (env: MOLBRIDGE_ENABLE_NORMALIZATION)."
        )
        prov_on = st.checkbox(
            "Provenance Hashing",
            value=p_default,
            help="Attach reproducibility hash & parameter signature (env: MOLBRIDGE_ENABLE_PROVENANCE)."
        )
        with st.expander("Advanced Data Layout (Columnar Mode)"):
            col_on = st.checkbox(
                "Columnar Storage",
                value=c_default,
                help="Store interactions in columnar arrays to reduce Python object overhead (env: MOLBRIDGE_ENABLE_COLUMNAR)."
            )
            direct_col_on = st.checkbox(
                "Direct Columnar JSON",
                value=dc_default and c_default,
                disabled=not col_on,
                help="Serialize columnar arrays directly to JSON (compact) (env: MOLBRIDGE_DIRECT_COLUMNAR_JSON)."
            )
            if not col_on and dc_default:
                st.info("Direct columnar JSON disabled because Columnar Storage is off.")

        if st.button("Apply Data Output Flags"):
            _dos.environ['MOLBRIDGE_ENABLE_NORMALIZATION'] = '1' if norm_on else '0'
            _dos.environ['MOLBRIDGE_ENABLE_PROVENANCE'] = '1' if prov_on else '0'
            _dos.environ['MOLBRIDGE_ENABLE_COLUMNAR'] = '1' if col_on else '0'
            _dos.environ['MOLBRIDGE_DIRECT_COLUMNAR_JSON'] = '1' if (col_on and direct_col_on) else '0'
            st.success("Data output flags applied. Re-run analysis to take effect.")

        # Advanced Interaction Parameters (auto-generated)
        st.subheader("🧪 Interaction Parameters (All)")
        st.caption("Unified parameter control generated from central registry. Adjust, then re-run analysis.")
        if not hasattr(self, '_interaction_parameter_registry'):
            # Ensure registry exists (user may open Settings before sidebar renders)
            self._render_interaction_strength_filters()
        registry = getattr(self, '_interaction_parameter_registry', {})
        icfg = self.config.interactions
        # Group interactions into logical categories for layout
        groups = {
            'Backbone & Polar': ['hydrogenbond', 'ionicinteraction', 'hydrophobiccontact'],
            'σ-hole': ['halogenbond', 'chalcogenbond', 'pnictogenbond', 'tetrelbond'],
            'π / Aromatic': ['pipi', 'chpi', 'anionpi', 'npistar', 'sulfur_pi', 'cation_pi'],
            'Special / Metal': ['salt_bridge', 'metal_coordination', 'dispersion']
        }
        for gname, ilist in groups.items():
            with st.expander(gname, expanded=False):
                for itype in ilist:
                    params = registry.get(itype, [])
                    if not params:
                        continue
                    st.markdown(f"**{get_interaction_display_names().get(itype, itype)}**")
                    cols = st.columns(min(3, max(1, len(params))))
                    for idx, meta in enumerate(params):
                        with cols[idx % len(cols)]:
                            # Support single-value slider (field) and range slider (fields)
                            if 'fields' in meta:
                                fmin, fmax = meta['fields']
                                cur_min = getattr(icfg, fmin, meta['min'])
                                cur_max = getattr(icfg, fmax, meta['max'])
                                if cur_min > cur_max:
                                    cur_min, cur_max = cur_max, cur_min
                                new_min, new_max = st.slider(meta['label'], meta['min'], meta['max'], (float(cur_min), float(cur_max)), meta['step'], key=f"st_{fmin}_{fmax}", help=meta.get('help'))
                                if new_min != cur_min and hasattr(icfg, fmin):
                                    setattr(icfg, fmin, float(new_min))
                                    st.session_state['parameter_mismatch'] = True
                                if new_max != cur_max and hasattr(icfg, fmax):
                                    setattr(icfg, fmax, float(new_max))
                                    st.session_state['parameter_mismatch'] = True
                            else:
                                field = meta['field']
                                current_val = getattr(icfg, field) if hasattr(icfg, field) else meta['min']
                                new_val = st.slider(meta['label'], meta['min'], meta['max'], float(current_val), meta['step'], key=f"st_{field}", help=meta.get('help'))
                                if new_val != current_val and hasattr(icfg, field):
                                    setattr(icfg, field, float(new_val))
                                    st.session_state['parameter_mismatch'] = True
                    st.markdown("---")
        if st.button("💾 Commit Parameter Changes (All)"):
            st.success("Parameters updated in configuration. Re-run analyses to apply.")
        # Global preset management for full parameter space
        st.markdown("### 📦 Global Parameter Presets")
        gp_cols = st.columns([2,2,2])
        with gp_cols[0]:
            new_preset_name = st.text_input("New Preset Name", key="global_preset_name", placeholder="e.g. hotspot_scan")
        with gp_cols[1]:
            if st.button("Save All As Preset", key="save_global_preset", disabled=not new_preset_name.strip()):
                snapshot = {}
                for pgroup in self._interaction_parameter_registry.values():
                    for meta in pgroup:
                        if 'fields' in meta:
                            f1, f2 = meta['fields']
                            snapshot[f1] = getattr(self.config.interactions, f1, None)
                            snapshot[f2] = getattr(self.config.interactions, f2, None)
                        else:
                            f = meta['field']
                            snapshot[f] = getattr(self.config.interactions, f, None)
                if new_preset_name.strip():
                    self.config.presets[new_preset_name.strip()] = snapshot
                    self.config.applied_preset = new_preset_name.strip()
                    st.success(f"Saved and applied preset '{new_preset_name.strip()}' (re-run analysis to refresh results).")
                    st.session_state['parameter_mismatch'] = False
        with gp_cols[2]:
            available_presets = list(self.config.presets.keys())
            selected_apply = st.selectbox("Apply Preset", options=["(choose)"] + available_presets, key="apply_global_preset_select")
            if st.button("Apply Selected", key="apply_global_preset_btn", disabled=(selected_apply=="(choose)")):
                if selected_apply != "(choose)":
                    pv = self.config.presets.get(selected_apply, {})
                    for k, v in pv.items():
                        if hasattr(self.config.interactions, k):
                            setattr(self.config.interactions, k, v)
                    self.config.applied_preset = selected_apply
                    st.session_state['parameter_mismatch'] = False
                    st.info(f"Applied preset '{selected_apply}'. Re-run analysis to apply to detection.")
    
    def _render_bookmarks_section(self):
        """Render bookmarks and notes section."""
        st.subheader("🔖 Bookmarks & Notes")
        
        # Add bookmark
        with st.expander("Add Bookmark"):
            if st.session_state.current_pdb:
                bookmark_note = st.text_input("Bookmark note:")
                
                if st.button("Add Bookmark"):
                    # This would add a bookmark for the current residue/interaction
                    # Implementation would depend on selection state
                    st.success("Bookmark added!")
        
        # Display existing bookmarks
        bookmarks = st.session_state.session_data.get('bookmarks', [])
        if bookmarks:
            for bookmark in bookmarks:
                with st.container():
                    st.write(f"**{bookmark['pdb_id']}** - {bookmark['note']}")
                    if st.button(f"Remove", key=f"remove_{bookmark['id']}"):
                        st.session_state.session_manager.remove_bookmark(bookmark['id'])
                        st.rerun()
        else:
            st.info("No bookmarks yet")
    
    def _create_interactions_dataframe(self, result: Dict[str, Any]) -> pd.DataFrame:
        """Create a DataFrame from interaction results with individual strength filtering."""
        interactions_data = []

        all_interactions = result.get('interactions', {})
        selected_interactions = [itype for itype, selected in st.session_state.selected_interactions.items() if selected]

        for interaction_type, interactions in all_interactions.items():
            resolved_key = self._resolve_filter_key(interaction_type)

            if resolved_key in selected_interactions:
                # Apply individual strength filtering for this interaction type
                strength_filter = st.session_state.individual_strength_filters.get(resolved_key, 'all')
                filtered_interactions = self._filter_interactions_by_strength(
                    interactions, strength_filter
                )
                
                for interaction in filtered_interactions:
                    # Handle different angle fields for different interaction types
                    angle_info = self._get_angle_display(interaction, interaction_type)
                    
                    # Extract chain information separately for filtering
                    chain1 = self._get_interaction_property(interaction, 'chain1', '')
                    chain2 = self._get_interaction_property(interaction, 'chain2', '')
                    
                    interactions_data.append({
                        "Type": get_interaction_display_names().get(resolved_key, resolved_key.title()),
                        "Residue 1": f"{self._get_interaction_property(interaction, 'residue1', '')} {chain1}",
                        "Residue 2": f"{self._get_interaction_property(interaction, 'residue2', '')} {chain2}",
                        "Chain 1": chain1,
                        "Chain 2": chain2,
                        "Distance (Å)": f"{self._get_interaction_property(interaction, 'distance', 0):.2f}",
                        "Angle (°)": angle_info,
                        "Strength": self._get_interaction_property(interaction, 'strength', 'Moderate')
                    })
        
        return pd.DataFrame(interactions_data)
    
    def _get_interaction_property(self, interaction: Any, property_name: str, default: Any = None) -> Any:
        """Safely get a property from an interaction object or dictionary."""
        if isinstance(interaction, dict):
            return interaction.get(property_name, default)
        
        # Handle different field naming conventions for different interaction types
        if property_name == 'residue1':
            # Try different field names based on interaction type
            for field_name in ['donor_residue', 'halogen_residue', 'chalcogen_residue', 'pnictogen_residue', 
                             'tetrel_residue', 'cation_residue', 'anion_residue', 'residue1']:
                if hasattr(interaction, field_name):
                    return getattr(interaction, field_name, default)
        elif property_name == 'residue2':
            for field_name in ['acceptor_residue', 'residue2']:
                if hasattr(interaction, field_name):
                    return getattr(interaction, field_name, default)
        elif property_name == 'chain1':
            for field_name in ['donor_chain', 'halogen_chain', 'chalcogen_chain', 'pnictogen_chain',
                             'tetrel_chain', 'cation_chain', 'anion_chain', 'chain1']:
                if hasattr(interaction, field_name):
                    return getattr(interaction, field_name, default)
        elif property_name == 'chain2':
            for field_name in ['acceptor_chain', 'chain2']:
                if hasattr(interaction, field_name):
                    return getattr(interaction, field_name, default)
        else:
            # For other properties, try the direct name first
            return getattr(interaction, property_name, default)
        
        return default

    def _get_angle_display(self, interaction: Any, interaction_type: str) -> str:
        """Get formatted angle display for different interaction types."""
        if interaction_type in ['chalcogenbond', 'chalcogen_bond']:
            # Chalcogen bonds have theta and delta angles
            theta = self._get_interaction_property(interaction, 'theta_angle')
            delta = self._get_interaction_property(interaction, 'delta_angle')
            if theta is not None and delta is not None:
                return f"θ:{theta:.1f}° δ:{delta:.1f}°"
            elif theta is not None:
                return f"θ:{theta:.1f}°"

        elif interaction_type in ['halogen_bond', 'halogenbond']:
            # Halogen bonds have C-X...A angle
            angle = self._get_interaction_property(interaction, 'angle')
            if angle is not None:
                return f"{angle:.1f}°"

        elif interaction_type in ['tetrel_bond', 'tetrelbond']:
            # Tetrel bonds have theta1 and theta2 angles
            theta1 = self._get_interaction_property(interaction, 'theta1_angle')
            theta2 = self._get_interaction_property(interaction, 'theta2_angle')
            if theta1 is not None and theta2 is not None:
                return f"θ₁:{theta1:.1f}° θ₂:{theta2:.1f}°"

        elif interaction_type in ['n_pi_star', 'npistar']:
            # n→π* interactions have Bürgi-Dunitz angle
            angle = self._get_interaction_property(interaction, 'angle')
            if angle is not None:
                return f"BD:{angle:.1f}°"

        elif interaction_type in ['hydrogen_bond', 'hbond', 'hydrogenbond']:
            # Hydrogen bonds have D-H...A angle
            angle = self._get_interaction_property(interaction, 'angle')
            if angle is not None:
                return f"DHA:{angle:.1f}°"

        # Generic fallback for other interaction types
        angle = self._get_interaction_property(interaction, 'angle')
        if angle is not None:
            return f"{angle:.1f}°"

        return "N/A"
    
    def _render_session_summary(self):
        """Render summary of the entire session's analysis results."""
        total_structures = len(st.session_state.analysis_results)
        total_interactions = sum(
            sum(len(v) for v in result.get('interactions', {}).values())
            for result in st.session_state.analysis_results.values()
        )
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Structures Analyzed", total_structures)
        
        with col2:
            st.metric("Total Interactions (All)", total_interactions)
        
        with col3:
            avg_interactions = total_interactions / total_structures if total_structures > 0 else 0
            st.metric("Avg. Interactions/Structure", f"{avg_interactions:.1f}")

    def _resolve_filter_key(self, itype: str) -> str:
        """Resolves interaction type key from analysis modules to match session state keys."""
        # Session state keys are clean (e.g., 'chalcogenbond'), from get_interaction_types()
        # Analysis result keys can be different (e.g., 'chalcogen_bond')

        # Direct match
        if itype in st.session_state.individual_strength_filters:
            return itype

        # Common variation: remove underscore
        key_no_underscore = itype.replace('_', '')
        if key_no_underscore in st.session_state.individual_strength_filters:
            return key_no_underscore
        
        # Fallback to original key, which will likely fail to find a filter and use 'all'
        return itype

    def _filter_interactions_by_strength(self, interactions: List[Any], strength_filter: str) -> List[Any]:
        """
        Filter interactions by strength level.

        Args:
            interactions: List of interaction objects or dictionaries
            strength_filter: 'strong', 'strong_moderate', or 'all'

        Returns:
            Filtered list of interactions
        """
        if strength_filter == 'all':
            return interactions

        # Define which strength levels to include
        if strength_filter == 'strong':
            allowed_strengths = ['strong']
        elif strength_filter == 'strong_moderate':
            allowed_strengths = ['strong', 'moderate']
        else:
            return interactions  # Default to all if unknown filter

        # Filter interactions
        filtered = []
        for interaction in interactions:
            strength = self._get_interaction_property(interaction, 'strength', '').lower()

            # Handle different strength naming conventions
            if strength in allowed_strengths:
                filtered.append(interaction)
            # Handle alternative naming (e.g., 'very_weak' -> 'weak')
            elif strength == 'very_weak' and 'weak' in allowed_strengths:
                filtered.append(interaction)
            elif strength == 'very_strong' and 'strong' in allowed_strengths:
                filtered.append(interaction)

        return filtered

    def _apply_individual_strength_filters(self, all_interactions: Dict[str, List[Any]]) -> Dict[str, List[Any]]:
        """
        Apply individual strength filters to each interaction type.

        Args:
            all_interactions: Dictionary mapping interaction types to lists of interactions

        Returns:
            Dictionary with filtered interactions for each type
        """
        filtered_interactions = {}

        # Apply filters and log counts for debugging
        for interaction_type, interactions in all_interactions.items():
            resolved_key = self._resolve_filter_key(interaction_type)
            strength_filter = st.session_state.individual_strength_filters.get(resolved_key, 'all')

            # Apply the filter using the existing filtering function
            filtered = self._filter_interactions_by_strength(interactions, strength_filter)
            filtered_interactions[interaction_type] = filtered

            # Debug info (visible in Streamlit logs) about filtering effect
            try:
                original_count = len(interactions)
                filtered_count = len(filtered)
                logger.info(f"Filter applied for '{interaction_type}' (resolved '{resolved_key}'): {strength_filter} -> {filtered_count}/{original_count} kept")
            except Exception:
                pass

        return filtered_interactions
