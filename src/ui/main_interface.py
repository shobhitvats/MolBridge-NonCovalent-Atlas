"""
Main user interface for Protein Interaction Explorer.
Enhanced with high-performance processing capabilities.
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import zipfile
from pathlib import Path
from typing import List, Dict, Any, Optional
import time
import uuid
import logging
import sys
import os

# Add src to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.config import AppConfig, get_interaction_types, get_interaction_display_names
from utils.pdb_handler import PDBHandler
from analysis.batch_processor import HighPerformanceBatchProcessor
from visualization.structure_viewer import StructureViewer
from visualization.plots import InteractionPlots
from reporting.report_generator import ReportGenerator
from performance.parallel_processor import get_global_processor

# Configure logging for performance monitoring
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class MainInterface:
    """Enhanced main interface controller with high-performance processing."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        self.pdb_handler = PDBHandler(config)
        
        # Use high-performance batch processor
        self.batch_processor = HighPerformanceBatchProcessor(config, use_parallel=True)
        
        self.structure_viewer = StructureViewer(config)
        self.plots = InteractionPlots(config)
        self.report_generator = ReportGenerator(config)
        
        # Get global performance processor for metrics
        self.performance_processor = get_global_processor()
        
        # Initialize session state with performance tracking
        if 'analysis_results' not in st.session_state:
            st.session_state.analysis_results = {}
        if 'current_pdb' not in st.session_state:
            st.session_state.current_pdb = None
        if 'interaction_filters' not in st.session_state:
            st.session_state.interaction_filters = get_interaction_types()
        # Initialize individual strength filters for each interaction type
        if 'individual_strength_filters' not in st.session_state:
            interaction_types = get_interaction_types()
            st.session_state.individual_strength_filters = {
                itype: 'all' for itype in interaction_types
            }
        if 'performance_metrics' not in st.session_state:
            st.session_state.performance_metrics = []
        if 'processing_mode' not in st.session_state:
            st.session_state.processing_mode = 'high_performance'

    def _update_individual_strength_filter(self, itype: str):
        """Callback to update the central strength filter dictionary from the widget's state."""
        widget_key = f"strength_{itype}"
        if widget_key in st.session_state:
            st.session_state.individual_strength_filters[itype] = st.session_state[widget_key]

    def _update_interaction_filter(self, itype: str, widget_key: str):
        """Callback to update the central interaction filter dictionary from the widget's state."""
        if widget_key in st.session_state:
            st.session_state.selected_interactions[itype] = st.session_state[widget_key]

    def _render_interaction_strength_filters(self):
        """Render the integrated interaction and strength filter UI in the sidebar."""
        st.subheader("🔬 Select Interactions & Strength")
        
        interaction_types = get_interaction_types()
        display_names = get_interaction_display_names()

        if 'selected_interactions' not in st.session_state:
            st.session_state.selected_interactions = {itype: False for itype in interaction_types}

        strength_options = {
            'strong': 'Strong',
            'strong_moderate': 'S+M',
            'all': 'All'
        }
        strength_keys = list(strength_options.keys())

        # Quick selection buttons
        col1, col2 = st.columns(2)
        with col1:
            if st.button("Select All", help="Select all interaction types"):
                for itype in interaction_types:
                    st.session_state.selected_interactions[itype] = True
                st.rerun()
        
        with col2:
            if st.button("Clear All", help="Deselect all interaction types"):
                for itype in interaction_types:
                    st.session_state.selected_interactions[itype] = False
                st.rerun()

        # Interaction selection with integrated strength filters
        for itype in interaction_types:
            col1, col2 = st.columns([3, 2])
            with col1:
                widget_key = f"sidebar_interaction_{itype}"
                st.checkbox(
                    display_names[itype],
                    value=st.session_state.selected_interactions.get(itype, False),
                    key=widget_key,
                    on_change=self._update_interaction_filter,
                    args=(itype, widget_key)
                )
            
            if st.session_state.selected_interactions.get(itype, False):
                with col2:
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
                        label_visibility="collapsed"
                    )


    def render(self):
        """Render the main interface with enhanced features."""
        # Accessibility & Theming
        with st.sidebar:
            theme = st.radio("Theme:", ["Light", "Dark"], index=0, help="Switch between light and dark mode for accessibility.")
            st.session_state["theme"] = theme
            st.markdown("<style>body { background-color: #222 if theme=='Dark' else '#fff'; color: #fff if theme=='Dark' else '#222'; }</style>", unsafe_allow_html=True)

        # Sidebar for navigation and settings
        self._render_sidebar()

        # In-app help/tutorial
        with st.expander("🛈 Getting Started & Help", expanded=False):
            st.markdown("""
            **Welcome to the Protein Interaction Explorer!**
            - Use the sidebar to upload or fetch PDB files, select interaction types, and adjust analysis parameters.
            - Switch between tabs for analysis, visualization, results, reports, and settings.
            - Hover over any control for tooltips and guidance.
            - Use the 'Bookmarks & Notes' section to annotate findings.
            - For more help, see the ℹ️ Info tab.
            """)

        # Main content area
        tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
            "🔍 Analysis", 
            "📊 Visualization", 
            "📋 Results", 
            "📄 Reports",
            "🔄 Batch Comparison",
            "⚙️ Settings",
            "ℹ️ Info"
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
            self._render_batch_comparison_tab()

        with tab6:
            self._render_settings_tab()

        with tab7:
            self._render_info_tab()

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
        """Render the info tab with interaction criteria."""
        st.header("ℹ️ Interaction Detection Criteria")
        st.write("This section details the criteria used to detect each type of noncovalent interaction, along with their sources.")

        # Hydrogen Bonds
        with st.expander("Hydrogen Bonds", expanded=True):
            st.write("""
            **Criteria:**
            - **Distance Cutoff:** 3.5 Å (preliminary check)
            - **Angle Cutoff:** 120.0° (D-H...A angle) 
            
            **Source:** IUPAC Recommendations 2011
            **Link:** [doi.org/10.1351/PAC-REP-10-01-01](https://doi.org/10.1351/PAC-REP-10-01-01)
            """)

        # Halogen Bonds
        with st.expander("Halogen Bonds", expanded=True):
            st.write("""
            **Criteria:**
            - **Distance Cutoff:** Sum of van der Waals radii of halogen and acceptor atoms.
            - **Angle Cutoff:** 160.0° (R-X...Y angle)
            
            **Source:** IUPAC Recommendations 2013
            **Link:** [doi.org/10.1351/PAC-REC-12-05-10](https://doi.org/10.1351/PAC-REC-12-05-10)
            """)

        # Pnictogen Bonds
        with st.expander("Pnictogen Bonds", expanded=True):
            st.write("""
            **Criteria:**
            - **Distance Cutoff:** Sum of van der Waals radii of pnictogen and acceptor atoms.
            - **Angle Cutoff:** 150.0° (R-Pn...Y angle)
            
            **Source:** IUPAC Recommendations 2023
            **Link:** [doi.org/10.1515/pac-2020-1002](https://doi.org/10.1515/pac-2020-1002)
            """)

        # Chalcogen Bonds
        with st.expander("Chalcogen Bonds", expanded=True):
            st.write("""
            **Criteria:**
            - Based on geometric parameters derived from statistical analysis of protein structures.
            
            **Source:** Sulfur-mediated chalcogen versus hydrogen bonds in proteins: a see-saw effect in the conformational space
            **Link:** [biorxiv.org/content/10.1101/2023.03.24.534045v1](https://www.biorxiv.org/content/10.1101/2023.03.24.534045v1)
            """)

        # Other Interactions
        with st.expander("Other Interactions (π-π, C-H···π, Anion-π, etc.)", expanded=True):
            st.write("""
            Criteria for other interactions like π-π Stacking, C-H···π, Anion-π, etc., are based on a comprehensive review of unconventional noncovalent interactions.
            
            **Source:** The Realm of Unconventional Noncovalent Interactions in Proteins: Their Significance in Structure and Function
            **Link:** [doi.org/10.1021/acsomega.3c02144](https://doi.org/10.1021/acsomega.3c02144)
            """)
    
    def _render_sidebar(self):
        """Render the sidebar with input options and quick actions."""
        with st.sidebar:
            st.header("🧬 Quick Analysis")
            
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
            
            # Interaction and strength selection
            self._render_interaction_strength_filters()
            
            st.divider()
            
            st.header("📥 Advanced Input")
            
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
            
            # Quick settings
            st.subheader("⚡ Quick Settings")
            
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
            
            # Session management
            st.subheader("💾 Session")
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Save Session", use_container_width=True):
                    session_id = st.session_state.session_manager.save_session()
                    if session_id:
                        st.success("Session saved!")
            
            with col2:
                if st.button("New Session", use_container_width=True):
                    st.session_state.session_data = st.session_state.session_manager._create_empty_session()
                    st.session_state.analysis_results = {}
                    st.rerun()
            
            st.divider()
            
            # Performance monitoring section
            st.subheader("⚡ Performance Monitor")
            
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
                use_container_width=True,
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
                    
                    if st.button("🚀 Analyze Batch", type="primary", use_container_width=True, disabled=(selected_count == 0)):
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
            
            if st.button("Load Session", use_container_width=True):
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
        
        try:
            # Step 1: Download/load structure
            status_text.text(f"Loading structure {pdb_id}...")
            progress_bar.progress(0.2)
            
            structure = self.pdb_handler.load_structure(
                pdb_id, 
                assembly=st.session_state.structure_settings["assembly"]
            )
            
            if not structure:
                st.error(f"Failed to load structure {pdb_id}")
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
            
            # Use enhanced processor for single protein with selected interactions
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
                        st.json(perf_metrics)
            else:
                st.error(f"❌ Analysis failed for {pdb_id}: {results.get('error', 'Unknown error')}")
            
            perf_info.empty()  # Clear the performance info
            
        except Exception as e:
            st.error(f"Analysis failed: {e}")
        finally:
            progress_bar.empty()
            status_text.empty()
    
    def _run_batch_analysis(self, pdb_ids: List[str], selected_interactions: List[str] = None):
        """Run high-performance batch analysis for multiple PDB structures with selected interactions."""
        if selected_interactions is None:
            selected_interactions = get_interaction_types()
            
        batch_id = str(uuid.uuid4())
        
        # Create progress tracking
        progress_bar = st.progress(0)
        status_text = st.empty()
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
                        st.dataframe(perf_df, use_container_width=True)
                
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
        
        # Re-analysis controls
        st.write("---")
        col1, col2 = st.columns(2)
        
        with col1:
            # Re-analysis button
            if st.button("🔄 Re-analyze with Current Settings", type="primary", use_container_width=True):
                if st.session_state.current_pdb:
                    selected_interactions = [itype for itype, selected in st.session_state.selected_interactions.items() if selected]
                    if selected_interactions:
                        with st.spinner("Re-analyzing with new parameters..."):
                            self._run_single_analysis(st.session_state.current_pdb, selected_interactions)
                        st.success("✅ Re-analysis complete!")
                    else:
                        st.error("❌ Please select at least one interaction type")
                else:
                    st.error("❌ No structure selected")
        
        with col2:
            # Apply filters button (for real-time filtering without re-analysis)
            if st.button("🎯 Apply Filters Only", type="secondary", use_container_width=True):
                st.success("✅ Filters applied! Check Visualization and Results tabs.")
                st.rerun()
    
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
        
        # Initialize parameter values in session state if not exists
        if 'custom_parameters' not in st.session_state:
            st.session_state.custom_parameters = {}
        
        # Organize parameters by interaction type
        parameter_tabs = st.tabs([
            "🔹 Basic Interactions", 
            "🔸 σ-hole Interactions", 
            "🔺 π-System Interactions",
            "🔻 Other Interactions"
        ])
        
        with parameter_tabs[0]:
            st.write("**Basic Non-Covalent Interactions**")
            col1, col2 = st.columns(2)
            
            with col1:
                # Hydrogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('hydrogenbond', True):
                    st.write("**Hydrogen Bonds**")
                    st.session_state.custom_parameters["hbond_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.0, max_value=5.0, value=3.5, step=0.1,
                        key="hbond_distance", help="Maximum donor-acceptor distance"
                    )
                    st.session_state.custom_parameters["hbond_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=90, max_value=180, value=120, step=5,
                        key="hbond_angle", help="Minimum donor-H-acceptor angle"
                    )
                
                # Ionic Interactions
                st.write("**Ionic Interactions**")
                st.session_state.custom_parameters["ionic_distance"] = st.slider(
                    "Distance Cutoff (Å)",
                    min_value=3.0, max_value=8.0, value=6.0, step=0.1,
                    key="ionic_distance", help="Maximum charge-charge distance"
                )
            
            with col2:
                # Hydrophobic Contacts (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('hydrophobiccontact', True):
                    st.write("**Hydrophobic Contacts**")
                    st.session_state.custom_parameters["hydrophobic_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=6.0, value=5.0, step=0.1,
                        key="hydrophobic_distance", help="Maximum carbon-carbon distance"
                    )
                
                # London Dispersion (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('dispersion', True):
                    st.write("**London Dispersion**")
                    st.session_state.custom_parameters["dispersion_min"] = st.slider(
                        "Min Distance (Å)",
                        min_value=3.0, max_value=4.0, value=3.5, step=0.1,
                        key="dispersion_min", help="Minimum contact distance"
                    )
                    st.session_state.custom_parameters["dispersion_max"] = st.slider(
                        "Max Distance (Å)",
                        min_value=4.0, max_value=6.0, value=5.0, step=0.1,
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
                        min_value=2.5, max_value=5.0, value=4.0, step=0.1,
                        key="halogen_distance", help="Maximum X···A distance"
                    )
                    st.session_state.custom_parameters["halogen_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=120, max_value=180, value=140, step=5,
                        key="halogen_angle", help="Minimum C-X···A angle"
                    )
                
                # Chalcogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('chalcogenbond', True):
                    st.write("**Chalcogen Bonds**")
                    st.session_state.custom_parameters["chalcogen_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=5.0, value=4.0, step=0.1,
                        key="chalcogen_distance", help="Maximum Ch···A distance"
                    )
                    st.session_state.custom_parameters["chalcogen_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=120, max_value=180, value=140, step=5,
                        key="chalcogen_angle", help="Minimum C-Ch···A angle"
                    )
            
            with col2:
                # Pnictogen Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('pnictogenbond', True):
                    st.write("**Pnictogen Bonds**")
                    st.session_state.custom_parameters["pnictogen_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=5.0, value=4.0, step=0.1,
                        key="pnictogen_distance", help="Maximum Pn···A distance"
                    )
                    st.session_state.custom_parameters["pnictogen_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=120, max_value=180, value=140, step=5,
                        key="pnictogen_angle", help="Minimum C-Pn···A angle"
                    )
                
                # Tetrel Bonds (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('tetrelbond', True):
                    st.write("**Tetrel Bonds**")
                    st.session_state.custom_parameters["tetrel_distance_min"] = st.slider(
                        "Min Distance (Å)",
                        min_value=2.0, max_value=3.0, value=2.5, step=0.1,
                        key="tetrel_distance_min", help="Minimum T···A distance"
                    )
                    st.session_state.custom_parameters["tetrel_distance_max"] = st.slider(
                        "Max Distance (Å)",
                        min_value=3.0, max_value=4.5, value=3.6, step=0.1,
                        key="tetrel_distance_max", help="Maximum T···A distance"
                    )
                    st.session_state.custom_parameters["tetrel_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=140, max_value=180, value=160, step=5,
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
                        min_value=3.0, max_value=8.0, value=5.5, step=0.1,
                        key="pi_pi_distance", help="Maximum centroid-centroid distance"
                    )
                    st.session_state.custom_parameters["pi_pi_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=0, max_value=45, value=30, step=5,
                        key="pi_pi_angle", help="Maximum dihedral angle between rings"
                    )

                # C-H···π Interactions (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('chpi', True):
                    st.write("**C-H···π Interactions**")
                    st.session_state.custom_parameters["ch_pi_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=6.0, value=4.5, step=0.1,
                        key="ch_pi_distance", help="Maximum C···centroid distance"
                    )
                    st.session_state.custom_parameters["ch_pi_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=60, max_value=120, value=90, step=5,
                        key="ch_pi_angle", help="Maximum C-H···centroid angle"
                    )
            
            with col2:
                # Anion-π Interactions (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('anionpi', True):
                    st.write("**Anion-π Interactions**")
                    st.session_state.custom_parameters["anion_pi_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=3.0, max_value=7.0, value=5.0, step=0.1,
                        key="anion_pi_distance", help="Maximum anion-centroid distance"
                    )

                # n→π* Interactions (render only if selected)
                if st.session_state.get('selected_interactions', {}).get('npistar', True):
                    st.write("**n→π* Interactions**")
                    st.session_state.custom_parameters["n_pi_star_distance"] = st.slider(
                        "Distance Cutoff (Å)",
                        min_value=2.5, max_value=4.5, value=3.5, step=0.1,
                        key="n_pi_star_distance", help="Maximum n···C=O distance"
                    )
                    st.session_state.custom_parameters["n_pi_star_angle"] = st.slider(
                        "Angle Cutoff (°)",
                        min_value=90, max_value=150, value=120, step=5,
                        key="n_pi_star_angle", help="Minimum n···C=O angle"
                    )
        
        with parameter_tabs[3]:
            st.write("**Specialized Parameters**")
            col1, col2 = st.columns(2)
            
            with col1:
                st.write("**Quality Thresholds**")
                st.session_state.custom_parameters["strong_threshold"] = st.slider(
                    "Strong Interaction Threshold",
                    min_value=0.7, max_value=0.95, value=0.85, step=0.05,
                    key="strong_threshold", help="Score threshold for strong classification"
                )
                st.session_state.custom_parameters["moderate_threshold"] = st.slider(
                    "Moderate Interaction Threshold", 
                    min_value=0.5, max_value=0.8, value=0.65, step=0.05,
                    key="moderate_threshold", help="Score threshold for moderate classification"
                )
            
            with col2:
                st.write("**Environmental Factors**")
                st.session_state.custom_parameters["exclude_intraresidue"] = st.checkbox(
                    "Exclude Intraresidue Interactions",
                    value=False,
                    key="exclude_intraresidue", 
                    help="Skip interactions within the same residue"
                )
                st.session_state.custom_parameters["min_sequence_separation"] = st.slider(
                    "Min Sequence Separation",
                    min_value=0, max_value=5, value=1, step=1,
                    key="min_sequence_separation", 
                    help="Minimum residue separation in sequence"
                )
        
        # Parameter reset options
        st.write("---")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("🔄 Reset to Defaults", help="Reset all parameters to default values"):
                st.session_state.custom_parameters = {}
                st.rerun()
        
        with col2:
            if st.button("📊 Conservative Preset", help="Set conservative parameters for high confidence"):
                st.session_state.custom_parameters.update({
                    "strong_threshold": 0.9, "moderate_threshold": 0.75,
                    "hbond_distance": 3.2, "halogen_distance": 3.8,
                    "pi_pi_distance": 5.0, "ionic_distance": 5.5
                })
                st.rerun()
        
        with col3:
            if st.button("🔍 Exploratory Preset", help="Set liberal parameters for discovery"):
                st.session_state.custom_parameters.update({
                    "strong_threshold": 0.7, "moderate_threshold": 0.5,
                    "hbond_distance": 4.0, "halogen_distance": 4.5,
                    "pi_pi_distance": 6.0, "ionic_distance": 7.0
                })
                st.rerun()

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

    def _render_analysis_summary(self):
        """Render the new point-wise analysis summary for the current protein."""
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

        # Accessibility: ARIA label for main visualization
        st.markdown('<div aria-label="Protein 3D Structure Visualization" tabindex="0"></div>', unsafe_allow_html=True)

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

        # External database links
        pdb_id = st.session_state.current_pdb
        st.markdown(f"[� View {pdb_id} on RCSB PDB](https://www.rcsb.org/structure/{pdb_id}) | [UniProt Search](https://www.uniprot.org/uniprot/?query={pdb_id}) | [Literature](https://pubmed.ncbi.nlm.nih.gov/?term={pdb_id})")

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

    

    def _render_results_tab(self):
        """Render the results summary tab with annotation and sharing features."""
        if not st.session_state.analysis_results:
            st.info("No analysis results available")
            return

        st.header("📋 Analysis Results")

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
            st.subheader(f"Interactions in {st.session_state.current_pdb}")
            interactions_df = self._create_interactions_dataframe(result)
            if not interactions_df.empty:
                total_interactions = sum(len(interactions) for interactions in result.get('interactions', {}).values())
                filtered_count = len(interactions_df)
                if active_filters:
                    st.metric("Filtered Interactions Shown", f"{filtered_count}", f"of {total_interactions} total")
                st.dataframe(interactions_df, use_container_width=True)
                # Export, copy, and share options
                col1, col2, col3, col4 = st.columns(4)
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
            else:
                st.warning("No interactions found with current filters")

        # Bookmarks and notes
        self._render_bookmarks_section()
    
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
                    
                    st.success("Report generated successfully!")
                    
                except Exception as e:
                    st.error(f"Report generation failed: {e}")
    
    def _render_settings_tab(self):
        """Render the settings and configuration tab."""
        st.header("⚙️ Settings & Configuration")
        
        # Deployment Section
        st.subheader("🚀 Deployment & Sharing")
        
        with st.expander("🌐 Deploy Online", expanded=False):
            st.write("**Deploy your Protein Interaction Analysis Server online for easy sharing and access.**")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.markdown("### Streamlit Cloud")
                st.write("**Free & Easy**")
                st.write("Perfect for sharing with collaborators")
                
                if st.button("🚀 Deploy to Streamlit Cloud", type="primary", use_container_width=True):
                    st.markdown("""
                    **Follow these steps:**
                    1. Go to [share.streamlit.io](https://share.streamlit.io)
                    2. Connect your GitHub account
                    3. Select this repository
                    4. Set main file path to: `server.py`
                    5. Click **Deploy**
                    
                    Your app will be live at: `https://your-app-name.streamlit.app`
                    """, unsafe_allow_html=True)
                    st.info("💡 Make sure your repository is public for free deployment!")
                
                st.markdown("[📖 Streamlit Cloud Docs](https://docs.streamlit.io/streamlit-cloud)")
            
            with col2:
                st.markdown("### Heroku")
                st.write("**Professional Hosting**")
                st.write("Great for production use")
                
                if st.button("🐘 Deploy to Heroku", use_container_width=True):
                    st.markdown("""
                    **Quick Setup:**
                    1. Create account at [heroku.com](https://heroku.com)
                    2. Install Heroku CLI
                    3. Run in terminal:
                    ```bash
                    heroku create your-app-name
                    git push heroku main
                    ```
                    
                    **Files needed:**
                    - `Procfile`: `web: streamlit run server.py --server.port $PORT`
                    - `requirements.txt`: Your dependencies
                    """, unsafe_allow_html=True)
                
                st.markdown("[📖 Heroku Deployment Guide](https://devcenter.heroku.com/articles/getting-started-with-python)")
            
            with col3:
                st.markdown("### Docker + Cloud")
                st.write("**Full Control**")
                st.write("Deploy anywhere with Docker")
                
                if st.button("🐳 Docker Setup", use_container_width=True):
                    st.markdown("""
                    **Dockerfile:**
                    ```dockerfile
                    FROM python:3.12-slim
                    WORKDIR /app
                    COPY requirements.txt .
                    RUN pip install -r requirements.txt
                    COPY . .
                    CMD ["streamlit", "run", "server.py", "--server.address", "0.0.0.0"]
                    ```
                    
                    **Deploy to:**
                    - AWS, Google Cloud, Azure
                    - DigitalOcean, Railway
                    - Your own server
                    """, unsafe_allow_html=True)
                
                st.markdown("[📖 Docker Deployment](https://docs.docker.com/get-started/)")
        
        # Quick Deployment Status
        st.write("---")
        st.subheader("📊 Deployment Status")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("🔍 Check Repository Status", use_container_width=True):
                st.info("**Repository Check:**")
                st.write("✅ Main file: `server.py`")
                st.write("✅ Dependencies: `requirements.txt`")
                st.write("✅ Source code: Organized in `src/`")
                st.write("✅ Ready for deployment!")
        
        with col2:
            if st.button("📋 Generate Deployment Checklist", use_container_width=True):
                st.success("**Deployment Checklist Generated:**")
                checklist = """
                □ Repository is public on GitHub
                □ Main file is `server.py`
                □ `requirements.txt` includes all dependencies
                □ App runs locally with `streamlit run server.py`
                □ No hardcoded local paths
                □ Sensitive data moved to environment variables
                □ README.md includes setup instructions
                """
                st.code(checklist, language="text")
        
        # Environment Variables Setup
        st.write("---")
        st.subheader("🔐 Environment Variables")
        
        with st.expander("Configure Environment Variables", expanded=False):
            st.write("**Set up environment variables for your deployment:**")
            
            env_vars = {
                "STREAMLIT_SERVER_PORT": "8501",
                "STREAMLIT_SERVER_ADDRESS": "0.0.0.0",
                "STREAMLIT_SERVER_HEADLESS": "true",
                "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false"
            }
            
            for var, default in env_vars.items():
                value = st.text_input(f"{var}", value=default, help=f"Default: {default}")
                if value != default:
                    st.code(f"export {var}={value}")
            
            st.info("💡 Copy these commands to your deployment platform's environment variables section.")
        
        st.write("---")
        
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
                    
                    interactions_data.append({
                        "Type": get_interaction_display_names().get(resolved_key, resolved_key.title()),
                        "Residue 1": f"{self._get_interaction_property(interaction, 'residue1', '')} {self._get_interaction_property(interaction, 'chain1', '')}",
                        "Residue 2": f"{self._get_interaction_property(interaction, 'residue2', '')} {self._get_interaction_property(interaction, 'chain2', '')}",
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
