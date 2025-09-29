"""
Analysis plots and visualizations for Protein Interaction Explorer.
"""

import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, Any, List, Optional
import os
import json
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from utils.config import AppConfig
from ui.components.color_legend import COLORBLIND_PALETTE, normalize_key

class InteractionPlots:
    """Generates plots and visualizations for interaction analysis."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        # Attempt to map configured order onto unified palette; fallback to palette values
        try:
            conf_map = config.visualization.interaction_colors
            self.color_palette = [
                COLORBLIND_PALETTE.get(normalize_key(k), v)
                for k, v in conf_map.items()
            ]
        except Exception:
            self.color_palette = list(COLORBLIND_PALETTE.values())
    
    def _get_interaction_property(self, interaction: Any, property_name: str, default_value: Any = None) -> Any:
        """
        Get a property from an interaction object, handling both dictionary and object styles.
        
        Args:
            interaction: The interaction object (dict or dataclass)
            property_name: Name of the property to retrieve
            default_value: Default value if property is not found
            
        Returns:
            The property value or default_value
        """
        if hasattr(interaction, property_name):
            # Object-style interaction (e.g., HydrogenBond dataclass)
            return getattr(interaction, property_name, default_value)
        elif isinstance(interaction, dict):
            # Dictionary-style interaction
            return interaction.get(property_name, default_value)
        else:
            # Unknown format, return default
            return default_value
    
    def render_chain_heatmap(self, analysis_result: Dict[str, Any]):
        """Render chain-vs-chain interaction heatmap with safer Lottie fallback."""
        st.subheader("🔥 Chain Interaction Heatmap")
        chain_matrix = self._create_chain_interaction_matrix(analysis_result)
        if chain_matrix.empty:
            try:
                from streamlit_lottie import st_lottie
                import requests
                import time
                unique_key = f"empty_heatmap_{int(time.time() * 1000)}"
                animation_data = requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json", timeout=5).json()
                st_lottie(animation_data, height=180, key=unique_key)
            except Exception:
                st.info("🔍 No interactions found to display")
            st.warning("No inter-chain interactions found")
            return
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(
                chain_matrix,
                annot=True,
                fmt='d',
                cmap='YlOrRd',
                ax=ax,
                cbar_kws={'label': 'Number of Interactions'}
            )
            ax.set_title('Inter-Chain Interaction Count')
            ax.set_xlabel('Chain')
            ax.set_ylabel('Chain')
            st.pyplot(fig)
            plt.close()
        except Exception as e:
            try:
                from streamlit_lottie import st_lottie
                import requests
                import time
                unique_key = f"error_heatmap_{int(time.time() * 1000)}"
                animation_data = requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json", timeout=5).json()
                st_lottie(animation_data, height=180, key=unique_key)
            except Exception:
                st.error("⚠️ Error generating visualization")
            st.error(f"Error rendering heatmap: {e}")
            return
        if st.checkbox("Show Interactive Heatmap"):
            try:
                self._render_interactive_heatmap(chain_matrix)
            except Exception as e:
                try:
                    from streamlit_lottie import st_lottie
                    import requests
                    import time
                    unique_key = f"error_interactive_heatmap_{int(time.time() * 1000)}"
                    animation_data = requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json", timeout=5).json()
                    st_lottie(animation_data, height=180, key=unique_key)
                except Exception:
                    st.error("⚠️ Error generating interactive visualization")
                st.error(f"Error rendering interactive heatmap: {e}")
    
    def _create_chain_interaction_matrix(self, analysis_result: Dict[str, Any]) -> pd.DataFrame:
        """Create chain interaction matrix."""
        interactions = analysis_result.get('interactions', {})
        
        # Get all chains
        chains = set()
        chain_pairs = []
        
        for interaction_type, interaction_list in interactions.items():
            for interaction in interaction_list:
                chain1 = self._get_interaction_property(interaction, 'chain1', '')
                chain2 = self._get_interaction_property(interaction, 'chain2', '')
                
                if chain1 and chain2 and chain1 != chain2:
                    chains.add(chain1)
                    chains.add(chain2)
                    chain_pairs.append((chain1, chain2))
        
        if not chains:
            return pd.DataFrame()
        
        # Create matrix
        chains = sorted(list(chains))
        matrix = pd.DataFrame(0, index=chains, columns=chains)
        
        for chain1, chain2 in chain_pairs:
            matrix.loc[chain1, chain2] += 1
            matrix.loc[chain2, chain1] += 1  # Make symmetric
        
        return matrix
    
    def _render_interactive_heatmap(self, chain_matrix: pd.DataFrame):
        """Render interactive heatmap with Plotly."""
        fig = px.imshow(
            chain_matrix,
            labels=dict(x="Chain", y="Chain", color="Interactions"),
            x=chain_matrix.columns,
            y=chain_matrix.index,
            color_continuous_scale="YlOrRd",
            title="Interactive Chain Interaction Heatmap"
        )
        
        fig.update_layout(
            width=600,
            height=600
        )
        
        st.plotly_chart(fig, width="stretch")
    
    def render_interaction_distribution(self, analysis_result: Dict[str, Any]):
        """Render interaction type distribution plots with safer Lottie fallback."""
        st.subheader("📊 Interaction Distribution")
        interactions = analysis_result.get('interactions', {})
        interaction_counts = {
            self._get_display_name(int_type): len(int_list)
            for int_type, int_list in interactions.items()
            if int_list
        }
        if not interaction_counts:
            try:
                from streamlit_lottie import st_lottie
                import requests
                import time
                unique_key = f"empty_dist_{int(time.time() * 1000)}"
                animation_data = requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json", timeout=5).json()
                st_lottie(animation_data, height=180, key=unique_key)
            except Exception:
                st.info("📊 No distribution data available")
            st.warning("No interactions found")
            return
        viz_type = st.radio(
            "Visualization Type:",
            ["Bar Chart", "Pie Chart", "Distance Distribution"],
            horizontal=True
        )
        try:
            if viz_type == "Bar Chart":
                self._render_bar_chart(interaction_counts)
            elif viz_type == "Pie Chart":
                self._render_pie_chart(interaction_counts)
            elif viz_type == "Distance Distribution":
                self._render_distance_distribution(interactions)
        except Exception as e:
            try:
                from streamlit_lottie import st_lottie
                import requests
                import time
                unique_key = f"error_dist_{int(time.time() * 1000)}"
                animation_data = requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json", timeout=5).json()
                st_lottie(animation_data, height=180, key=unique_key)
            except Exception:
                st.error("⚠️ Error generating distribution chart")
            st.error(f"Error rendering distribution: {e}")
    
    def _render_bar_chart(self, interaction_counts: Dict[str, int]):
        """Render bar chart of interaction counts using unified palette."""
        keys = list(interaction_counts.keys())
        colors = [COLORBLIND_PALETTE.get(normalize_key(k), '#999999') for k in keys]
        fig = go.Figure(data=[go.Bar(x=keys, y=[interaction_counts[k] for k in keys], marker_color=colors)])
        
        fig.update_layout(
            title="Interaction Type Distribution",
            xaxis_title="Interaction Type",
            yaxis_title="Count",
            xaxis_tickangle=-45
        )
        
        st.plotly_chart(fig, width="stretch")
    
    def _render_pie_chart(self, interaction_counts: Dict[str, int]):
        """Render pie chart of interaction counts with consistent palette."""
        keys = list(interaction_counts.keys())
        colors = [COLORBLIND_PALETTE.get(normalize_key(k), '#999999') for k in keys]
        fig = go.Figure(data=[go.Pie(labels=keys, values=[interaction_counts[k] for k in keys], hole=0.3, marker_colors=colors)])
        fig.update_layout(title="Interaction Type Distribution")
        st.plotly_chart(fig, width="stretch")
    
    def _render_distance_distribution(self, interactions: Dict[str, List[Dict]]):
        """Render distance distribution for different interaction types."""
        # Collect distance data
        distance_data = []
        
        for interaction_type, interaction_list in interactions.items():
            for interaction in interaction_list:
                if self._get_interaction_property(interaction, 'distance', None) is not None:
                    distance_data.append({
                        'Interaction Type': self._get_display_name(interaction_type),
                        'Distance (Å)': self._get_interaction_property(interaction, 'distance', 0)
                    })
        
        if not distance_data:
            st.warning("No distance data available")
            return
        
        df = pd.DataFrame(distance_data)
        
        # Create violin plot
        fig = px.violin(
            df,
            x='Interaction Type',
            y='Distance (Å)',
            box=True,
            title="Distance Distribution by Interaction Type",
            color='Interaction Type',
            color_discrete_map={k: COLORBLIND_PALETTE.get(normalize_key(k), '#999999') for k in df['Interaction Type'].unique()}
        )
        
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")
    
    def render_ramachandran_plot(self, analysis_result: Dict[str, Any]):
        """Render Ramachandran plot with interaction highlights."""
        st.subheader("📐 Ramachandran Plot")
        # Always use an expander for controls (compact by default)
        with st.expander("Plot Controls", expanded=False):
            ui_col1, ui_col2, ui_col3, ui_col4 = st.columns([1,1,1,1])
            with ui_col1:
                show_density = st.checkbox("Density", value=st.session_state.get('rama_show_density', True), key='rama_show_density', help="Show KDE-like density background")
            with ui_col2:
                hide_polygons = st.checkbox("Hide Polygons", value=st.session_state.get('rama_hide_polygons', True), key='rama_hide_polygons', help="Hide favored/allowed region boundary polygons.")
            with ui_col3:
                filter_outliers = st.checkbox("Hide Outliers", value=False)
            with ui_col4:
                hotspots_only = st.checkbox("Hotspots Only", value=False)
            # Multi-structure overlay toggle
            overlay_row1, overlay_row2 = st.columns([1,3])
            with overlay_row1:
                multi_overlay = st.checkbox("Multi-Overlay", value=st.session_state.get('rama_multi_overlay', False), key='rama_multi_overlay', help="Overlay phi/psi from other analyzed structures (basic union; current structure highlighted)")
            selected_overlay_ids = []
            if multi_overlay:
                # Gather available structures from session (excluding current if present)
                all_results = getattr(st.session_state, 'analysis_results', {}) or {}
                current_id = analysis_result.get('pdb_id') or st.session_state.get('current_pdb')
                other_ids = [pid for pid in all_results.keys() if pid != current_id]
                if other_ids:
                    with overlay_row2:
                        selected_overlay_ids = st.multiselect("Structures", other_ids, default=other_ids[:3], help="Select structures to overlay (limited to first 5)")[:5]
                else:
                    st.caption("No additional analyzed structures for overlay yet.")
            class_mode = st.radio(
                "Residue Classes",
                ["All","Custom"],
                horizontal=True,
                help="Choose whether to display all residue classes or select a subset.",
                key='rama_class_mode'
            )
            if class_mode == "Custom":
                class_filter = st.multiselect(
                    "Select classes",
                    ["General","Pre-Pro","Gly","Pro"],
                    default=st.session_state.get('rama_class_filter', ["General","Pre-Pro","Gly","Pro"]),
                    help="Classes included in plot (data + stats).",
                    key='rama_class_filter'
                )
            else:
                class_filter = ["General","Pre-Pro","Gly","Pro"]
            adv_col1, adv_col2, adv_col3, adv_col4 = st.columns([1,1,1,1])
            with adv_col1:
                outlier_focus = st.checkbox("Outlier Focus", value=st.session_state.get('rama_outlier_focus', False), key='rama_outlier_focus', help="Emphasize outliers with halo and dim others")
            with adv_col2:
                interaction_color_scale = st.checkbox("Color by Interactions", value=st.session_state.get('rama_color_by_interactions', False), key='rama_color_by_interactions', help="Override region colors with a color scale based on interaction count")
            with adv_col3:
                density_tune = st.checkbox("Density Settings", value=st.session_state.get('rama_density_tune', False), key='rama_density_tune', help="Show advanced bandwidth/grid controls and optional quantile contours")
            with adv_col4:
                hide_density_scale = st.checkbox("Hide Density Scale", value=st.session_state.get('rama_hide_density_scale', False), key='rama_hide_density_scale', help="Hide the density colorbar even when density layer is shown.")
            legend_col, reset_col = st.columns([4,1])
            with legend_col:
                minimal_legend = st.checkbox("Minimal Legend", value=st.session_state.get('rama_minimal_legend', False), key='rama_minimal_legend', help="Hide size reference points; keep only region color entries.")
            with reset_col:
                if st.button("Reset", help="Restore default plot control settings"):
                    for k in [
                        'rama_show_density','rama_hide_polygons','rama_class_mode','rama_class_filter',
                        'rama_outlier_focus','rama_color_by_interactions','rama_density_tune','rama_hide_density_scale',
                        'rama_minimal_legend','rama_overlay_mode','rama_overlay_classes'
                    ]:
                        if k in st.session_state: del st.session_state[k]
                    st.experimental_rerun()
            # Optional region focus to explicitly isolate without relying on legend gestures
            region_focus = st.radio(
                "Region Focus",
                ["All","Favored","Allowed","Outlier"],
                index=["All","Favored","Allowed","Outlier"].index(st.session_state.get('rama_region_focus','All')),
                horizontal=True,
                key='rama_region_focus',
                help="Filter points and density by a single Ramachandran region."
            )
            # Optional: generate a sharable link that encodes current control state
            if st.button("Copy Sharable Link"):
                # Collect a minimal set of params
                params = {
                    'density': str(int(show_density)),
                    'hidePoly': str(int(hide_polygons)),
                    'classMode': class_mode,
                    'classes': ",".join(class_filter),
                    'outlierFocus': str(int(outlier_focus)),
                    'colorByInter': str(int(interaction_color_scale)),
                    'regionFocus': region_focus,
                    'minLegend': str(int(minimal_legend)),
                    'hideDensityScale': str(int(hide_density_scale))
                }
                st.experimental_set_query_params(**params)
                st.info("URL updated with current controls. Use your browser to copy the address.")
        bw_method = None
        grid_points = None
        quantile_contours = []
        if density_tune and show_density:
            with st.expander("Density Parameters"):
                bw_method = st.slider("KDE Bandwidth Multiplier", 0.15, 0.75, 0.25, 0.05)
                grid_points = st.slider("Grid Resolution", 60, 220, 160, 20)
                qc_on = st.checkbox("Quantile Contours", value=False, help="Overlay highest-density regions at selected probability mass levels.")
                if qc_on:
                    picked = st.multiselect("Contour levels (percent)", [50,90,99], default=[50,90,99])
                    quantile_contours = sorted({p/100 for p in picked})
        rama_df = None
        cache_key = None
        try:
            # Attempt to use stored structure from session state if available
            structure = getattr(st.session_state, 'rama_source_structure', None)
            if structure is not None:
                # Build a cache key from chain and residue count
                chains = [chain.id for model in structure for chain in model]
                residue_total = sum(len(chain) for model in structure for chain in model)
                cache_key = f"rama|{'-'.join(sorted(set(chains)))}|{residue_total}"
                if 'rama_cache' not in st.session_state:
                    st.session_state.rama_cache = {}
                if cache_key in st.session_state.rama_cache:
                    rama_df = st.session_state.rama_cache[cache_key]
                else:
                    rama_df = self._compute_phi_psi_dataframe(structure, analysis_result)
                    if rama_df is not None and not rama_df.empty:
                        st.session_state.rama_cache[cache_key] = rama_df
            if rama_df is None or rama_df.empty:
                st.info("Using placeholder Ramachandran data (structure angles unavailable)")
                rama_df = self._placeholder_rama_dataframe(analysis_result)
        except Exception as e:
            st.warning(f"Ramachandran extraction failed ({e}); using placeholder data")
            rama_df = self._placeholder_rama_dataframe(analysis_result)
        
        # Multi-structure overlay merging
        overlay_df = None
        if multi_overlay and selected_overlay_ids:
            overlay_rows = []
            # Lazy import to avoid upstream circular imports
            try:
                from utils.pdb_handler import PDBHandler  # type: ignore
                pdb_loader = PDBHandler(self.config)
            except Exception:
                pdb_loader = None
            for pid in selected_overlay_ids:
                base_res = getattr(st.session_state, 'analysis_results', {}).get(pid, {})
                existing = base_res.get('ramachandran_df')
                if existing is not None:
                    df_pid = existing.copy()
                else:
                    structure_variant = None
                    if pdb_loader is not None:
                        try:
                            structure_variant = pdb_loader.fetch_structure_variant(pid, assembly=self.config.default_assembly, include_ligands=self.config.include_ligands, exclude_waters=self.config.exclude_waters)
                        except Exception:
                            structure_variant = None
                    if structure_variant is not None:
                        try:
                            df_pid = self._compute_phi_psi_dataframe(structure_variant, base_res)
                        except Exception:
                            df_pid = None
                    else:
                        df_pid = None
                if df_pid is None or df_pid.empty:
                    continue
                df_pid = df_pid[['Phi','Psi','InteractionCount','Type','Chain','ResSeq','Residue']].copy()
                df_pid['SourcePDB'] = pid
                overlay_rows.append(df_pid)
            if overlay_rows:
                overlay_df = pd.concat(overlay_rows, ignore_index=True)
        # Mark current structure source for clarity
        if rama_df is not None and not rama_df.empty and multi_overlay:
            rama_df['SourcePDB'] = analysis_result.get('pdb_id') or st.session_state.get('current_pdb') or 'current'
        # Merge overlay if present
        if overlay_df is not None and not overlay_df.empty:
            rama_df = pd.concat([rama_df, overlay_df], ignore_index=True)

        rama_df = self._augment_rama_classification(rama_df, analysis_result)
        # Polygon overlay class selection (optional distinction from data filter)
        show_regions = not hide_polygons
        polygon_class_selection = []
        if show_regions:
            overlays_expanded = class_mode != 'All'
            with st.expander("Region Overlays", expanded=overlays_expanded):
                overlay_mode = st.radio(
                    "Overlay Classes",
                    ["Match Data","Custom"],
                    horizontal=True,
                    help="Use the same classes as data points or choose a custom subset for polygon overlays.",
                    key='rama_overlay_mode'
                )
                if overlay_mode == "Custom":
                    polygon_class_selection = st.multiselect(
                        "Polygon classes",
                        ["General","Pre-Pro","Gly","Pro"],
                        default=st.session_state.get('rama_overlay_classes', class_filter),
                        help="Classes whose favored/allowed boundaries will be drawn.",
                        key='rama_overlay_classes'
                    )
                else:
                    polygon_class_selection = class_filter
        # Apply authoritative polygon classification
        rama_df = self._apply_polygon_classification(rama_df)
        # Aggregate interaction counts
        rama_df = self._apply_interaction_counts(rama_df, analysis_result)
        pre_filter_count = len(rama_df)
        if filter_outliers:
            rama_df = rama_df[rama_df['Region'] != 'Outlier']
        if hotspots_only:
            rama_df = rama_df[rama_df['Type'] == 'Interaction Hotspot']
        if class_filter:
            rama_df = rama_df[rama_df['ResidueClass'].isin(class_filter)]

        # Class counts caption & empty state handling
        if pre_filter_count:
            counts = rama_df['ResidueClass'].value_counts().to_dict()
            st.caption(" | ".join([f"{k}: {v}" for k,v in counts.items()]) if counts else "No residues after filters")
        if rama_df is None or rama_df.empty:
            st.warning("No residues remain after current filters. Adjust filters or disable 'Hide Outliers'/'Hotspots Only'.")
            if st.button("Reset Ramachandran Filters"):
                st.experimental_rerun()
            return

        density_traces = []
        effective_show_density = show_density
        # Apply region focus for density and points if selected
        focus_region = st.session_state.get('rama_region_focus', 'All')
        df_for_density = rama_df if focus_region == 'All' else rama_df[rama_df['Region'] == focus_region]
        residue_count = len(df_for_density)
        if residue_count > 5000 and show_density:
            st.info("Density automatically disabled for >5000 residues (performance).")
            effective_show_density = False
        if effective_show_density:
            density_traces = self._kde_density_traces(
                df_for_density,
                high_res=residue_count < 2000 if grid_points is None else (grid_points >= 140),
                bw_override=bw_method,
                npts_override=grid_points,
                quantiles=quantile_contours
            )

        # Scatter with interaction frequency sizing
        display_df = df_for_density
        original_count = residue_count
        target_cap = 8000
        if residue_count > target_cap:
            # Stratified downsample by Region to preserve distribution
            rng = np.random.default_rng(42)
            parts = []
            for region, group in display_df.groupby('Region'):
                proportion = len(group) / residue_count
                take = max(1, int(proportion * target_cap))
                if take >= len(group):
                    parts.append(group)
                else:
                    idx = rng.choice(group.index, size=take, replace=False)
                    parts.append(group.loc[idx])
            display_df = pd.concat(parts, ignore_index=True)
            st.info(f"Downsampled residues for rendering: showing {len(display_df)} of {original_count} (stratified by region) for performance.")

        size_series = display_df['InteractionCount'].astype(float)
        max_val = size_series.max() if len(size_series) else 0
        if max_val <= 0 or np.isnan(max_val):
            max_val = 1.0
        # Build separate traces per region so legend double-click isolates correctly
        region_traces = []
        show_colorbar = interaction_color_scale  # show colorbar on first region trace only
        # Region color scheme (classic: Favored=green, Allowed=gold, Outlier=red)
        REGION_COLORS = {
            'Favored': '#2ca02c',    # green
            'Allowed': '#d4aa00',    # golden yellow
            'Outlier': '#d62728'     # red
        }
        # Prepare symbol mapping if multi-overlay (different structure sources)
        symbol_map = {}
        if multi_overlay and 'SourcePDB' in display_df.columns:
            unique_sources = list(display_df['SourcePDB'].unique())
            # Cap number of distinct symbols for performance/readability
            SYMBOLS = ['circle','square','diamond','cross','x','triangle-up','triangle-down','triangle-left','triangle-right','pentagon']
            for i, src in enumerate(unique_sources):
                symbol_map[src] = SYMBOLS[i % len(SYMBOLS)]
        for region_name, subset in display_df.groupby('Region'):
            sizes = 6 + (np.sqrt(subset['InteractionCount'].astype(float) / max_val)) * 14
            # If multi-overlay: split further by SourcePDB to provide legend entries with symbols
            if multi_overlay and 'SourcePDB' in subset.columns and len(symbol_map) <= 15:
                for source_id, src_group in subset.groupby('SourcePDB'):
                    marker = dict(
                        size=6 + (np.sqrt(src_group['InteractionCount'].astype(float) / max_val)) * 14,
                        line=dict(width=0.5, color='rgba(20,20,20,0.6)'),
                        opacity=0.35 if outlier_focus else 0.9,
                        symbol=symbol_map.get(source_id, 'circle')
                    )
                    if interaction_color_scale:
                        marker.update(dict(
                            color=src_group['InteractionCount'],
                            colorscale='Turbo',
                            colorbar=dict(title='Interactions') if show_colorbar else None
                        ))
                        show_colorbar = False
                    else:
                        marker.update(dict(color=REGION_COLORS.get(region_name, '#4c78a8')))
                    hover_txt = src_group.apply(lambda r: f"{r['Residue']}{r['ResSeq']} {r['Chain']}<br>Class: {r['ResidueClass']}<br>Region: {r['Region']}<br>Interactions: {r['InteractionCount']}<br>Structure: {r.get('SourcePDB','')}" , axis=1)
                    region_traces.append(go.Scattergl(
                        x=src_group['Phi'],
                        y=src_group['Psi'],
                        mode='markers',
                        marker=marker,
                        text=hover_txt,
                        hoverinfo='text',
                        name=f"{region_name} · {source_id}"
                    ))
                continue  # handled splitting
            marker = dict(
                size=sizes,
                line=dict(width=0.5, color='rgba(20,20,20,0.6)'),
                opacity=0.35 if outlier_focus else 0.9
            )
            if interaction_color_scale:
                marker.update(dict(
                    color=subset['InteractionCount'],
                    colorscale='Turbo',
                    colorbar=dict(title='Interactions') if show_colorbar else None
                ))
                # ensure only one colorbar
                show_colorbar = False
            else:
                marker.update(dict(color=REGION_COLORS.get(region_name, '#4c78a8')))
            # If multi overlay and SourcePDB present, encode symbol by SourcePDB to differentiate
            hover_txt = subset.apply(lambda r: f"{r['Residue']}{r['ResSeq']} {r['Chain']}<br>Class: {r['ResidueClass']}<br>Region: {r['Region']}<br>Interactions: {r['InteractionCount']}" + (f"<br>Structure: {r['SourcePDB']}" if 'SourcePDB' in subset.columns else ''), axis=1)
            region_traces.append(go.Scattergl(
                x=subset['Phi'],
                y=subset['Psi'],
                mode='markers',
                marker=marker,
                text=hover_txt,
                hoverinfo='text',
                name=region_name
            ))

        # Outlier halo layer
        halo_traces = []
        if outlier_focus:
            outliers_df = display_df[display_df['Region'] == 'Outlier']
            if not outliers_df.empty:
                halo_traces.append(go.Scattergl(
                    x=outliers_df['Phi'],
                    y=outliers_df['Psi'],
                    mode='markers',
                    marker=dict(
                        size= (6 + (np.sqrt(outliers_df['InteractionCount'].astype(float) / max_val)) * 14) * 1.35,
                        color='rgba(214,39,40,0.55)',
                        line=dict(width=1.5, color='rgba(255,255,255,0.8)')
                    ),
                    hoverinfo='skip',
                    name='Outlier Halo',
                    showlegend=False
                ))

        # Legend below (Plotly native); removed duplicate manual legend per user request

        shapes = self._polygon_region_overlays(class_filter=polygon_class_selection) if show_regions else []
        # Optionally hide density colorbar
        if hide_density_scale:
            for tr in density_traces:
                if hasattr(tr, 'showscale'):
                    tr.showscale = False
        fig = go.Figure(data=[*density_traces, *region_traces, *halo_traces])
        fig.update_layout(
            title='Ramachandran Plot',
            xaxis=dict(range=[-180, 180], dtick=60, title='φ (degrees)'),
            yaxis=dict(range=[-180, 180], dtick=60, title='ψ (degrees)'),
            width=770,
            height=670,
            legend=dict(
                orientation='h',
                yanchor='bottom', y=-0.30,
                xanchor='center', x=0.5,
                title=dict(text=''),
                font=dict(size=11),
                itemclick='toggle',
                itemdoubleclick='toggleothers'
            ),
            shapes=shapes,
            margin=dict(l=60, r=70, t=70, b=210),
            plot_bgcolor='rgba(10,10,30,0.85)',
            paper_bgcolor='rgba(0,0,0,0)'
        )
        fig.update_xaxes(title_standoff=18)
        fig.update_yaxes(title_standoff=18)
        # Updated for Streamlit API deprecation: use width parameter instead of use_container_width
        st.plotly_chart(fig, width="stretch")
        st.caption("Point size ∝ √(interaction count). Use legend to reveal reference markers.")
        # Export
        csv = rama_df.to_csv(index=False).encode()
        st.download_button("💾 Download Ramachandran CSV", csv, file_name="ramachandran_data.csv", mime="text/csv")
        # Region summary export
        summary = rama_df.groupby(['ResidueClass','Region']).size().unstack(fill_value=0)
        summary['Total'] = summary.sum(axis=1)
        pct = (summary.div(summary['Total'], axis=0) * 100).round(2)
        pct.columns = [f"{c}_pct" for c in pct.columns]
        region_summary_df = pd.concat([summary, pct], axis=1).reset_index()
        st.download_button(
            "📊 Download Region Summary CSV",
            region_summary_df.to_csv(index=False).encode(),
            file_name="ramachandran_region_summary.csv",
            mime="text/csv"
        )

        # Show statistics
        if st.checkbox("Show Ramachandran Statistics"):
            self._render_rama_statistics(rama_df)

    def _placeholder_rama_dataframe(self, analysis_result: Dict[str, Any]) -> pd.DataFrame:
        np.random.seed(42)
        n_residues = 250
        phi_angles = np.concatenate([
            np.random.normal(-60, 30, 200),
            np.random.normal(-120, 20, 50)
        ])
        psi_angles = np.concatenate([
            np.random.normal(-45, 30, 200),
            np.random.normal(120, 20, 50)
        ])
        residue_types = ['Normal'] * len(phi_angles)
        hotspots = analysis_result.get('hotspots', [])
        for i in range(min(20, len(hotspots))):
            residue_types[i] = 'Interaction Hotspot'
        df = pd.DataFrame({'Phi': phi_angles, 'Psi': psi_angles, 'Type': residue_types})
        df['Chain'] = 'A'
        df['ResSeq'] = range(1, len(df)+1)
        df['Residue'] = 'ALA'
        df['InteractionCount'] = 0
        return df

    def _compute_phi_psi_dataframe(self, structure, analysis_result: Dict[str, Any]) -> Optional[pd.DataFrame]:
        try:
            from Bio.PDB.Polypeptide import PPBuilder
            ppb = PPBuilder()
            rows = []
            hotspots = analysis_result.get('hotspots', [])
            hotspot_keys = set()
            for h in hotspots:
                chain_id = h.get('chain') or h.get('chain_id') or h.get('chain1')
                resnum = h.get('residue_number') or h.get('resnum') or h.get('residue')
                if chain_id is not None and resnum is not None:
                    hotspot_keys.add((chain_id, resnum))
            for poly in ppb.build_peptides(structure):
                phi_psi = poly.get_phi_psi_list()
                for residue, (phi, psi) in zip(poly, phi_psi):
                    if phi is None or psi is None:
                        continue
                    chain_id = residue.get_parent().id
                    resseq = residue.id[1]
                    rows.append({
                        'Chain': chain_id,
                        'ResSeq': resseq,
                        'Residue': residue.get_resname(),
                        'Phi': np.degrees(phi),
                        'Psi': np.degrees(psi),
                        'Type': 'Interaction Hotspot' if (chain_id, resseq) in hotspot_keys else 'Normal'
                    })
            if not rows:
                return None
            df = pd.DataFrame(rows)
            # Add interaction count if hotspot list includes frequency info (placeholder: mark hotspots as 3 interactions)
            df['InteractionCount'] = df.apply(lambda r: 3 if r['Type'] == 'Interaction Hotspot' else 0, axis=1)
            return df
        except Exception as e:
            return None

    def _apply_interaction_counts(self, df: pd.DataFrame, analysis_result: Dict[str, Any]) -> pd.DataFrame:
        if df is None or df.empty:
            return df
        interactions = analysis_result.get('interactions', {})
        counter = {}
        for itype, lst in interactions.items():
            for inter in lst:
                chain1 = self._get_interaction_property(inter, 'chain1', None)
                res1 = self._get_interaction_property(inter, 'residue1', None)
                chain2 = self._get_interaction_property(inter, 'chain2', None)
                res2 = self._get_interaction_property(inter, 'residue2', None)
                # residue strings like ALA123; extract numeric part
                def parse_res(res):
                    if not res: return None
                    # find trailing digits
                    num = ''.join(ch for ch in res if ch.isdigit())
                    return int(num) if num.isdigit() else None
                if chain1 and res1:
                    rnum = parse_res(res1)
                    if rnum is not None:
                        counter[(chain1, rnum)] = counter.get((chain1, rnum), 0) + 1
                if chain2 and res2:
                    rnum = parse_res(res2)
                    if rnum is not None:
                        counter[(chain2, rnum)] = counter.get((chain2, rnum), 0) + 1
        if counter:
            df['InteractionCount'] = df.apply(lambda r: counter.get((r['Chain'], r['ResSeq']), 0), axis=1)
        return df

    def _kde_density_traces(self, df: pd.DataFrame, high_res: bool = True, bw_override: Optional[float] = None, npts_override: Optional[int] = None, quantiles: Optional[List[float]] = None):
        try:
            from scipy.stats import gaussian_kde
            values = np.vstack([df['Phi'], df['Psi']])
            bw = bw_override if bw_override is not None else (0.25 if high_res else 0.35)
            npts = npts_override if npts_override is not None else (160 if high_res else 90)
            kde = gaussian_kde(values, bw_method=bw)
            grid_x = np.linspace(-180, 180, npts)
            grid_y = np.linspace(-180, 180, npts)
            xx, yy = np.meshgrid(grid_x, grid_y)
            positions = np.vstack([xx.ravel(), yy.ravel()])
            zz = np.reshape(kde(positions), xx.shape)
            heat = go.Contour(
                x=grid_x,
                y=grid_y,
                z=zz,
                colorscale='Viridis',
                contours=dict(showlines=False),
                opacity=0.55,
                showscale=True,
                name='Density',
                showlegend=False
            )
            traces = [heat]
            if quantiles:
                # Convert density surface to cumulative levels, then draw one contour per requested quantile
                z_flat = zz.ravel()
                z_sorted = np.sort(z_flat)[::-1]
                csum = np.cumsum(z_sorted)
                csum /= csum[-1] if csum[-1] != 0 else 1.0
                levels = []
                for q in quantiles:
                    # clamp q to [0,1]
                    q = max(0.0, min(1.0, q))
                    idx = np.searchsorted(csum, q)
                    thr = z_sorted[min(idx, len(z_sorted)-1)]
                    levels.append(thr)
                for thr in levels:
                    traces.append(go.Contour(
                        x=grid_x,
                        y=grid_y,
                        z=zz,
                        autocontour=False,
                        contours=dict(
                            start=thr,
                            end=thr,
                            size=1e-9,  # effectively draw a single level
                            coloring='none',
                            showlines=True
                        ),
                        line=dict(color='white', width=1.2),
                        showscale=False,
                        name='HDR Contour',
                        hoverinfo='skip',
                        showlegend=False
                    ))
            return traces
        except Exception:
            hist = go.Histogram2d(
                x=df['Phi'],
                y=df['Psi'],
                colorscale='Viridis',
                showscale=True,
                nbinsx=60,
                nbinsy=60,
                opacity=0.5,
                name='Density',
                showlegend=False
            )
            return [hist]

    def _augment_rama_classification(self, df: pd.DataFrame, analysis_result: Dict[str, Any]) -> pd.DataFrame:
        if df is None or df.empty:
            return df
        # Residue class determination (simple): GLY, PRO, PRE-PRO (residue preceding a proline), OTHER
        # Build mapping (chain, resseq)->resname for lookahead
        chain_map = {}
        for _, row in df.iterrows():
            chain_map.setdefault(row['Chain'], {})[row['ResSeq']] = row['Residue']
        residue_classes = []
        for _, row in df.iterrows():
            resname = row['Residue']
            chain = row['Chain']
            resseq = row['ResSeq']
            if resname == 'GLY':
                residue_classes.append('Gly')
            elif resname == 'PRO':
                residue_classes.append('Pro')
            else:
                # Pre-Pro check (next residue is PRO)
                next_res = chain_map.get(chain, {}).get(resseq + 1, '')
                if next_res == 'PRO':
                    residue_classes.append('Pre-Pro')
                else:
                    residue_classes.append('General')
        df['ResidueClass'] = residue_classes
        return df

    # Polygon helpers
    def _load_external_polygons(self) -> Optional[Dict[str, Dict[str, List[List[List[float]]]]]]:
        """Load Ramachandran polygons from packaged JSON if available.
        Returns the 'polygons' map or None on failure.
        """
        try:
            base_dir = os.path.dirname(__file__)
            json_path = os.path.join(base_dir, 'data', 'ramachandran_polygons.json')
            if not os.path.exists(json_path):
                return None
            with open(json_path, 'r') as f:
                data = json.load(f)
            polys = data.get('polygons')
            # Basic validation shape
            if not isinstance(polys, dict):
                return None
            return polys
        except Exception:
            return None

    def _default_polygons(self) -> Dict[str, Dict[str, List[List[List[float]]]]]:
        """Fallback minimal polygon set when JSON isn't available."""
        return {
            'General': {
                'favored': [[[-160,-40], [-140,-60], [-110,-70], [-80,-55], [-60,-40], [-50,-20], [-60,10], [-80,30], [-120,20], [-150,-10]]],
                'allowed': [[[-180,-80], [-170,-120], [-140,-140], [-110,-135], [-70,-120], [-40,-90], [-30,-40], [-40,40], [-70,70], [-110,60], [-150,40]]]
            },
            'Pre-Pro': {
                'favored': [[[-160,-50], [-140,-65], [-120,-70], [-95,-60], [-80,-40], [-90,-15], [-120,-5], [-150,-25]]],
                'allowed': [[[-175,-85], [-160,-110], [-130,-125], [-100,-120], [-70,-100], [-55,-65], [-60,-10], [-100,15], [-140,5], [-165,-20]]]
            },
            'Gly': {
                'favored': [[[-80,140], [-60,120], [-40,100], [-20,80], [0,70], [20,80], [30,100], [25,130], [0,150], [-40,155], [-70,150]],
                            [[40,-10], [60,-30], [80,-40], [120,-30], [140,0], [130,30], [100,40], [70,30], [50,10]]],
                'allowed': [[[-100,180], [-40,180], [10,160], [40,140], [50,110], [45,70], [20,40], [-20,30], [-60,40], [-90,70], [-110,110]],
                            [[20,-60], [50,-80], [90,-90], [150,-60], [180,-10], [170,40], [130,60], [90,55], [60,40], [30,10]]]
            },
            'Pro': {
                'favored': [[[-85,140], [-70,130], [-55,115], [-45,95], [-40,70], [-45,40], [-65,20], [-90,25], [-105,50], [-110,85], [-100,120]]],
                'allowed': [[[-100,155], [-60,155], [-35,130], [-30,90], [-35,50], [-55,15], [-90,5], [-115,30], [-120,70]]]
            }
        }

    def _ramachandran_polygons(self) -> Dict[str, Dict[str, List[List[List[float]]]]]:
        polys = self._load_external_polygons()
        if polys is None:
            polys = self._default_polygons()
        return polys

    def _point_in_polygon(self, x: float, y: float, poly: List[List[float]]) -> bool:
        """Ray casting algorithm for point-in-polygon.
        poly: list of [x, y] vertices
        """
        inside = False
        n = len(poly)
        for i in range(n):
            x1, y1 = poly[i]
            x2, y2 = poly[(i + 1) % n]
            if ((y1 > y) != (y2 > y)):
                xinters = (x2 - x1) * (y - y1) / ((y2 - y1) if (y2 - y1) != 0 else 1e-12) + x1
                if x < xinters:
                    inside = not inside
        return inside

    def _apply_polygon_classification(self, df: pd.DataFrame) -> pd.DataFrame:
        if df is None or df.empty:
            return df
        poly_map = self._ramachandran_polygons()
        labels = []
        for _, row in df.iterrows():
            pdata = poly_map.get(row['ResidueClass'], poly_map['General'])
            phi, psi = row['Phi'], row['Psi']
            region = 'Outlier'
            for poly in pdata['favored']:
                if self._point_in_polygon(phi, psi, poly):
                    region = 'Favored'
                    break
            if region == 'Outlier':
                for poly in pdata['allowed']:
                    if self._point_in_polygon(phi, psi, poly):
                        region = 'Allowed'
                        break
            labels.append(region)
        df['Region'] = labels
        cmap = {'Favored':'#2ca02c','Allowed':'#ffbf00','Outlier':'#d62728'}
        df['CategoryColor'] = df['Region'].map(cmap).fillna('#1f77b4')
        return df

    def _polygon_region_overlays(self, class_filter: Optional[List[str]] = None):
        poly_map = self._ramachandran_polygons()
        shapes = []
        for cls, pdata in poly_map.items():
            if class_filter and cls not in class_filter:
                continue
            for poly in pdata['favored']:
                shapes.append(dict(type='path', path='M '+' L '.join(f'{x},{y}' for x,y in poly)+' Z', line=dict(color='rgba(0,128,0,0.35)', width=1), fillcolor='rgba(0,128,0,0.04)'))
            for poly in pdata['allowed']:
                shapes.append(dict(type='path', path='M '+' L '.join(f'{x},{y}' for x,y in poly)+' Z', line=dict(color='rgba(255,191,0,0.3)', width=1), fillcolor='rgba(255,191,0,0.03)'))
        shapes.append(dict(type='rect', x0=-180, y0=-180, x1=180, y1=180, line=dict(color='lightgray', width=1)))
        return shapes

    # Legacy rectangular classification helpers removed in favor of polygon-based approach
    
    def _render_rama_statistics(self, rama_df: pd.DataFrame):
        """Render Ramachandran plot statistics with residue class details."""
        if rama_df is None or rama_df.empty:
            st.info("No Ramachandran data available")
            return
        total = len(rama_df)
        favored = (rama_df['Region'] == 'Favored').sum()
        allowed = (rama_df['Region'] == 'Allowed').sum()
        outliers = (rama_df['Region'] == 'Outlier').sum()
        hotspot_count = (rama_df['Type'] == 'Interaction Hotspot').sum()

        c1, c2, c3, c4 = st.columns(4)
        with c1: st.metric("Total Residues", total)
        with c2: st.metric("Favored", f"{favored/total*100:.1f}%")
        with c3: st.metric("Allowed", f"{allowed/total*100:.1f}%")
        with c4: st.metric("Outliers", f"{outliers/total*100:.1f}%")

        st.caption(f"Hotspot residues (interaction-enriched): {hotspot_count}")

        with st.expander("Residue Class Breakdown"):
            class_summary = rama_df.groupby(['ResidueClass','Region']).size().unstack(fill_value=0)
            class_summary['Total'] = class_summary.sum(axis=1)
            class_summary_pct = (class_summary.div(class_summary['Total'], axis=0) * 100).round(1)
            st.write("Counts")
            st.dataframe(class_summary)
            st.write("Percentages")
            st.dataframe(class_summary_pct)
    
    def render_interaction_network(self, analysis_result: Dict[str, Any]):
        """Render network graph of residue interactions."""
        st.subheader("🕸️ Interaction Network")
        ctrl_col1, ctrl_col2, ctrl_col3 = st.columns(3)
        with ctrl_col1:
            min_edge_count = st.number_input("Min Edge Count", min_value=1, value=1, help="Hide edges with fewer than this number of occurrences across interaction types")
        with ctrl_col2:
            collapse_types = st.checkbox("Collapse Interaction Types", value=False, help="Aggregate all interaction types into weighted edges")
        with ctrl_col3:
            show_labels = st.checkbox("Show Node Labels", value=True)
        # Attempt force-directed layout using networkx (spring layout). Fallback to circular.
        interactions = analysis_result.get('interactions', {})
        total_interactions = sum(len(int_list) for int_list in interactions.values())
        if total_interactions == 0:
            st.info("No interactions available for network graph")
            return
        try:
            import networkx as nx, math
            G = nx.Graph()
            edges = []
            for itype, lst in interactions.items():
                for rec in lst:
                    r1 = getattr(rec, 'residue1', None) if hasattr(rec, 'residue1') else rec.get('residue1') if isinstance(rec, dict) else None
                    r2 = getattr(rec, 'residue2', None) if hasattr(rec, 'residue2') else rec.get('residue2') if isinstance(rec, dict) else None
                    if not r1 or not r2:
                        continue
                    G.add_node(r1)
                    G.add_node(r2)
                    G.add_edge(r1, r2, interaction=itype)
                    edges.append((r1, r2, itype))
            if G.number_of_edges() == 0:
                st.info("Not enough structured interaction data for network rendering")
                return
            try:
                pos = nx.spring_layout(G, seed=42, k=None)
            except Exception:
                # fallback circular
                pos = {node: (math.cos(2*math.pi*i/G.number_of_nodes()), math.sin(2*math.pi*i/G.number_of_nodes())) for i, node in enumerate(G.nodes())}
            # Build traces grouped by interaction type for colored edges
            from collections import defaultdict
            grouped = defaultdict(list)
            for a,b,itype in edges:
                grouped[itype].append((a,b))
            edge_traces = []
            # Aggregate counts per (a,b,itype) since duplicates may exist
            from collections import Counter
            aggregated = Counter()
            for a,b,itype in edges:
                key = tuple(sorted([a,b])) + (itype,)
                aggregated[key] += 1
            type_group_map = defaultdict(list)
            for (a,b,itype), cnt in aggregated.items():
                type_group_map[itype].append((a,b,cnt))
            if collapse_types:
                from collections import defaultdict
                collapsed = defaultdict(int)
                for (a,b,itype), cnt in aggregated.items():
                    if cnt < min_edge_count: continue
                    key = tuple(sorted([a,b]))
                    collapsed[key] += cnt
                xs, ys, hovertexts, widths = [], [], [], []
                for (a,b), cnt in collapsed.items():
                    xs.extend([pos[a][0], pos[b][0], None])
                    ys.extend([pos[a][1], pos[b][1], None])
                    hovertexts.extend([f"{a}–{b}: {cnt} interactions", f"{a}–{b}: {cnt} interactions", None])
                    widths.append(cnt)
                edge_traces.append(go.Scatter(x=xs, y=ys, mode='lines', line=dict(width=2, color='#888'), hoverinfo='text', text=hovertexts, name='All Types'))
            else:
                for itype, items in type_group_map.items():
                    xs, ys, hovertexts = [], [], []
                    for a,b,cnt in items:
                        if cnt < min_edge_count: continue
                        xs.extend([pos[a][0], pos[b][0], None])
                        ys.extend([pos[a][1], pos[b][1], None])
                        hovertexts.append(f"{a}–{b}: {cnt} {itype} interactions")
                        hovertexts.append(hovertexts[-1])
                        hovertexts.append(None)
                    if xs:
                        edge_traces.append(go.Scatter(x=xs, y=ys, mode='lines', line=dict(width=1, color=COLORBLIND_PALETTE.get(normalize_key(itype), '#999999')), hoverinfo='text', text=hovertexts, name=itype))
            # Nodes
            degrees = dict(G.degree())
            node_trace = go.Scatter(
                x=[pos[n][0] for n in G.nodes()],
                y=[pos[n][1] for n in G.nodes()],
                mode='markers'+('+text' if show_labels else ''),
                text=[n for n in G.nodes()],
                textposition='top center',
                marker=dict(size=[6 + degrees[n]*0.8 for n in G.nodes()], color='#ffffff', line=dict(width=1,color='#333')),
                hovertext=[f"{n} (deg {degrees[n]})" for n in G.nodes()],
                hoverinfo='text',
                showlegend=False
            )
            fig = go.Figure(data=edge_traces + [node_trace])
            fig.update_layout(title=f"Interaction Network (nodes={G.number_of_nodes()}, edges={G.number_of_edges()})", showlegend=True, xaxis=dict(visible=False), yaxis=dict(visible=False), height=550)
            st.plotly_chart(fig, use_container_width=True)
        except Exception as e:
            st.warning(f"Network rendering unavailable: {e}")
    
    def render_comparative_analysis(self, batch_results: Dict[str, Dict[str, Any]]):
        """Render comparative analysis across multiple structures."""
        st.subheader("🔄 Comparative Analysis")
        
        if len(batch_results) < 2:
            st.warning("Need at least 2 structures for comparison")
            return
        
        # Create comparison DataFrame
        comparison_data = []
        
        for pdb_id, result in batch_results.items():
            if not result:
                continue
            
            interactions = result.get('interactions', {})
            
            row = {'PDB ID': pdb_id}
            
            for interaction_type, interaction_list in interactions.items():
                display_name = self._get_display_name(interaction_type)
                row[display_name] = len(interaction_list)
            
            comparison_data.append(row)
        
        if not comparison_data:
            st.warning("No data available for comparison")
            return
        
        df = pd.DataFrame(comparison_data).fillna(0)
        
        # Create comparison plot
        comparison_type = st.selectbox(
            "Comparison Type:",
            ["Stacked Bar", "Heatmap", "Correlation Matrix"]
        )
        
        if comparison_type == "Stacked Bar":
            self._render_stacked_comparison(df)
        
        elif comparison_type == "Heatmap":
            self._render_comparison_heatmap(df)
        
        elif comparison_type == "Correlation Matrix":
            self._render_correlation_matrix(df)
    
    def _render_stacked_comparison(self, df: pd.DataFrame):
        """Render stacked bar chart for comparison."""
        # Prepare data for stacked bar
        df_melted = df.melt(
            id_vars=['PDB ID'],
            var_name='Interaction Type',
            value_name='Count'
        )
        
        fig = px.bar(
            df_melted,
            x='PDB ID',
            y='Count',
            color='Interaction Type',
            title="Interaction Count Comparison Across Structures"
        )
        
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, width="stretch")
    
    def _render_comparison_heatmap(self, df: pd.DataFrame):
        """Render heatmap for structure comparison."""
        # Set PDB ID as index
        df_indexed = df.set_index('PDB ID')
        
        fig = px.imshow(
            df_indexed.T,
            labels=dict(x="PDB ID", y="Interaction Type", color="Count"),
            title="Interaction Pattern Comparison"
        )
        
        st.plotly_chart(fig, width="stretch")
    
    def _render_correlation_matrix(self, df: pd.DataFrame):
        """Render correlation matrix of interaction types."""
        # Calculate correlation matrix (excluding PDB ID)
        numeric_df = df.select_dtypes(include=[np.number])
        corr_matrix = numeric_df.corr()
        
        fig = px.imshow(
            corr_matrix,
            labels=dict(color="Correlation"),
            title="Interaction Type Correlation Matrix"
        )
        
        st.plotly_chart(fig, width="stretch")
    
    def _get_display_name(self, interaction_type: str) -> str:
        """Get human-readable name for interaction type."""
        display_names = {
            'hydrogen_bond': 'Hydrogen Bonds',
            'halogen_bond': 'Halogen Bonds',
            'pi_pi': 'π-π Stacking',
            'ionic': 'Ionic Interactions',
            'hydrophobic': 'Hydrophobic Contacts',
            'ch_pi': 'C-H···π Interactions',
            'chalcogen_bond': 'Chalcogen Bonds',
            'pnictogen_bond': 'Pnictogen Bonds',
            'tetrel_bond': 'Tetrel Bonds',
            'anion_pi': 'Anion-π Interactions',
            'n_pi_star': 'n→π* Interactions',
            'dispersion': 'London Dispersion'
        }
        return display_names.get(interaction_type, interaction_type.replace('_', ' ').title())
    
    def render_summary_dashboard(self, analysis_result: Dict[str, Any]):
        """Render a summary dashboard with key metrics."""
        st.subheader("📈 Analysis Summary Dashboard")
        
        # Key metrics
        interactions = analysis_result.get('interactions', {})
        total_interactions = sum(len(int_list) for int_list in interactions.values())
        hotspots = analysis_result.get('hotspots', [])
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Interactions", total_interactions)
        
        with col2:
            st.metric("Interaction Types", len([t for t, l in interactions.items() if l]))
        
        with col3:
            st.metric("Hotspot Residues", len(hotspots))
        
        with col4:
            processing_time = analysis_result.get('metadata', {}).get('analysis_time', 0)
            st.metric("Analysis Time", f"{processing_time:.2f}s")
        
        # Quick visualization
        if total_interactions > 0:
            interaction_counts = {
                self._get_display_name(int_type): len(int_list)
                for int_type, int_list in interactions.items()
                if int_list
            }
            
            # Mini bar chart
            fig = go.Figure(data=[
                go.Bar(
                    x=list(interaction_counts.values()),
                    y=list(interaction_counts.keys()),
                    orientation='h',
                    marker_color=self.color_palette[:len(interaction_counts)]
                )
            ])
            
            fig.update_layout(
                title="Interaction Count Summary",
                height=300,
                margin=dict(l=0, r=0, t=30, b=0)
            )
            
            st.plotly_chart(fig, width="stretch")
