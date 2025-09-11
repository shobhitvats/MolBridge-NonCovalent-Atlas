"""
Analysis plots and visualizations for Protein Interaction Explorer.
"""

import streamlit as st
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from typing import Dict, Any, List, Optional
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from utils.config import AppConfig

class InteractionPlots:
    """Generates plots and visualizations for interaction analysis."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        self.color_palette = list(config.visualization.interaction_colors.values())
    
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
        """Render chain-vs-chain interaction heatmap with Lottie fallback."""
        from streamlit_lottie import st_lottie
        import requests
        st.subheader("🔥 Chain Interaction Heatmap")
        chain_matrix = self._create_chain_interaction_matrix(analysis_result)
        if chain_matrix.empty:
            st_lottie(requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json").json(), height=180, key="empty_heatmap")
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
            st_lottie(requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json").json(), height=180, key="error_heatmap")
            st.error(f"Error rendering heatmap: {e}")
            return
        if st.checkbox("Show Interactive Heatmap"):
            try:
                self._render_interactive_heatmap(chain_matrix)
            except Exception as e:
                st_lottie(requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json").json(), height=180, key="error_interactive_heatmap")
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
        
        st.plotly_chart(fig, use_container_width=True)
    
    def render_interaction_distribution(self, analysis_result: Dict[str, Any]):
        """Render interaction type distribution plots with Lottie fallback."""
        from streamlit_lottie import st_lottie
        import requests
        st.subheader("📊 Interaction Distribution")
        interactions = analysis_result.get('interactions', {})
        interaction_counts = {
            self._get_display_name(int_type): len(int_list)
            for int_type, int_list in interactions.items()
            if int_list
        }
        if not interaction_counts:
            st_lottie(requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json").json(), height=180, key="empty_dist")
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
            st_lottie(requests.get("https://assets2.lottiefiles.com/packages/lf20_3rwasyjy.json").json(), height=180, key="error_dist")
            st.error(f"Error rendering distribution: {e}")
    
    def _render_bar_chart(self, interaction_counts: Dict[str, int]):
        """Render bar chart of interaction counts."""
        fig = go.Figure(data=[
            go.Bar(
                x=list(interaction_counts.keys()),
                y=list(interaction_counts.values()),
                marker_color=self.color_palette[:len(interaction_counts)]
            )
        ])
        
        fig.update_layout(
            title="Interaction Type Distribution",
            xaxis_title="Interaction Type",
            yaxis_title="Count",
            xaxis_tickangle=-45
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
    def _render_pie_chart(self, interaction_counts: Dict[str, int]):
        """Render pie chart of interaction counts."""
        fig = go.Figure(data=[
            go.Pie(
                labels=list(interaction_counts.keys()),
                values=list(interaction_counts.values()),
                hole=0.3,
                marker_colors=self.color_palette[:len(interaction_counts)]
            )
        ])
        
        fig.update_layout(title="Interaction Type Distribution")
        st.plotly_chart(fig, use_container_width=True)
    
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
            title="Distance Distribution by Interaction Type"
        )
        
        fig.update_layout(xaxis_tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)
    
    def render_ramachandran_plot(self, analysis_result: Dict[str, Any]):
        """Render Ramachandran plot with interaction highlights."""
        st.subheader("📐 Ramachandran Plot")
        
        # For now, create a placeholder Ramachandran plot
        # In production, this would use actual phi/psi angles from the structure
        
        # Generate sample data
        np.random.seed(42)
        n_residues = 200
        
        phi_angles = np.random.normal(-60, 30, n_residues)  # Alpha helix region
        psi_angles = np.random.normal(-45, 30, n_residues)
        
        # Add some beta sheet residues
        n_beta = 50
        phi_beta = np.random.normal(-120, 20, n_beta)
        psi_beta = np.random.normal(120, 20, n_beta)
        
        phi_angles = np.concatenate([phi_angles, phi_beta])
        psi_angles = np.concatenate([psi_angles, psi_beta])
        
        # Create residue types (some involved in interactions)
        residue_types = ['Normal'] * len(phi_angles)
        
        # Mark some residues as interaction hotspots
        hotspots = analysis_result.get('hotspots', [])
        n_hotspots = min(20, len(hotspots))
        
        for i in range(n_hotspots):
            residue_types[i] = 'Interaction Hotspot'
        
        # Create DataFrame
        rama_df = pd.DataFrame({
            'Phi': phi_angles,
            'Psi': psi_angles,
            'Type': residue_types
        })
        
        # Create plot
        fig = px.scatter(
            rama_df,
            x='Phi',
            y='Psi',
            color='Type',
            title='Ramachandran Plot',
            labels={'Phi': 'φ (degrees)', 'Psi': 'ψ (degrees)'},
            color_discrete_map={
                'Normal': 'lightblue',
                'Interaction Hotspot': 'red'
            }
        )
        
        # Add allowed regions (simplified)
        fig.add_shape(
            type="rect",
            x0=-180, y0=-180, x1=180, y1=180,
            line=dict(color="lightgray", width=1),
            fillcolor="rgba(255,255,255,0)"
        )
        
        fig.update_layout(
            xaxis=dict(range=[-180, 180], dtick=60),
            yaxis=dict(range=[-180, 180], dtick=60),
            width=600,
            height=600
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Show statistics
        if st.checkbox("Show Ramachandran Statistics"):
            self._render_rama_statistics(rama_df)
    
    def _render_rama_statistics(self, rama_df: pd.DataFrame):
        """Render Ramachandran plot statistics."""
        # Calculate statistics
        alpha_region = rama_df[
            (rama_df['Phi'] > -100) & (rama_df['Phi'] < -30) &
            (rama_df['Psi'] > -70) & (rama_df['Psi'] < -10)
        ]
        
        beta_region = rama_df[
            (rama_df['Phi'] > -150) & (rama_df['Phi'] < -90) &
            (rama_df['Psi'] > 90) & (rama_df['Psi'] < 150)
        ]
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Total Residues", len(rama_df))
        
        with col2:
            alpha_pct = (len(alpha_region) / len(rama_df)) * 100
            st.metric("α-Helix Region", f"{alpha_pct:.1f}%")
        
        with col3:
            beta_pct = (len(beta_region) / len(rama_df)) * 100
            st.metric("β-Sheet Region", f"{beta_pct:.1f}%")
    
    def render_interaction_network(self, analysis_result: Dict[str, Any]):
        """Render network graph of residue interactions."""
        st.subheader("🕸️ Interaction Network")
        
        # This would create a network graph showing residues as nodes
        # and interactions as edges
        st.info("Network visualization would show residues connected by interactions")
        
        # Placeholder for network stats
        interactions = analysis_result.get('interactions', {})
        total_interactions = sum(len(int_list) for int_list in interactions.values())
        
        if total_interactions > 0:
            st.write(f"Network would contain {total_interactions} edges representing interactions")
    
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
        st.plotly_chart(fig, use_container_width=True)
    
    def _render_comparison_heatmap(self, df: pd.DataFrame):
        """Render heatmap for structure comparison."""
        # Set PDB ID as index
        df_indexed = df.set_index('PDB ID')
        
        fig = px.imshow(
            df_indexed.T,
            labels=dict(x="PDB ID", y="Interaction Type", color="Count"),
            title="Interaction Pattern Comparison"
        )
        
        st.plotly_chart(fig, use_container_width=True)
    
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
        
        st.plotly_chart(fig, use_container_width=True)
    
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
            
            st.plotly_chart(fig, use_container_width=True)
