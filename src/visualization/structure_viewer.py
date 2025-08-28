"""
3D structure visualization using py3Dmol for Protein Interaction Explorer.
"""

import streamlit as st
import py3Dmol
import numpy as np
from typing import Dict, Any, List, Optional
from io import StringIO
import json
from loguru import logger

from utils.config import AppConfig
from utils.pdb_handler import PDBHandler

class StructureViewer:
    """Handles 3D visualization of protein structures and interactions."""
    
    def __init__(self, config: AppConfig):
        self.config = config
        self.viewer_config = config.visualization
        self.pdb_handler = PDBHandler(config)
    
    def render_structure(self, 
                        pdb_id: str,
                        analysis_result: Dict[str, Any],
                        selected_interactions: List[str]) -> None:
        """
        Render 3D structure with interactions in Streamlit.
        
        Args:
            pdb_id: PDB identifier
            analysis_result: Results from interaction analysis
            selected_interactions: List of interaction types to display
        """
    # Viewer header removed as per UI cleanup request
        
        # Viewer controls
        col1, col2, col3 = st.columns(3)
        
        with col1:
            style = st.selectbox(
                "Protein Style:",
                ["cartoon", "stick", "sphere", "line"],
                index=0
            )
        
        with col2:
            color_scheme = st.selectbox(
                "Color Scheme:",
                ["spectrum", "chain", "residue", "element"],
                index=0
            )
        
        with col3:
            show_ligands = st.checkbox("Show Ligands", value=True)
        
        # Interaction display controls
        st.write("**Interaction Display:**")
        
        interaction_cols = st.columns(4)
        interaction_visibility = {}
        
        for i, interaction_type in enumerate(selected_interactions):
            col_idx = i % 4
            with interaction_cols[col_idx]:
                display_name = self._get_interaction_display_name(interaction_type)
                color = self.viewer_config.interaction_colors.get(interaction_type, "#888888")
                
                # Create colored checkbox label
                st.markdown(f'<span style="color: {color}">●</span> {display_name}', 
                           unsafe_allow_html=True)
                
                interaction_visibility[interaction_type] = st.checkbox(
                    f"Show {display_name}",
                    value=True,
                    key=f"show_{interaction_type}",
                    label_visibility="collapsed"
                )
        
        # Display viewer
        st.components.v1.html(
            self._create_main_viewer(pdb_id, analysis_result, style, color_scheme, show_ligands, interaction_visibility),
            height=self.viewer_config.viewer_height + 50
        )

    def _create_main_viewer(self, pdb_id: str, analysis_result: Dict[str, Any], style: str, color_scheme: str, show_ligands: bool, interaction_visibility: Dict[str, bool]) -> str:
        """Create the main viewer HTML."""
        # Debug: print analysis result keys
        print(f"DEBUG: Analysis result keys: {list(analysis_result.keys())}")
        
        # Get PDB data from analysis result
        pdb_data = analysis_result.get('pdb_data', '')
        print(f"DEBUG: PDB data from analysis_result: {len(pdb_data) if pdb_data else 0} characters")
        
        if not pdb_data or len(pdb_data) < 100:
            # Fallback: try to get PDB data using the PDB handler
            print(f"DEBUG: Fetching PDB data for {pdb_id} using PDB handler")
            pdb_data = self.pdb_handler.get_pdb_content(pdb_id)
            print(f"DEBUG: PDB data from handler: {len(pdb_data) if pdb_data else 0} characters")
        
        if not pdb_data:
            return f'<div style="color: red; padding: 20px;">❌ Failed to load PDB data for {pdb_id}</div>'
        
        if len(pdb_data) < 100:
            return f'<div style="color: red; padding: 20px;">❌ PDB data too short for {pdb_id} ({len(pdb_data)} chars)</div>'
        
        # Handle PDB data more efficiently - use base64 encoding for large files
        print(f"DEBUG: Original PDB data length: {len(pdb_data)}")
        
        # For very large PDB files, use base64 encoding to avoid JavaScript string issues
        import base64
        pdb_data_b64 = base64.b64encode(pdb_data.encode('utf-8')).decode('ascii')
        print(f"DEBUG: Base64 encoded PDB data length: {len(pdb_data_b64)}")
        
        # Get interactions for visualization
        interactions = analysis_result.get('interactions', {})
        print(f"DEBUG: Available interaction types: {list(interactions.keys())}")
        
        # Debug: print detailed interaction structure
        for interaction_type, interaction_list in interactions.items():
            if interaction_list:
                print(f"DEBUG: Sample {interaction_type} structure: {type(interaction_list[0])}")
                for key, value in interaction_list[0].items():
                    print(f"  {key}: {type(value)} = {value}")
                break
        
        # Create interaction data for JavaScript
        interaction_data = []
        for interaction_type, interaction_list in interactions.items():
            if interaction_visibility.get(interaction_type, True) and interaction_list:
                for interaction in interaction_list:
                    # Convert all numpy types to Python native types for JSON serialization
                    interaction_item = {}
                    
                    # Safely convert all values to JSON-serializable types
                    for key, value in interaction.items():
                        if hasattr(value, 'item'):  # numpy scalar
                            interaction_item[key] = value.item()
                        elif isinstance(value, (int, str, bool, type(None))):
                            interaction_item[key] = value
                        elif isinstance(value, float):
                            interaction_item[key] = float(value)
                        else:
                            interaction_item[key] = str(value)
                    
                    # Ensure required fields exist with safe defaults
                    safe_interaction = {
                        'type': str(interaction_item.get('type', interaction_type)),
                        'residue1': str(interaction_item.get('residue1', '')),
                        'chain1': str(interaction_item.get('chain1', '')),
                        'residue2': str(interaction_item.get('residue2', '')),
                        'chain2': str(interaction_item.get('chain2', '')),
                        'distance': float(interaction_item.get('distance', 0)),
                        'strength': str(interaction_item.get('strength', 'unknown')),
                        'atoms': []
                    }
                    
                    # For now, we'll create placeholder coordinates based on residue positions
                    # In a full implementation, we'd extract actual atomic coordinates from the PDB
                    # This is a simplified approach for basic visualization
                    if safe_interaction['residue1'] and safe_interaction['residue2']:
                        # Extract residue numbers for basic positioning
                        try:
                            res1_num = int(''.join(filter(str.isdigit, safe_interaction['residue1'])))
                            res2_num = int(''.join(filter(str.isdigit, safe_interaction['residue2'])))
                            
                            # Create approximate positions (this would normally come from actual PDB parsing)
                            safe_interaction['atoms'] = [
                                {'x': float(res1_num * 3.8), 'y': 0.0, 'z': 0.0},  # Approximate CA position
                                {'x': float(res2_num * 3.8), 'y': 5.0, 'z': 2.0}   # Approximate CA position
                            ]
                        except:
                            # If residue number parsing fails, use default positions
                            safe_interaction['atoms'] = [
                                {'x': 0.0, 'y': 0.0, 'z': 0.0},
                                {'x': 10.0, 'y': 0.0, 'z': 0.0}
                            ]
                    
                    interaction_data.append(safe_interaction)
        
        import json
        import numpy as np
        
        # Custom JSON encoder to handle numpy types
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                return super().default(obj)
        
        interaction_data_json = json.dumps(interaction_data, cls=NumpyEncoder)
        print(f"DEBUG: Interaction data prepared: {len(interaction_data)} interactions")
        
        # Simple, robust viewer HTML
        main_viewer = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                #main-viewer {{
                    width: 100%;
                    height: {self.viewer_config.viewer_height}px;
                    background: {self.viewer_config.background_color};
                    border: 1px solid #ddd;
                    border-radius: 5px;
                    position: relative;
                }}
                .status {{
                    position: absolute;
                    top: 10px;
                    left: 10px;
                    padding: 8px;
                    border-radius: 4px;
                    z-index: 100;
                    font-size: 14px;
                    font-family: Arial, sans-serif;
                }}
                .loading {{ background: rgba(255,255,0,0.9); color: #000; }}
                .success {{ background: rgba(0,255,0,0.9); color: #000; }}
                .error {{ background: rgba(255,0,0,0.9); color: #fff; }}
            </style>
        </head>
        <body>
            <div id="main-viewer">
                <div id="status" class="status loading">Loading structure...</div>
            </div>
            
            <script>
            (function() {{
                let viewer = null;
                const statusEl = document.getElementById('status');
                const interactionData = {interaction_data_json};
                
                // Interaction type colors
                const interactionColors = {{
                    'hydrogen_bonds': '#FF6B6B',
                    'ionic_interactions': '#4ECDC4', 
                    'hydrophobic_contacts': '#45B7D1',
                    'pi_pi_stacking': '#96CEB4',
                    'chalcogen_bonds': '#FECA57',
                    'halogen_bonds': '#FF9FF3',
                    'ch_pi_interactions': '#54A0FF',
                    'anion_pi_interactions': '#5F27CD',
                    'n_pi_star_interactions': '#00D2D3',
                    'pnictogen_bonds': '#FF6348',
                    'tetrel_bonds': '#2ED573',
                    'london_dispersion': '#A4B0BE'
                }};
                
                function setStatus(msg, type) {{
                    console.log('Status:', msg);
                    statusEl.textContent = msg;
                    statusEl.className = 'status ' + type;
                }}
                
                function addInteractionVisualization() {{
                    console.log('🔗 Adding interaction visualizations...');
                    console.log('Interaction data:', interactionData);
                    console.log('Interaction data length:', interactionData ? interactionData.length : 'undefined');
                    
                    if (!interactionData || interactionData.length === 0) {{
                        console.warn('❌ No interaction data available for visualization');
                        return 0;
                    }}
                    
                    let addedCount = 0;
                    
                    interactionData.forEach((interaction, index) => {{
                        try {{
                            const color = interactionColors[interaction.type] || '#FECA57'; // Default yellow
                            console.log(`Processing interaction ${{index + 1}}: ${{interaction.type}}`);
                            
                            // Get actual atomic coordinates from the 3D model
                            const model = viewer.getModel(0);
                            if (!model) {{
                                console.warn('No model available for coordinate lookup');
                                return;
                            }}
                            
                            // Parse residue information
                            const res1 = parseInt(interaction.residue1.replace(/[^0-9]/g, ''));
                            const res2 = parseInt(interaction.residue2.replace(/[^0-9]/g, ''));
                            const chain1 = interaction.chain1;
                            const chain2 = interaction.chain2;
                            
                            console.log(`Looking for atoms: Chain ${{chain1}} Res ${{res1}} <-> Chain ${{chain2}} Res ${{res2}}`);
                            
                            // Get atoms from the actual residues
                            const atoms1 = model.selectedAtoms({{chain: chain1, resi: res1}});
                            const atoms2 = model.selectedAtoms({{chain: chain2, resi: res2}});
                            
                            if (atoms1.length > 0 && atoms2.length > 0) {{
                                // Find representative atoms (CA for backbone, or first available)
                                let atom1 = atoms1.find(a => a.atom === 'CA') || atoms1[0];
                                let atom2 = atoms2.find(a => a.atom === 'CA') || atoms2[0];
                                
                                console.log(`Found atoms: ${{atom1.atom}} at (${{atom1.x}}, ${{atom1.y}}, ${{atom1.z}}) and ${{atom2.atom}} at (${{atom2.x}}, ${{atom2.y}}, ${{atom2.z}})`);
                                
                                // Add dashed line between actual atomic positions
                                viewer.addCylinder({{
                                    start: {{x: atom1.x, y: atom1.y, z: atom1.z}},
                                    end: {{x: atom2.x, y: atom2.y, z: atom2.z}},
                                    radius: 0.15,
                                    color: color,
                                    dashed: true,
                                    opacity: 0.8
                                }});
                                
                                // Calculate actual midpoint between real atoms
                                const midpoint = {{
                                    x: (atom1.x + atom2.x) / 2,
                                    y: (atom1.y + atom2.y) / 2,
                                    z: (atom1.z + atom2.z) / 2
                                }};
                                
                                // Create detailed label with donor/acceptor information
                                const interactionName = interaction.type.replace('_', ' ').replace('bond', '');
                                const donorInfo = `${{interaction.residue1}} (${{chain1}})`;
                                const acceptorInfo = `${{interaction.residue2}} (${{chain2}})`;
                                const labelText = `${{interactionName}}\\n${{donorInfo}} ↔ ${{acceptorInfo}}\\n${{interaction.distance.toFixed(2)}}Å`;
                                
                                console.log(`Adding detailed label "${{labelText}}" at real midpoint (${{midpoint.x.toFixed(2)}}, ${{midpoint.y.toFixed(2)}}, ${{midpoint.z.toFixed(2)}})`);
                                
                                // Add label at the real midpoint
                                try {{
                                    viewer.addLabel(labelText, {{
                                        position: midpoint,
                                        fontSize: 11,
                                        fontColor: 'white',
                                        backgroundColor: 'black',
                                        backgroundOpacity: 0.95,
                                        showBackground: true,
                                        alignment: 'center',
                                        borderThickness: 1,
                                        borderColor: 'white'
                                    }});
                                    console.log(`✅ Successfully added label for interaction ${{index + 1}}`);
                                }} catch (labelError) {{
                                    console.error(`❌ Failed to add label for interaction ${{index + 1}}:`, labelError);
                                }}
                                
                                // Highlight the specific interacting residues
                                viewer.addStyle({{chain: chain1, resi: res1}}, {{
                                    stick: {{
                                        color: color,
                                        radius: 0.3,
                                        opacity: 0.8
                                    }}
                                }});
                                
                                viewer.addStyle({{chain: chain2, resi: res2}}, {{
                                    stick: {{
                                        color: color,
                                        radius: 0.3,
                                        opacity: 0.8
                                    }}
                                }});
                                
                                addedCount++;
                            }} else {{
                                console.warn(`Could not find atoms for interaction: Chain ${{chain1}} Res ${{res1}} (${{atoms1.length}} atoms) <-> Chain ${{chain2}} Res ${{res2}} (${{atoms2.length}} atoms)`);
                                
                                // Fallback: still highlight the residues even if we can't draw lines
                                if (atoms1.length > 0) {{
                                    viewer.addStyle({{chain: chain1, resi: res1}}, {{
                                        stick: {{color: color, radius: 0.4}}
                                    }});
                                }}
                                if (atoms2.length > 0) {{
                                    viewer.addStyle({{chain: chain2, resi: res2}}, {{
                                        stick: {{color: color, radius: 0.4}}
                                    }});
                                }}
                                addedCount++;
                            }}
                        }} catch (e) {{
                            console.error('Failed to add interaction visualization:', e, interaction);
                        }}
                    }});
                    
                    console.log(`🎯 Summary: Added ${{addedCount}} interaction visualizations out of ${{interactionData.length}} interactions`);
                    
                    if (addedCount === 0) {{
                        console.warn('⚠️ No interactions were successfully visualized. Check console for errors.');
                    }} else {{
                        console.log(`✅ Successfully visualized ${{addedCount}} interactions with labels and bonds`);
                    }}
                    
                    return addedCount;
                }}
                
                // Remove complex functions that might be causing issues
                // Keep only the essential interaction colors and simple functions
                
                function initViewer() {{
                    try {{
                        setStatus('Checking 3Dmol...', 'loading');
                        
                        if (typeof $3Dmol === 'undefined') {{
                            throw new Error('3Dmol.js library not loaded');
                        }}
                        
                        setStatus('Creating viewer...', 'loading');
                        viewer = $3Dmol.createViewer('main-viewer', {{
                            defaultcolors: $3Dmol.rasmolElementColors
                        }});
                        
                        setStatus('Decoding PDB data ({len(pdb_data_b64)} chars)...', 'loading');
                        const pdbData = atob('{pdb_data_b64}');
                        console.log('Decoded PDB data length:', pdbData.length);
                        
                        if (pdbData.length < 100) {{
                            throw new Error('Invalid PDB data - too short');
                        }}
                        
                        setStatus('Loading molecular structure...', 'loading');
                        viewer.addModel(pdbData, 'pdb');
                        
                        setStatus('Applying styling...', 'loading');
                        // Simple but effective styling
                        viewer.setStyle({{}}, {{}}); // Clear all
                        viewer.setStyle({{}}, {{
                            cartoon: {{
                                color: 'spectrum'
                            }}
                        }});
                        
                        setStatus('Adding interactions...', 'loading');
                        const interactionCount = addInteractionVisualization();
                        
                        setStatus('Rendering...', 'loading');
                        viewer.zoomTo();
                        viewer.render();
                        
                        setStatus(`✅ Structure loaded! (${{interactionCount}} interactions)`, 'success');
                        setTimeout(() => statusEl.style.display = 'none', 4000);
                        
                        console.log('✅ Viewer initialization complete!');
                        
                    }} catch (error) {{
                        console.error('❌ Viewer error:', error);
                        setStatus('❌ Error: ' + error.message, 'error');
                    }}
                }}
                
                // Initialize when ready
                document.addEventListener('DOMContentLoaded', function() {{
                    console.log('DOM ready, checking for 3Dmol...');
                    if (typeof $3Dmol !== 'undefined') {{
                        console.log('3Dmol available immediately');
                        initViewer();
                    }} else {{
                        console.log('Waiting for 3Dmol to load...');
                        let attempts = 0;
                        const check = setInterval(() => {{
                            attempts++;
                            if (typeof $3Dmol !== 'undefined') {{
                                console.log('3Dmol loaded after', attempts * 100, 'ms');
                                clearInterval(check);
                                initViewer();
                            }} else if (attempts > 100) {{
                                clearInterval(check);
                                setStatus('❌ 3Dmol.js failed to load after 10 seconds', 'error');
                            }}
                        }}, 100);
                    }}
                }});
            }})();
            </script>
        </body>
        </html>
        """
        
        return main_viewer

    def _render_viewer(self, analysis_results: Dict[str, Any]) -> None:
        """Render the main 3D protein structure viewer."""
        try:
            # Extract data
            pdb_id = analysis_results['pdb_id']
            pdb_data = analysis_results['pdb_data']
            interaction_lines = analysis_results.get('interaction_lines', [])
            
            # Create viewer HTML
            main_viewer = self._create_main_viewer(pdb_id, pdb_data, interaction_lines)
            
            # Render the viewer
            st.components.v1.html(main_viewer, height=self.viewer_config.viewer_height + 50)
            
        except Exception as e:
            st.error(f"Error creating main viewer: {e}")
        
        # Export options
        st.subheader("📥 Export Options")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("📷 PNG Snapshot"):
                st.info("Snapshot functionality would capture current view")
        
        with col2:
            if st.button("🔬 PyMOL Session"):
                pymol_session = self._generate_pymol_session(pdb_id, analysis_results, interaction_lines)
                st.download_button(
                    "Download .pse",
                    pymol_session,
                    f"{pdb_id}_interactions.pse",
                    "application/octet-stream"
                )
        
        with col3:
            if st.button("🧪 ChimeraX Session"):
                chimera_script = self._generate_chimera_script(pdb_id, analysis_results, interaction_lines)
                st.download_button(
                    "Download .cxc",
                    chimera_script,
                    f"{pdb_id}_interactions.cxc",
                    "text/plain"
                )
        
        # Residue search and selection
        self._render_residue_search(analysis_results)
    
    def _create_viewer_html(self,
                           pdb_id: str,
                           analysis_result: Dict[str, Any],
                           style: str,
                           color_scheme: str,
                           show_ligands: bool,
                           interaction_visibility: Dict[str, bool]) -> str:
        """Create HTML for py3Dmol viewer."""
        
        # Get PDB data
        pdb_data = self._get_pdb_data(pdb_id)
        
        # Create interaction lines data
        interaction_lines = self._create_interaction_lines(analysis_result, interaction_visibility)
        
        # Simple JavaScript escaping (like the working simple viewer)
        pdb_data_escaped = pdb_data.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n')
        
        html_template = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
            <style>
                .mol-container {{
                    width: 100%;
                    height: {self.viewer_config.viewer_height}px;
                    position: relative;
                    background-color: {self.viewer_config.background_color};
                    border: 1px solid #ddd;
                    border-radius: 5px;
                }}
                .viewer-controls {{
                    position: absolute;
                    top: 10px;
                    right: 10px;
                    z-index: 100;
                    background: rgba(255,255,255,0.9);
                    padding: 10px;
                    border-radius: 5px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                }}
                .viewer-controls button {{
                    margin: 2px;
                    padding: 5px 10px;
                    border: none;
                    border-radius: 3px;
                    background-color: #007bff;
                    color: white;
                    cursor: pointer;
                }}
                .viewer-controls button:hover {{
                    background-color: #0056b3;
                }}
                .loading-message {{
                    position: absolute;
                    top: 50%;
                    left: 50%;
                    transform: translate(-50%, -50%);
                    z-index: 50;
                    background: rgba(0,0,0,0.7);
                    color: white;
                    padding: 20px;
                    border-radius: 5px;
                }}
                .error-message {{
                    position: absolute;
                    top: 50%;
                    left: 50%;
                    transform: translate(-50%, -50%);
                    z-index: 50;
                    background: rgba(255,0,0,0.8);
                    color: white;
                    padding: 20px;
                    border-radius: 5px;
                    max-width: 80%;
                }}
            </style>
        </head>
        <body>
            <div id="viewer" class="mol-container">
                <div id="loading" class="loading-message">Loading protein structure...</div>
            </div>
            <div class="viewer-controls">
                <button onclick="resetView()">Reset View</button>
                <button onclick="toggleSpin()">Toggle Spin</button>
                <button onclick="centerStructure()">Center</button>
                <button onclick="toggleFullscreen()">Fullscreen</button>
            </div>
            
            <script>
                console.log("Main viewer initializing for {pdb_id}...");
                
                let viewer = null;
                let spinning = false;
                
                // Hide loading message
                function hideLoading() {{
                    const loading = document.getElementById('loading');
                    if (loading) loading.style.display = 'none';
                }}
                
                // Show error message
                function showError(message) {{
                    hideLoading();
                    const errorDiv = document.createElement('div');
                    errorDiv.className = 'error-message';
                    errorDiv.innerHTML = '<strong>Error loading 3D structure:</strong><br>' + message + '<br><br>Check browser console for details.';
                    document.getElementById('viewer').appendChild(errorDiv);
                }}
                
                // Set timeout for initialization
                initializationTimeout = setTimeout(function() {{
                    console.error("TIMEOUT: 3Dmol initialization took too long");
                    showError("Initialization timeout - 3Dmol.js may not have loaded properly");
                }}, 10000); // 10 second timeout
                
                
                // Function to initialize when 3Dmol is ready
                function initialize3Dmol() {{
                    try {{
                        console.log("STEP 1: Starting 3Dmol initialization...");
                        
                        // Check if 3Dmol is loaded
                        if (typeof $3Dmol === 'undefined') {{
                            throw new Error("3Dmol.js library not loaded");
                        }}
                        console.log("STEP 2: 3Dmol library confirmed loaded");
                        
                        console.log("STEP 3: Creating 3Dmol viewer...");
                        viewer = $3Dmol.createViewer("viewer", {{
                            defaultcolors: $3Dmol.rasmolElementColors
                        }});
                        
                        if (!viewer) {{
                            throw new Error("Failed to create viewer");
                        }}
                        console.log("STEP 4: Viewer created successfully");
                        
                        console.log("Loading PDB data...");
                        const pdbData = "{pdb_data_escaped}";
                        console.log("PDB data length:", pdbData.length);
                        
                        console.log("STEP 8: Validating PDB data content...");
                        if (pdbData.length < 100) {{
                            throw new Error("PDB data seems too short (" + pdbData.length + " chars), might be placeholder");
                        }}
                        
                        // Check if it's valid PDB format
                        if (!pdbData.includes("ATOM") && !pdbData.includes("HETATM")) {{
                            throw new Error("PDB data doesn't contain ATOM records");
                        }}
                        
                        console.log("STEP 9: Adding PDB model to viewer...");
                        const numModels = viewer.addModel(pdbData, "pdb");
                        console.log("STEP 10: Model added successfully, number of models:", numModels);
                        
                        // Set style based on user selection
                        console.log("STEP 11: Setting style:", "{style}", "with color scheme:", "{color_scheme}");
                        
                        // Apply color scheme properly based on type
                        var colorScheme = "{color_scheme}";
                        var styleConfig = {{}};
                        
                        if (colorScheme === "spectrum") {{
                            styleConfig = {{spectrum: true}};
                        }} else if (colorScheme === "chain") {{
                            styleConfig = {{colorscheme: "chain"}};
                        }} else if (colorScheme === "element") {{
                            styleConfig = {{colorscheme: $3Dmol.elementColors.rasmol}};
                        }} else if (colorScheme === "residue") {{
                            styleConfig = {{colorscheme: "residue"}};
                        }} else {{
                            styleConfig = {{colorscheme: colorScheme}};
                        }}
                        
                        var mainStyle = "{style}";
                        var finalStyleConfig = {{}};
                        finalStyleConfig[mainStyle] = styleConfig;
                        
                        viewer.setStyle({{}}, finalStyleConfig);
                        console.log("STEP 12: Style applied successfully with config:", finalStyleConfig);
                        
                        // Add ligands if requested
                        {self._get_ligand_style_js(show_ligands) if show_ligands else "// No ligands requested"}
                        
                        // Set initial view and render first
                        console.log("STEP 13: Setting initial view...");
                        viewer.zoomTo();
                        console.log("STEP 14: Zooming complete, now rendering...");
                        
                        viewer.render();
                        console.log("STEP 15: Initial rendering complete!");
                        
                        // Add interaction lines after initial render
                        setTimeout(function() {{
                            console.log("STEP 16: Adding interaction lines...");
                            {self._get_interaction_lines_js(interaction_lines, pdb_data, interaction_visibility)}
                            viewer.render();
                            console.log("STEP 17: Interaction lines added and rendered!");
                        }}, 500);
                        
                        hideLoading();
                        
                        // Clear the timeout since we succeeded
                        if (initializationTimeout) {{
                            clearTimeout(initializationTimeout);
                        }}
                        
                        console.log("STEP 18: 3Dmol viewer initialization FULLY COMPLETE!");
                        
                    }} catch (error) {{
                        console.error("FAILED at step:", error.message);
                        console.error("Full error:", error);
                        showError(error.message || "Unknown error occurred");
                    }}
                }}
                
                // Initialize viewer when page loads
                window.addEventListener('load', function() {{
                    // Wait a bit for 3Dmol to load if needed
                    if (typeof $3Dmol === 'undefined') {{
                        console.log("3Dmol not ready yet, waiting...");
                        let checkCount = 0;
                        const checkInterval = setInterval(function() {{
                            checkCount++;
                            if (typeof $3Dmol !== 'undefined') {{
                                console.log("3Dmol ready after", checkCount * 100, "ms");
                                clearInterval(checkInterval);
                                initialize3Dmol();
                            }} else if (checkCount >= 50) {{ // 5 second max wait
                                clearInterval(checkInterval);
                                showError("3Dmol.js library failed to load from CDN");
                            }}
                        }}, 100);
                    }} else {{
                        initialize3Dmol();
                    }}
                }});
                function resetView() {{
                    viewer.zoomTo();
                    viewer.render();
                }}
                
                function toggleSpin() {{
                    if (spinning) {{
                        viewer.spin(false);
                        spinning = false;
                    }} else {{
                        viewer.spin(true);
                        spinning = true;
                    }}
                }}
                
                function centerStructure() {{
                    viewer.center();
                    viewer.render();
                }}
                
                function toggleFullscreen() {{
                    const elem = document.getElementById('viewer');
                    if (!document.fullscreenElement) {{
                        elem.requestFullscreen().catch(err => {{
                            console.log('Error attempting to enable fullscreen:', err.message);
                        }});
                    }} else {{
                        document.exitFullscreen();
                    }}
                }}
                
                // Interaction highlighting
                function highlightInteraction(residue1, residue2) {{
                    viewer.addStyle({{resi: residue1}}, {{stick: {{colorscheme: "redCarbon"}}}});
                    viewer.addStyle({{resi: residue2}}, {{stick: {{colorscheme: "redCarbon"}}}});
                    viewer.render();
                }}
                
                // Click handler for residues
                viewer.setClickCallback(function(atom, viewer) {{
                    if (atom) {{
                        let info = `Residue: ${{atom.resn}}${{atom.resi}} Chain: ${{atom.chain}}`;
                        alert(info);
                    }}
                }});
                
            </script>
        </body>
        </html>
        """
        
        return html_template
    
    def _get_pdb_data(self, pdb_id: str) -> str:
        """Get PDB data for viewer."""
        try:
            # First try to load from the structure that was already analyzed
            structure = self.pdb_handler.load_structure(pdb_id, assembly="biological")
            if structure:
                # Convert structure back to PDB format using Bio.PDB
                from io import StringIO
                from Bio.PDB import PDBIO
                
                pdb_io = PDBIO()
                pdb_io.set_structure(structure)
                
                pdb_string = StringIO()
                pdb_io.save(pdb_string)
                pdb_content = pdb_string.getvalue()
                
                if pdb_content and len(pdb_content) > 100:  # Ensure it's not just header
                    logger.info(f"Successfully loaded PDB structure for {pdb_id} ({len(pdb_content)} characters)")
                    return pdb_content
            
            # Fallback: Try to get from cache
            if hasattr(self.pdb_handler, 'cache_manager') and self.pdb_handler.cache_manager:
                pdb_content = self.pdb_handler.cache_manager.get_pdb_file(pdb_id, "biological")
                if pdb_content and len(pdb_content) > 100:
                    logger.info(f"Successfully loaded PDB from cache for {pdb_id}")
                    return pdb_content
            
            # Last resort: Download from RCSB
            pdb_content = self.pdb_handler._download_pdb_file(pdb_id, "biological")
            if pdb_content and len(pdb_content) > 100:
                logger.info(f"Successfully downloaded PDB for {pdb_id}")
                return pdb_content
            else:
                logger.warning(f"Using placeholder structure for {pdb_id} - failed to get real PDB data")
                # Return a placeholder if download fails
                return f"HEADER    PLACEHOLDER FOR {pdb_id}\nATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00 20.00           C\nEND"
                
        except Exception as e:
            logger.error(f"Failed to get PDB data for {pdb_id}: {e}")
            # Return a minimal placeholder structure
            return f"HEADER    PLACEHOLDER FOR {pdb_id}\nATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00 20.00           C\nEND"
    
    def _create_interaction_lines(self, 
                                 analysis_result: Dict[str, Any],
                                 interaction_visibility: Dict[str, bool]) -> List[Dict[str, Any]]:
        """Create line data for visualizing interactions."""
        lines = []
        
        interactions = analysis_result.get('interactions', {})
        
        for interaction_type, interaction_list in interactions.items():
            if not interaction_visibility.get(interaction_type, False):
                continue
            
            color = self.viewer_config.interaction_colors.get(interaction_type, "#888888")
            
            for interaction in interaction_list:
                # Extract residue information for coordinate lookup
                residue1 = interaction.get('residue1', '')
                residue2 = interaction.get('residue2', '')
                chain1 = interaction.get('chain1', '')
                chain2 = interaction.get('chain2', '')
                
                # Create line data with residue info for coordinate extraction in JS
                line_data = {
                    'type': interaction_type,
                    'color': color,
                    'residue1': residue1,
                    'residue2': residue2,
                    'chain1': chain1,
                    'chain2': chain2,
                    'distance': interaction.get('distance', 0),
                    'strength': interaction.get('strength', 'moderate'),
                    'label': f"{interaction_type}: {residue1}({chain1}) - {residue2}({chain2})"
                }
                lines.append(line_data)
        
        return lines
    
    def _get_ligand_style_js(self, show_ligands: bool) -> str:
        """Generate JavaScript for ligand styling."""
        if not show_ligands:
            return ""
        
        return """
        // Style ligands
        viewer.addStyle({hetflag: true}, {
            stick: {
                colorscheme: "yellowCarbon",
                radius: 0.3
            }
        });
        """
    
    def _get_interaction_lines_js(self, interaction_lines: List[Dict[str, Any]], pdb_data: str, interaction_visibility: Dict[str, bool]) -> str:
        """Generate JavaScript for interaction lines."""
        if not interaction_lines:
            return ""
        
        # Filter interactions based on visibility settings
        filtered_interactions = []
        for interaction in interaction_lines:
            interaction_type = interaction.get('type', '')
            if interaction_visibility.get(interaction_type, True):
                filtered_interactions.append(interaction)
        
        if not filtered_interactions:
            return ""
        
        # Clean, production-ready interaction visualization
        js_code = f"""
        // Function to find residue center coordinates
        function findResidue(pdbData, resName, resNum, chain) {{
            var atoms = [];
            var lines = pdbData.split('\\n');
            
            for (var i = 0; i < lines.length; i++) {{
                var line = lines[i];
                if (line.startsWith('ATOM') || line.startsWith('HETATM')) {{
                    var lineResNum = line.substring(22, 26).trim();
                    var lineResName = line.substring(17, 20).trim();
                    var lineChain = line.substring(21, 22).trim();
                    
                    if (lineChain === chain && lineResNum === resNum && lineResName === resName) {{
                        var x = parseFloat(line.substring(30, 38));
                        var y = parseFloat(line.substring(38, 46));
                        var z = parseFloat(line.substring(46, 54));
                        if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {{
                            atoms.push({{x: x, y: y, z: z}});
                        }}
                    }}
                }}
            }}
            
            if (atoms.length > 0) {{
                var center = {{x: 0, y: 0, z: 0}};
                for (var j = 0; j < atoms.length; j++) {{
                    center.x += atoms[j].x;
                    center.y += atoms[j].y;
                    center.z += atoms[j].z;
                }}
                center.x /= atoms.length;
                center.y /= atoms.length;
                center.z /= atoms.length;
                return center;
            }}
            return null;
        }}
        
        // Use the already decoded PDB data that's available in the viewer scope
        """
        
        # Add each interaction line
        for i, line in enumerate(filtered_interactions):
            # Parse residue info safely
            res1_name = line['residue1'][:3] if len(line['residue1']) >= 3 else line['residue1']
            res1_num = line['residue1'][3:] if len(line['residue1']) > 3 else "1"
            res2_name = line['residue2'][:3] if len(line['residue2']) >= 3 else line['residue2']
            res2_num = line['residue2'][3:] if len(line['residue2']) > 3 else "1"
            
            # Prepare data for JavaScript (escape and format properly)
            interaction_type = line['type'].upper()
            distance = line.get('distance', 'N/A')
            strength = line.get('strength', 'N/A')
            
            js_code += f"""
        try {{
            var coord1 = findResidue(pdbData, "{res1_name}", "{res1_num}", "{line['chain1']}");
            var coord2 = findResidue(pdbData, "{res2_name}", "{res2_num}", "{line['chain2']}");
            
            if (coord1 && coord2) {{
                // Calculate distance for verification
                var actualDistance = Math.sqrt(
                    Math.pow(coord2.x - coord1.x, 2) + 
                    Math.pow(coord2.y - coord1.y, 2) + 
                    Math.pow(coord2.z - coord1.z, 2)
                ).toFixed(2);
                
                // Store interaction data for hover lookup
                addInteractionData(
                    "{interaction_type}",
                    "{res1_name}{res1_num}",
                    "{line['chain1']}",
                    "{res2_name}{res2_num}",
                    "{line['chain2']}",
                    parseFloat(actualDistance),
                    []  // angles - can be added later
                );
                
                // Add interaction line
                viewer.addCylinder({{
                    start: coord1,
                    end: coord2,
                    radius: 0.08,
                    color: "{line['color']}",
                    alpha: 0.7
                }});
                
                // Add hover-enabled label at midpoint
                var midpoint = {{
                    x: (coord1.x + coord2.x) / 2,
                    y: (coord1.y + coord2.y) / 2,
                    z: (coord1.z + coord2.z) / 2
                }};
                
                // Create interactive sphere for hover and click functionality  
                const sphere = viewer.addSphere({{
                    center: midpoint,
                    radius: 0.2,
                    color: "{line['color']}",
                    alpha: 0.8
                }});
                
                // Add a label with hover functionality
                const label = viewer.addLabel("{res1_name}{res1_num}-{res2_name}{res2_num}", {{
                    position: midpoint,
                    fontSize: 14,
                    fontColor: 'white',
                    backgroundColor: 'rgba(0,0,0,0.9)',
                    showBackground: true,
                    alignment: 'center',
                    borderColor: "{line['color']}",
                    borderWidth: 2
                }});
                
                // Add hover functionality to the label area
                viewer.setHoverable({{}}, true, function(atom, viewer, event, container) {{
                    // Check if we're hovering near this interaction point/label
                    if (atom && atom.x !== undefined && atom.y !== undefined && atom.z !== undefined) {{
                        const dx = Math.abs(atom.x - midpoint.x);
                        const dy = Math.abs(atom.y - midpoint.y);
                        const dz = Math.abs(atom.z - midpoint.z);
                        const distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
                        
                        if (distance < 1.5) {{ // Within 1.5 Angstroms of label/interaction center
                            // Remove any existing hover tooltip
                            const existingHover = document.getElementById('label-hover-tooltip');
                            if (existingHover) existingHover.remove();
                            
                            // Create hover tooltip for this specific interaction
                            const hoverDiv = document.createElement('div');
                            hoverDiv.id = 'label-hover-tooltip';
                            hoverDiv.style.cssText = `
                                position: fixed;
                                background: rgba(0,0,0,0.95);
                                color: white;
                                padding: 12px 16px;
                                border-radius: 8px;
                                font-size: 13px;
                                font-family: 'Segoe UI', Arial, sans-serif;
                                border: 2px solid {line['color']};
                                z-index: 10000;
                                pointer-events: none;
                                left: ${{event.pageX + 15}}px;
                                top: ${{event.pageY - 10}}px;
                                max-width: 300px;
                                box-shadow: 0 6px 12px rgba(0,0,0,0.4);
                            `;
                            
                            hoverDiv.innerHTML = `
                                <div style="font-weight: bold; color: {line['color']}; margin-bottom: 8px; font-size: 14px;">
                                    🔗 {interaction_type}
                                </div>
                                <div style="margin-bottom: 4px;">
                                    <span style="color: #FFC107; font-weight: bold;">Acceptor:</span> {res1_name}{res1_num} (Chain {line['chain1']})
                                </div>
                                <div style="margin-bottom: 4px;">
                                    <span style="color: #FF5722; font-weight: bold;">Donor:</span> {res2_name}{res2_num} (Chain {line['chain2']})
                                </div>
                                <div style="margin-bottom: 4px;">
                                    <span style="color: #2196F3; font-weight: bold;">Distance:</span> " + actualDistance + " Å
                                </div>
                                <div style="margin-bottom: 8px;">
                                    <span style="color: #9C27B0; font-weight: bold;">Strength:</span> {strength}
                                </div>
                                <div style="font-size: 11px; color: #aaa; text-align: center;">
                                    Click for more details
                                </div>
                            `;
                            
                            document.body.appendChild(hoverDiv);
                        }}
                    }}
                }}, function(atom, viewer) {{
                    // Mouse out - remove hover tooltip
                    const existingHover = document.getElementById('label-hover-tooltip');
                    if (existingHover) existingHover.remove();
                }});
                
                // Add click handler for detailed information
                viewer.setClickable({{}}, true, function(atom, viewer, event) {{
                    // Check if click is near this interaction point
                    if (atom && atom.x !== undefined) {{
                        const dx = Math.abs(atom.x - midpoint.x);
                        const dy = Math.abs(atom.y - midpoint.y);
                        const dz = Math.abs(atom.z - midpoint.z);
                        const distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
                        
                        if (distance < 1.5) {{ // Within 1.5 Angstrom of interaction center
                            const details = "🔗 {interaction_type}\\n" +
                                          "═══════════════════════════════════════\\n" +
                                          "Acceptor: {res1_name}{res1_num} (Chain {line['chain1']})\\n" +
                                          "Donor: {res2_name}{res2_num} (Chain {line['chain2']})\\n" +
                                          "Distance: " + actualDistance + " Å\\n" +
                                          "Strength: {strength}\\n" +
                                          "═══════════════════════════════════════";
                            alert(details);
                        }}
                    }}
                }});
                
                }}
            }}
        }} catch (error) {{
            console.log("Error adding interaction:", error);
        }}
"""
        
        return js_code
    
    def _render_residue_search(self, analysis_result: Dict[str, Any]):
        """Render residue search and selection interface."""
        st.subheader("🔍 Residue Navigator")
        
        col1, col2 = st.columns(2)
        
        with col1:
            search_residue = st.text_input(
                "Search Residue:",
                placeholder="e.g., ALA123 or A:123",
                help="Enter residue name and number, optionally with chain"
            )
            
            if search_residue:
                st.info(f"Would jump to {search_residue} in viewer")
        
        with col2:
            # Show hotspots
            hotspots = analysis_result.get('hotspots', [])
            if hotspots:
                selected_hotspot = st.selectbox(
                    "Jump to Hotspot:",
                    [f"{h['residue']} ({h['interaction_count']} interactions)" for h in hotspots[:10]]
                )
                
                if st.button("Go to Hotspot"):
                    st.info(f"Would jump to {selected_hotspot} in viewer")
        
        # Interaction summary table
        if st.checkbox("Show Interaction Summary"):
            self._render_interaction_table(analysis_result)
    
    def _render_interaction_table(self, analysis_result: Dict[str, Any]):
        """Render a summary table of interactions."""
        interactions = analysis_result.get('interactions', {})
        
        summary_data = []
        for interaction_type, interaction_list in interactions.items():
            summary_data.append({
                'Interaction Type': self._get_interaction_display_name(interaction_type),
                'Count': len(interaction_list),
                'Avg Distance': f"{np.mean([i.get('distance', 0) for i in interaction_list]):.2f} Å" if interaction_list else "N/A"
            })
        
        if summary_data:
            import pandas as pd
            df = pd.DataFrame(summary_data)
            st.dataframe(df, use_container_width=True)
    
    def _get_interaction_display_name(self, interaction_type: str) -> str:
        """Get human-readable name for interaction type."""
        display_names = {
            'hydrogenbond': 'Hydrogen Bonds',
            'halogenbond': 'Halogen Bonds',
            'pipi': 'π-π Stacking',
            'ionicinteraction': 'Ionic Interactions',
            'hydrophobiccontact': 'Hydrophobic Contacts',
            'chpi': 'C-H···π Interactions',
            'chalcogenbond': 'Chalcogen Bonds',
            'pnictogenbond': 'Pnictogen Bonds',
            'tetrelbond': 'Tetrel Bonds',
            'anionpi': 'Anion-π Interactions',
            'npistar': 'n→π* Interactions',
            'dispersion': 'London Dispersion'
        }
        return display_names.get(interaction_type, interaction_type.replace('_', ' ').title())
    
    def _generate_pymol_session(self, 
                               pdb_id: str,
                               analysis_result: Dict[str, Any],
                               selected_interactions: List[str]) -> bytes:
        """Generate PyMOL session file (.pse)."""
        # This would generate a complete PyMOL session
        # For now, return placeholder
        pymol_script = f"""
# PyMOL session for {pdb_id}
# Generated by Protein Interaction Explorer

# Load structure
fetch {pdb_id}

# Basic styling
cartoon
color spectrum

# Show interactions
{self._generate_pymol_interaction_commands(analysis_result, selected_interactions)}

# Set view
orient
zoom
"""
        
        return pymol_script.encode('utf-8')
    
    def _generate_chimera_script(self,
                                pdb_id: str, 
                                analysis_result: Dict[str, Any],
                                selected_interactions: List[str]) -> str:
        """Generate ChimeraX script (.cxc)."""
        chimera_script = f"""# ChimeraX script for {pdb_id}
# Generated by Protein Interaction Explorer

# Open structure
open {pdb_id}

# Basic styling
cartoon
color bychain

# Show interactions
{self._generate_chimera_interaction_commands(analysis_result, selected_interactions)}

# Set view
view
"""
        
        return chimera_script
    
    def _generate_pymol_interaction_commands(self,
                                           analysis_result: Dict[str, Any],
                                           selected_interactions: List[str]) -> str:
        """Generate PyMOL commands for showing interactions."""
        commands = []
        
        interactions = analysis_result.get('interactions', {})
        
        for interaction_type in selected_interactions:
            if interaction_type not in interactions:
                continue
            
            color = self.viewer_config.interaction_colors.get(interaction_type, "yellow")
            
            for i, interaction in enumerate(interactions[interaction_type]):
                res1 = interaction.get('residue1', '')
                res2 = interaction.get('residue2', '')
                chain1 = interaction.get('chain1', '')
                chain2 = interaction.get('chain2', '')
                
                obj_name = f"{interaction_type}_{i}"
                
                commands.append(f"distance {obj_name}, chain {chain1} and resi {res1}, chain {chain2} and resi {res2}")
                commands.append(f"color {color}, {obj_name}")
        
        return "\n".join(commands)
    
    def _generate_chimera_interaction_commands(self,
                                             analysis_result: Dict[str, Any],
                                             selected_interactions: List[str]) -> str:
        """Generate ChimeraX commands for showing interactions."""
        commands = []
        
        interactions = analysis_result.get('interactions', {})
        
        for interaction_type in selected_interactions:
            if interaction_type not in interactions:
                continue
            
            color = self.viewer_config.interaction_colors.get(interaction_type, "yellow")
            
            for interaction in interactions[interaction_type]:
                res1 = interaction.get('residue1', '')
                res2 = interaction.get('residue2', '')
                chain1 = interaction.get('chain1', '')
                chain2 = interaction.get('chain2', '')
                
                commands.append(f"distance /{chain1}:{res1} /{chain2}:{res2}")
                commands.append(f"color {color} distances")
        
        return "\n".join(commands)
