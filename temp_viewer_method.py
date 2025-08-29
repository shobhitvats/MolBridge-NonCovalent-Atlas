from typing import Any, List, Dict

def _get_interaction_property(interaction: Any, property_name: str, default_value: Any = None) -> Any:
    """
    Get a property from an interaction object, handling both dictionary and object styles.
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

def _create_main_viewer(self, pdb_id: str, pdb_data: str, interaction_lines: List[Dict[str, Any]]) -> str:
    """Create a simple, working main viewer HTML."""
    # Ensure we have PDB data
    if not pdb_data or len(pdb_data) < 100:
        pdb_data = self.pdb_handler.get_pdb_content(pdb_id)
    
    if not pdb_data:
        return f'<div style="color: red; padding: 20px;">❌ Failed to load PDB data for {pdb_id}</div>'
    
    # Set up default values
    interaction_visibility = {_get_interaction_property(interaction, 'type', 'unknown'): True for interaction in interaction_lines}
    
    # Escape PDB data for JavaScript
    pdb_data_escaped = pdb_data.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n').replace('\r', '')
    
    # Create simple, working viewer
    main_viewer = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                #main-viewer {{ width: 100%; height: 500px; border: 2px solid #ddd; border-radius: 8px; position: relative; }}
                #loading {{ position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); font-size: 18px; color: #666; z-index: 100; }}
                .error-message {{ position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); color: red; font-weight: bold; text-align: center; z-index: 100; }}
            </style>
        </head>
        <body>
            <div id="main-viewer">
                <div id="loading">Loading structure...</div>
            </div>
            
            <script>
            let viewer;
            
            function hideLoading() {{
                const loading = document.getElementById('loading');
                if (loading) loading.style.display = 'none';
            }}
            
            function showError(message) {{
                hideLoading();
                const errorDiv = document.createElement('div');
                errorDiv.className = 'error-message';
                errorDiv.innerHTML = '<strong>Error:</strong><br>' + message;
                document.getElementById('main-viewer').appendChild(errorDiv);
            }}
            
            function init() {{
                try {{
                    console.log("Initializing viewer...");
                    
                    if (typeof $3Dmol === 'undefined') {{
                        throw new Error("3Dmol.js not loaded");
                    }}
                    
                    viewer = $3Dmol.createViewer("main-viewer");
                    
                    const pdbData = "{pdb_data_escaped}";
                    if (pdbData.length < 100) {{
                        throw new Error("No PDB data");
                    }}
                    
                    viewer.addModel(pdbData, "pdb");
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: "spectrum"}}}});
                    
                    viewer.setClickable({{}}, true, function(atom) {{
                        if (atom.resn && atom.resi && atom.chain) {{
                            alert(atom.resn + atom.resi + ":" + atom.chain);
                        }}
                    }});
                    
                    viewer.zoomTo();
                    viewer.render();
                    
                    // Add interactions
                    setTimeout(() => {{
                        {self._get_interaction_lines_js(interaction_lines, pdb_data, interaction_visibility)}
                        viewer.render();
                    }}, 500);
                    
                    hideLoading();
                    console.log("Viewer ready!");
                    
                }} catch (error) {{
                    console.error("Error:", error);
                    showError(error.message);
                }}
            }}
            
            // Start initialization
            if (document.readyState === 'loading') {{
                document.addEventListener('DOMContentLoaded', init);
            }} else {{
                setTimeout(init, 100);
            }}
            </script>
        </body>
        </html>
        """
    
    return main_viewer
