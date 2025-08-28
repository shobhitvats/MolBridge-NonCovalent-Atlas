#!/usr/bin/env python3

def create_simple_viewer(pdb_data_b64):
    """Create a simple 3Dmol viewer for testing"""
    
    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        <style>
            #viewer {{ width: 100%; height: 600px; border: 1px solid #ddd; }}
            .status {{ position: absolute; top: 10px; left: 10px; padding: 8px; border-radius: 4px; z-index: 100; }}
            .loading {{ background: rgba(255,255,0,0.8); }}
            .success {{ background: rgba(0,255,0,0.8); }}
            .error {{ background: rgba(255,0,0,0.8); }}
        </style>
    </head>
    <body>
        <div id="viewer">
            <div id="status" class="status loading">Loading structure...</div>
        </div>
        
        <script>
        (function() {{
            const statusEl = document.getElementById('status');
            
            function setStatus(msg, type) {{
                statusEl.textContent = msg;
                statusEl.className = 'status ' + type;
            }}
            
            function init() {{
                try {{
                    setStatus('Checking 3Dmol...', 'loading');
                    
                    if (typeof $3Dmol === 'undefined') {{
                        throw new Error('3Dmol.js not loaded');
                    }}
                    
                    setStatus('Creating viewer...', 'loading');
                    const viewer = $3Dmol.createViewer('viewer');
                    
                    setStatus('Decoding PDB data...', 'loading');
                    const pdbData = atob('{pdb_data_b64}');
                    
                    setStatus('Loading structure...', 'loading');
                    viewer.addModel(pdbData, 'pdb');
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'spectrum'}}}});
                    viewer.zoomTo();
                    viewer.render();
                    
                    setStatus('✅ Success!', 'success');
                    setTimeout(() => statusEl.style.display = 'none', 2000);
                    
                }} catch (error) {{
                    setStatus('❌ Error: ' + error.message, 'error');
                    console.error('Error:', error);
                }}
            }}
            
            document.addEventListener('DOMContentLoaded', function() {{
                if (typeof $3Dmol !== 'undefined') {{
                    init();
                }} else {{
                    let attempts = 0;
                    const check = setInterval(() => {{
                        if (typeof $3Dmol !== 'undefined') {{
                            clearInterval(check);
                            init();
                        }} else if (++attempts > 50) {{
                            clearInterval(check);
                            setStatus('❌ 3Dmol.js failed to load', 'error');
                        }}
                    }}, 100);
                }}
            }});
        }})();
        </script>
    </body>
    </html>
    """
    
    return html

# Test the function
if __name__ == "__main__":
    import base64
    
    # Create some test PDB data
    test_pdb = """
ATOM      1  N   ALA A   1      20.154  16.967  27.462  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.277  26.797  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.708  17.018  26.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.755  18.177  27.297  1.00 20.00           O  
ATOM      5  CB  ALA A   1      19.370  16.047  25.323  1.00 20.00           C  
END
    """
    
    pdb_b64 = base64.b64encode(test_pdb.encode()).decode()
    html = create_simple_viewer(pdb_b64)
    
    print("Generated HTML length:", len(html))
    print("Base64 data length:", len(pdb_b64))
