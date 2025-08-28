# Protein-Interaction-Explorer

A comprehensive, production-ready Streamlit application for analyzing noncovalent interactions in protein structures. This tool provides state-of-the-art detection algorithms for 11 different types of molecular interactions, interactive 3D visualization, and extensive reporting capabilities.

## 🚀 Features

### Interaction Detection
- **Hydrogen Bonds**: Conventional, low-barrier, and C5-type hydrogen bonds
- **Halogen Bonds**: Cl, Br, I, and F-mediated interactions with directional sigma holes
- **π-π Stacking**: Face-to-face and edge-to-face aromatic interactions
- **Ionic Interactions**: Salt bridges between charged residues
- **Hydrophobic Contacts**: Van der Waals interactions between nonpolar groups
- **C-H···π Interactions**: Weak hydrogen bonds to aromatic systems
- **Chalcogen Bonds**: S, Se, Te-mediated interactions
- **Pnictogen Bonds**: N, P, As-mediated interactions
- **Tetrel Bonds**: C, Si, Ge-mediated interactions
- **Anion-π Interactions**: Negative charges interacting with π-systems
- **n→π* Interactions**: Orbital overlap interactions
- **London Dispersion**: Heuristic detection of dispersion forces

### Visualization & Analysis
- **Interactive 3D Viewer**: py3Dmol-based structure visualization with interaction overlays
- **Comprehensive Plots**: Heatmaps, distribution plots, Ramachandran plots with hotspot highlighting
- **Comparative Analysis**: Side-by-side comparison of multiple structures
- **Hotspot Detection**: Identification of interaction-rich regions
- **Export Capabilities**: PyMOL and ChimeraX session generation

### Data Management
- **Intelligent Caching**: Persistent storage with 7-day expiration
- **Session Management**: Save, load, and compare analysis sessions
- **Batch Processing**: Parallel analysis of multiple structures
- **Quality Control**: Structure validation and completeness assessment

### Reports & Export
- **PDF Reports**: Comprehensive analysis reports with figures and tables
- **Excel Workbooks**: Multi-sheet data export with summary statistics
- **PowerPoint Presentations**: Publication-ready slide decks
- **CSV Data**: Raw interaction data for further analysis
- **LaTeX Export**: Publication-ready snippets for manuscripts

### REST API
- **Programmatic Access**: Full REST API with FastAPI backend
- **Background Processing**: Asynchronous job execution with progress tracking
- **Multiple Formats**: JSON, CSV, PDF, Excel output formats
- **File Upload**: Custom structure analysis support

## 📋 Requirements

### Core Dependencies
```
streamlit>=1.28.0
biopython>=1.83
biopandas>=0.4.0
numpy>=1.24.0
scipy>=1.12.0
pandas>=2.2.0
plotly>=5.18.0
py3dmol>=2.0.4
matplotlib>=3.8.0
seaborn>=0.13.0
```

### Extended Dependencies
```
fastapi>=0.104.0
uvicorn>=0.24.0
reportlab>=4.0.0
openpyxl>=3.1.0
python-pptx>=0.6.23
diskcache>=5.6.0
scikit-learn>=1.3.0
networkx>=3.2
```

## 🛠️ Installation

### Quick Start with Docker (Recommended)
```bash
git clone https://github.com/yourusername/protein-interaction-explorer.git
cd protein-interaction-explorer
docker build -t protein-explorer .
docker run -p 8501:8501 -p 8000:8000 protein-explorer
```

### Manual Installation
```bash
git clone https://github.com/yourusername/protein-interaction-explorer.git
cd protein-interaction-explorer
pip install -r requirements.txt
```

## 🚀 Usage

### Streamlit Web Interface
```bash
streamlit run server.py
```
Access the application at `http://localhost:8501`

### REST API Server
```bash
cd src
python api/endpoints.py
```
API documentation available at `http://localhost:8000/docs`

### Command Line Analysis
```python
from src.analysis.batch_processor import BatchProcessor
from src.utils.config import AppConfig

config = AppConfig()
processor = BatchProcessor(config)
results = processor.process_batch(['1LYZ', '2LYZ'])
```

## 📖 User Guide

### Basic Workflow

1. **Input Structures**
   - Enter PDB IDs (e.g., "1LYZ, 2LYZ, 3LYZ")
   - Upload custom PDB/CIF files
   - Paste structure data directly

2. **Configure Analysis**
   - Choose interaction types to analyze
   - Select parameter preset (Conservative/Literature Default/Exploratory)
   - Customize individual parameters if needed

3. **Run Analysis**
   - Click "Start Analysis" to begin processing
   - Monitor progress in real-time
   - Review results as they complete

4. **Explore Results**
   - **Structure Tab**: Interactive 3D visualization with interaction overlays
   - **Analysis Tab**: Detailed interaction tables and statistics
   - **Visualization Tab**: Heatmaps, distribution plots, and comparative charts
   - **Reports Tab**: Generate and download comprehensive reports

### Configuration Presets

#### Conservative
- Stricter distance and angle cutoffs
- Higher confidence interactions only
- Suitable for high-confidence structural biology

#### Literature Default (Recommended)
- Standard literature-based parameters
- Balanced sensitivity and specificity
- Most commonly used in publications

#### Exploratory
- Relaxed cutoffs for comprehensive screening
- Higher sensitivity, lower specificity
- Useful for discovering weak interactions

### Advanced Features

#### Session Management
```python
# Save current session
session_manager.save_session("my_analysis", {
    "pdb_ids": ["1LYZ"],
    "results": results,
    "config": config.to_dict()
})

# Load previous session
session_data = session_manager.load_session("my_analysis")
```

#### Custom Configuration
```python
# Customize hydrogen bond parameters
config.interaction_config.hbond_distance_cutoff = 3.2
config.interaction_config.hbond_angle_cutoff = 25.0

# Set processing options
config.processing_config.max_workers = 8
config.processing_config.use_cache = True
```

#### Batch Processing
```python
# Process multiple structures with custom parameters
results = batch_processor.process_batch(
    pdb_ids=["1LYZ", "2LYZ", "3LYZ"],
    interaction_types=["hydrogen_bond", "pi_pi", "ionic"],
    use_cache=True,
    include_metadata=True
)
```

## 🔬 Scientific Background

### Interaction Detection Algorithms

All interaction detection algorithms implement literature-validated geometric criteria:

- **Distance cutoffs** based on sum of van der Waals radii plus tolerance
- **Angle constraints** for directional interactions (hydrogen bonds, halogen bonds)
- **Planarity requirements** for π-π stacking interactions
- **Charge-distance relationships** for ionic interactions

### Validation & Quality Control

- **Structure completeness** assessment
- **Missing atom detection** and handling
- **Coordinate validation** and outlier detection
- **Chain break identification**
- **Alternative conformer handling**

### Performance Optimization

- **KD-tree spatial indexing** for efficient neighbor searches
- **Multiprocessing** for parallel structure analysis
- **Intelligent caching** to avoid redundant calculations
- **Memory-efficient** data structures for large complexes

## 📊 Example Results

### Hydrogen Bond Analysis
```
Structure: 1LYZ (Lysozyme)
Total hydrogen bonds: 89
Average distance: 2.85 ± 0.23 Å
Average angle: 162.3 ± 12.7°
Strength distribution:
  Strong: 45 bonds (50.6%)
  Moderate: 32 bonds (36.0%)
  Weak: 12 bonds (13.5%)
```

### Interaction Hotspots
```
Top interaction hotspots:
1. ASP52A - 12 interactions (4 HB, 3 ionic, 5 others)
2. GLU35A - 9 interactions (3 HB, 2 ionic, 4 others)
3. TRP62A - 8 interactions (2 HB, 4 π-π, 2 others)
```

## 🌐 REST API Usage

### Start Analysis Job
```bash
curl -X POST "http://localhost:8000/analyze" 
     -H "Content-Type: application/json" 
     -d '{
       "pdb_ids": ["1LYZ", "2LYZ"],
       "interaction_types": ["hydrogen_bond", "pi_pi"],
       "config_preset": "literature_default"
     }'
```

### Check Job Status
```bash
curl "http://localhost:8000/jobs/{job_id}"
```

### Download Results
```bash
curl "http://localhost:8000/results/{job_id}" > results.json
```

### Generate Reports
```bash
curl -X POST "http://localhost:8000/reports" 
     -H "Content-Type: application/json" 
     -d '{
       "pdb_ids": ["1LYZ"],
       "format": "pdf",
       "include_methodology": true
     }' 
     --output report.pdf
```

## 🧪 Development & Testing

### Running Tests
```bash
python -m pytest tests/ -v
```

### Development Setup
```bash
pip install -r requirements-dev.txt
pre-commit install
```

### Contributing
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📚 References

### Key Publications
1. Jeffrey, G.A. (1997). "An Introduction to Hydrogen Bonding"
2. Desiraju, G.R. et al. (2013). "Definition of the halogen bond" Pure Appl. Chem. 85, 1711-1713
3. McGaughey, G.B. et al. (1998). "π-Stacking interactions" J. Biol. Chem. 273, 15458-15463
4. Kumar, S. & Nussinov, R. (2002). "Close-range electrostatic interactions in proteins" ChemBioChem 3, 604-617

### Algorithm Sources
- Hydrogen bonds: Baker & Hubbard (1984) criteria with modern extensions
- Halogen bonds: Auffinger et al. (2004) geometric definitions
- π-π stacking: Janiak (2000) classification system
- Ionic interactions: Barlow & Thornton (1983) distance criteria

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Support

- **Documentation**: [User Guide](docs/user_guide.md) | [API Reference](docs/api_reference.md)
- **Examples**: [Jupyter Notebooks](examples/) | [Example Scripts](examples/scripts/)
- **Issues**: [GitHub Issues](https://github.com/yourusername/protein-interaction-explorer/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/protein-interaction-explorer/discussions)

## 🏆 Citation

If you use Protein Interaction Explorer in your research, please cite:

```bibtex
@software{protein_interaction_explorer,
  title={Protein Interaction Explorer: Comprehensive Analysis of Noncovalent Interactions in Protein Structures},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/protein-interaction-explorer}
}
```

## 🔄 Version History

### v1.0.0 (Current)
- Initial release with 11 interaction types
- Complete Streamlit interface
- REST API implementation
- Comprehensive reporting system
- Interactive 3D visualization
- Batch processing capabilities

### Planned Features (v1.1.0)
- Machine learning-based interaction prediction
- Enhanced visualization options
- Additional export formats
- Performance optimizations
- Extended API endpoints

---

**Built with ❤️ for the structural biology community**