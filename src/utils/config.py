"""
Configuration management for Protein Interaction Explorer.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

@dataclass
class InteractionConfig:
    """Configuration for interaction detection parameters."""
    
    # Hydrogen bonds
    hbond_distance_cutoff: float = 3.5
    hbond_angle_cutoff: float = 120.0
    
    # Halogen bonds
    halogen_distance_cutoff: float = 4.0
    halogen_angle_cutoff: float = 140.0
    
    # Chalcogen bonds
    chalcogen_distance_cutoff: float = 4.0
    chalcogen_angle_cutoff: float = 140.0
    
    # Pnictogen bonds
    pnictogen_distance_cutoff: float = 4.0
    pnictogen_angle_cutoff: float = 140.0
    
    # Tetrel bonds
    tetrel_distance_cutoff: float = 4.0
    tetrel_angle_cutoff: float = 140.0
    
    # C-H...π interactions
    ch_pi_distance_cutoff: float = 4.5
    ch_pi_angle_cutoff: float = 90.0
    ch_pi_max_distance: float = 4.5
    ch_pi_min_distance: float = 2.0
    ch_pi_max_angle: float = 90.0
    ch_pi_max_height: float = 2.5
    
    # π-π stacking
    pi_pi_distance_cutoff: float = 5.5
    pi_pi_angle_cutoff: float = 30.0
    
    # Anion-π interactions
    anion_pi_distance_cutoff: float = 5.0
    
    # n→π* interactions
    n_pi_star_distance_cutoff: float = 3.5
    n_pi_star_angle_cutoff: float = 120.0
    
    # London dispersion (heuristic)
    dispersion_distance_cutoff: float = 4.5
    dispersion_min_distance: float = 3.5
    
    # Ionic interactions
    ionic_distance_cutoff: float = 6.0
    
    # Hydrophobic contacts
    hydrophobic_distance_cutoff: float = 5.0

@dataclass
class VisualizationConfig:
    """Configuration for visualization settings."""
    
    default_style: str = "cartoon"
    protein_color: str = "spectrum"
    ligand_color: str = "yellow"
    interaction_colors: Dict[str, str] = field(default_factory=lambda: {
        "hydrogen_bond": "#FF6B6B",
        "halogen_bond": "#4ECDC4", 
        "chalcogen_bond": "#45B7D1",
        "pnictogen_bond": "#96CEB4",
        "tetrel_bond": "#FFEAA7",
        "ch_pi": "#DDA0DD",
        "pi_pi": "#F39C12",
        "anion_pi": "#E74C3C",
        "n_pi_star": "#9B59B6",
        "dispersion": "#95A5A6",
        "ionic": "#3498DB",
        "hydrophobic": "#E67E22"
    })
    background_color: str = "white"
    viewer_width: int = 800
    viewer_height: int = 600

@dataclass
class ProcessingConfig:
    """Configuration for processing and performance settings."""
    
    max_batch_size: int = 200
    max_workers: int = 4
    chunk_size: int = 10
    timeout_seconds: int = 300
    memory_limit_gb: float = 8.0
    cache_expiry_days: int = 7

@dataclass
class AppConfig:
    """Main application configuration."""
    
    # Basic settings
    app_name: str = "Protein Interaction Explorer"
    version: str = "1.0.0"
    debug: bool = False
    
    # Directories
    base_dir: Path = field(default_factory=lambda: Path(__file__).parent.parent)
    cache_dir: Path = field(init=False)
    templates_dir: Path = field(init=False)
    examples_dir: Path = field(init=False)
    
    # API settings
    api_enabled: bool = True
    api_host: str = "0.0.0.0"
    api_port: int = 8000
    
    # External API keys (optional)
    openai_api_key: Optional[str] = None
    anthropic_api_key: Optional[str] = None
    
    # Component configurations
    interactions: InteractionConfig = field(default_factory=InteractionConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    
    # PDB settings
    pdb_base_url: str = "https://files.rcsb.org/download/"
    include_ligands: bool = True
    exclude_waters: bool = True
    default_assembly: str = "biological"  # or "asymmetric"
    
    # Interaction presets
    presets: Dict[str, Dict[str, float]] = field(default_factory=lambda: {
        "conservative": {
            "hbond_distance_cutoff": 3.2,
            "halogen_distance_cutoff": 3.5,
            "pi_pi_distance_cutoff": 4.5,
        },
        "literature_default": {
            "hbond_distance_cutoff": 3.5,
            "halogen_distance_cutoff": 4.0,
            "pi_pi_distance_cutoff": 5.5,
        },
        "exploratory": {
            "hbond_distance_cutoff": 4.0,
            "halogen_distance_cutoff": 4.5,
            "pi_pi_distance_cutoff": 6.0,
        }
    })
    
    def __post_init__(self):
        """Initialize computed fields after object creation."""
        self.cache_dir = self.base_dir / "cache"
        self.templates_dir = self.base_dir / "templates" 
        self.examples_dir = self.base_dir / "examples"
        
        # Create directories if they don't exist
        self.cache_dir.mkdir(exist_ok=True)
        self.templates_dir.mkdir(exist_ok=True)
        self.examples_dir.mkdir(exist_ok=True)
        
        # Load API keys from environment
        self.openai_api_key = os.getenv("OPENAI_API_KEY")
        self.anthropic_api_key = os.getenv("ANTHROPIC_API_KEY")

def load_config() -> AppConfig:
    """Load application configuration."""
    return AppConfig()

def get_interaction_types() -> List[str]:
    """Get list of all supported interaction types."""
    return [
        "hydrogenbond",
        "halogenbond", 
        "chalcogenbond",
        "pnictogenbond",
        "tetrelbond",
        "chpi",
        "pipi",
        "anionpi", 
        "npistar",
        "dispersion",
        "ionicinteraction",
        "hydrophobiccontact"
    ]

def get_interaction_display_names() -> Dict[str, str]:
    """Get human-readable names for interaction types."""
    return {
        "hydrogenbond": "Hydrogen Bonds",
        "halogenbond": "Halogen Bonds",
        "chalcogenbond": "Chalcogen Bonds", 
        "pnictogenbond": "Pnictogen Bonds",
        "tetrelbond": "Tetrel Bonds",
        "chpi": "C-H···π Interactions",
        "pipi": "π-π Stacking",
        "anionpi": "Anion-π Interactions",
        "npistar": "n→π* Interactions", 
        "dispersion": "London Dispersion",
        "ionicinteraction": "Ionic Interactions",
        "hydrophobiccontact": "Hydrophobic Contacts"
    }
