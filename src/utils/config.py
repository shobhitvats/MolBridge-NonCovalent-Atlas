"""
Configuration management for Protein Interaction Explorer.
"""

from dataclasses import dataclass, field, asdict
import yaml
from pathlib import Path
from typing import Dict, List, Optional
import os
try:  # Optional dependency: PyYAML
    import yaml  # type: ignore
except Exception:  # pragma: no cover - handled gracefully
    yaml = None  # type: ignore
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

import hashlib
import json

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
    # Replace single angle with explicit theta range + phi (delta) max deviation
    chalcogen_theta_min: float = 115.0
    chalcogen_theta_max: float = 155.0
    chalcogen_phi_max: float = 50.0
    
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

    # Newly added detectors (explicit parameters)
    cation_pi_distance_cutoff: float = 6.0
    salt_bridge_distance_cutoff: float = 4.0
    salt_bridge_centroid_cutoff: float = 5.0
    sulfur_pi_distance_cutoff: float = 6.0
    sulfur_pi_max_perp: float = 3.0
    metal_coordination_primary_cutoff: float = 2.6
    metal_coordination_extended_cutoff: float = 3.0

    # Internal version & hash fields (auto-managed)
    _version: int = 1  # bump manually when semantic meaning of any param changes
    _param_hash: str = field(default="", init=False, repr=False)

    def compute_hash(self) -> str:
        """Compute a stable hash of all public interaction parameters.

        Excludes private / cache fields (those starting with underscore).
        Produces a short 10-char hex digest for compact cache keys.
        """
        data = {k: v for k, v in asdict(self).items() if not k.startswith('_')}
        # Stable ordering via sorted keys
        payload = json.dumps(data, sort_keys=True, separators=(",", ":"))
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:10]
        self._param_hash = digest
        return digest

    @property
    def param_hash(self) -> str:
        if not self._param_hash:
            return self.compute_hash()
        return self._param_hash

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
        ,"cation_pi": "#1ABC9C"
        ,"salt_bridge": "#2E86C1"
        ,"sulfur_pi": "#AF7AC5"
        ,"metal_coordination": "#17A589"
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

    # Feature toggles (protein interaction focused)
    enable_residue_profiles: bool = True
    enable_interface_analysis: bool = True
    enable_outlier_detection: bool = True
    enable_provenance_panel: bool = True
    enable_motif_detection: bool = False  # experimental

    # Structural quality / annotation feature toggles
    enable_secondary_structure: bool = True
    enable_sasa_bsa: bool = True
    enable_geometry_quality: bool = True
    enable_disulfide_analysis: bool = True
    enable_pocket_detection: bool = True
    enable_conservation: bool = False  # will start as stub / off by default
    enable_pi_pi_refinement: bool = True
    enable_hbond_subtypes: bool = True  # hydrogen bond subtype classification

    # Automation flags
    auto_compute_extensions: bool = True  # compute all enabled extensions after analysis
    force_package_extensions: bool = True  # ensure extensions are computed during package export
    
    # PDB settings
    pdb_base_url: str = "https://files.rcsb.org/download/"
    include_ligands: bool = True
    exclude_waters: bool = True
    default_assembly: str = "biological"  # or "asymmetric"
    
    # Interaction presets
    presets: Dict[str, Dict[str, float]] = field(default_factory=lambda: {
        # All-inclusive presets covering every interaction parameter for consistency across UI sections
        "conservative": {
            "hbond_distance_cutoff": 3.2, "hbond_angle_cutoff": 130.0,
            "halogen_distance_cutoff": 3.5, "halogen_angle_cutoff": 150.0,
            "chalcogen_distance_cutoff": 3.6, "chalcogen_theta_min": 120.0, "chalcogen_theta_max": 150.0, "chalcogen_phi_max": 40.0,
            "pnictogen_distance_cutoff": 3.6, "pnictogen_angle_cutoff": 150.0,
            "tetrel_distance_cutoff": 3.6, "tetrel_angle_cutoff": 165.0,
            "ch_pi_distance_cutoff": 4.2, "ch_pi_angle_cutoff": 85.0,
            "ch_pi_min_distance": 2.0, "ch_pi_max_distance": 4.5, "ch_pi_max_angle": 90.0, "ch_pi_max_height": 2.5,
            "pi_pi_distance_cutoff": 5.0, "pi_pi_angle_cutoff": 25.0,
            "anion_pi_distance_cutoff": 4.6,
            "n_pi_star_distance_cutoff": 3.3, "n_pi_star_angle_cutoff": 125.0,
            "dispersion_min_distance": 3.5, "dispersion_distance_cutoff": 4.3,
            "ionic_distance_cutoff": 5.5,
            "hydrophobic_distance_cutoff": 4.8,
            "cation_pi_distance_cutoff": 5.5,
            "salt_bridge_distance_cutoff": 3.8, "salt_bridge_centroid_cutoff": 4.8,
            "sulfur_pi_distance_cutoff": 5.5, "sulfur_pi_max_perp": 2.8,
            "metal_coordination_primary_cutoff": 2.5, "metal_coordination_extended_cutoff": 2.9
        },
        "literature_default": {
            "hbond_distance_cutoff": 3.5, "hbond_angle_cutoff": 120.0,
            "halogen_distance_cutoff": 4.0, "halogen_angle_cutoff": 140.0,
            "chalcogen_distance_cutoff": 4.0, "chalcogen_theta_min": 115.0, "chalcogen_theta_max": 155.0, "chalcogen_phi_max": 50.0,
            "pnictogen_distance_cutoff": 4.0, "pnictogen_angle_cutoff": 140.0,
            "tetrel_distance_cutoff": 4.0, "tetrel_angle_cutoff": 140.0,
            "ch_pi_distance_cutoff": 4.5, "ch_pi_angle_cutoff": 90.0,
            "ch_pi_min_distance": 2.0, "ch_pi_max_distance": 4.5, "ch_pi_max_angle": 90.0, "ch_pi_max_height": 2.5,
            "pi_pi_distance_cutoff": 5.5, "pi_pi_angle_cutoff": 30.0,
            "anion_pi_distance_cutoff": 5.0,
            "n_pi_star_distance_cutoff": 3.5, "n_pi_star_angle_cutoff": 120.0,
            "dispersion_min_distance": 3.5, "dispersion_distance_cutoff": 4.5,
            "ionic_distance_cutoff": 6.0,
            "hydrophobic_distance_cutoff": 5.0,
            "cation_pi_distance_cutoff": 6.0,
            "salt_bridge_distance_cutoff": 4.0, "salt_bridge_centroid_cutoff": 5.0,
            "sulfur_pi_distance_cutoff": 6.0, "sulfur_pi_max_perp": 3.0,
            "metal_coordination_primary_cutoff": 2.6, "metal_coordination_extended_cutoff": 3.0
        },
        "exploratory": {
            "hbond_distance_cutoff": 3.9, "hbond_angle_cutoff": 110.0,
            "halogen_distance_cutoff": 4.4, "halogen_angle_cutoff": 130.0,
            "chalcogen_distance_cutoff": 4.4, "chalcogen_theta_min": 110.0, "chalcogen_theta_max": 160.0, "chalcogen_phi_max": 60.0,
            "pnictogen_distance_cutoff": 4.4, "pnictogen_angle_cutoff": 130.0,
            "tetrel_distance_cutoff": 4.3, "tetrel_angle_cutoff": 135.0,
            "ch_pi_distance_cutoff": 4.8, "ch_pi_angle_cutoff": 95.0,
            "ch_pi_min_distance": 1.8, "ch_pi_max_distance": 4.8, "ch_pi_max_angle": 95.0, "ch_pi_max_height": 2.8,
            "pi_pi_distance_cutoff": 6.2, "pi_pi_angle_cutoff": 35.0,
            "anion_pi_distance_cutoff": 5.5,
            "n_pi_star_distance_cutoff": 3.7, "n_pi_star_angle_cutoff": 115.0,
            "dispersion_min_distance": 3.3, "dispersion_distance_cutoff": 4.8,
            "ionic_distance_cutoff": 6.8,
            "hydrophobic_distance_cutoff": 5.5,
            "cation_pi_distance_cutoff": 6.5,
            "salt_bridge_distance_cutoff": 4.3, "salt_bridge_centroid_cutoff": 5.3,
            "sulfur_pi_distance_cutoff": 6.3, "sulfur_pi_max_perp": 3.2,
            "metal_coordination_primary_cutoff": 2.7, "metal_coordination_extended_cutoff": 3.2
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

        # Load scenario profiles if present
        self.scenario_profiles: dict[str, dict] = {}
        profiles_path = self.templates_dir / "scenario_profiles.yaml"
        if profiles_path.exists():
            if yaml is None:
                print("Warning: PyYAML not installed; scenario profiles disabled.")
            else:
                try:
                    with profiles_path.open("r", encoding="utf-8") as fh:
                        data = yaml.safe_load(fh) or {}
                        if isinstance(data, dict):
                            # Normalize keys to snake_case
                            norm = {}
                            for k, v in data.items():
                                key = k.strip().lower().replace(" ", "_")
                                if isinstance(v, dict):
                                    norm[key] = v
                            self.scenario_profiles = norm
                except Exception as e:
                    # Non-fatal; just log via print (avoid loguru import cycle)
                    print(f"Warning: failed loading scenario profiles: {e}")

    # Convenience accessor
    def get_scenario_profile(self, name: str) -> Optional[dict]:
        if not name:
            return None
        key = name.strip().lower().replace(" ", "_")
        return self.scenario_profiles.get(key)

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
        "hydrophobiccontact",
        "cation_pi",
        "salt_bridge"
        ,"sulfur_pi"
        ,"metal_coordination"
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
        "hydrophobiccontact": "Hydrophobic Contacts",
        "cation_pi": "Cation–π Interactions"
        ,"salt_bridge": "Salt Bridges"
        ,"sulfur_pi": "Sulfur–π Interactions"
        ,"metal_coordination": "Metal Coordination"
    }
