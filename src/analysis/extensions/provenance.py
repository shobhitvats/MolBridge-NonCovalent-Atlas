"""Provenance helper: captures configuration and runtime context for a structure analysis."""
from typing import Dict, Any
from datetime import datetime

def compute(result: Dict[str, Any], config) -> Dict[str, Any]:
    meta = result.get('metadata', {})
    return {
        'generated': datetime.utcnow().isoformat() + 'Z',
        'app_version': config.version,
        'assembly_mode': config.default_assembly,
        'include_ligands': config.include_ligands,
        'exclude_waters': config.exclude_waters,
        'parameter_snapshot': {
            'hbond_distance_cutoff': config.interactions.hbond_distance_cutoff,
            'hbond_angle_cutoff': config.interactions.hbond_angle_cutoff,
            'pi_pi_distance_cutoff': config.interactions.pi_pi_distance_cutoff,
            'pi_pi_angle_cutoff': config.interactions.pi_pi_angle_cutoff,
            'ionic_distance_cutoff': config.interactions.ionic_distance_cutoff,
            'hydrophobic_distance_cutoff': config.interactions.hydrophobic_distance_cutoff,
        },
        'timing': {
            'analysis_time': meta.get('analysis_time'),
            'timestamp': meta.get('timestamp')
        }
    }
