"""Lightweight test harness for extension compute functions.

Run manually (or integrate with pytest later):
python -m analysis.extensions.test_extension_harness
"""
from analysis.extensions import (
    compute_residue_profiles,
    compute_interface_analysis,
    compute_outliers,
    compute_provenance,
    compute_motifs,
)
from utils.config import load_config

# Minimal synthetic analysis_results-like structure for a single PDB
MOCK_RESULT = {
    'pdb_id': 'XXXX',
    'interactions': {
        'hydrogen_bond': [
            {'residue1': 'ALA10', 'chain1': 'A', 'residue2': 'LYS45', 'chain2': 'A', 'distance': 2.9, 'angle': 155, 'strength': 'strong'},
            {'residue1': 'ASN12', 'chain1': 'A', 'residue2': 'GLU80', 'chain2': 'B', 'distance': 3.2, 'angle': 140, 'strength': 'moderate'},
        ],
        'pi_pi': [
            {'residue1': 'PHE30', 'chain1': 'A', 'residue2': 'TYR55', 'chain2': 'B', 'distance': 5.1, 'angle': 20, 'strength': 'moderate'}
        ],
        'ionic': [
            {'residue1': 'LYS45', 'chain1': 'A', 'residue2': 'ASP90', 'chain2': 'B', 'distance': 4.8, 'strength': 'strong'}
        ],
        'hydrophobic': [
            {'residue1': 'LEU25', 'chain1': 'A', 'residue2': 'VAL60', 'chain2': 'B', 'distance': 4.6, 'strength': 'weak'}
        ]
    },
    'metadata': {'analysis_time': 0.42, 'timestamp': '2025-09-27T00:00:00Z'}
}


def main():
    cfg = load_config()
    print("Running extension harness on synthetic data...\n")
    print("Residue Profiles:")
    print(compute_residue_profiles(MOCK_RESULT, cfg))

    print("\nInterface Analysis:")
    print(compute_interface_analysis(MOCK_RESULT, cfg))

    print("\nOutliers:")
    print(compute_outliers(MOCK_RESULT, cfg))

    print("\nProvenance:")
    print(compute_provenance(MOCK_RESULT, cfg))

    # Motifs will likely return empty due to missing coordinate centroids
    print("\nMotifs (experimental):")
    print(compute_motifs(MOCK_RESULT, cfg))


if __name__ == "__main__":
    main()
