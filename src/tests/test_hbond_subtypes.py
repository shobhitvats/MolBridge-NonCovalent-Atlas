"""Smoke test for hydrogen bond subtype extension.

Creates a minimal fake analysis result containing a handful of hydrogen bond
records spanning different subtypes (backbone/backbone, backbone/sidechain,
sidechain/sidechain, ligand, water) and validates classification counts and
fractions sum.
"""
from analysis.extensions.hbond_subtypes import compute as compute_hbond_subtypes


def _make_hbond(donor_atom, acceptor_atom, res1, res2):
    return {
        'donor_atom': donor_atom,
        'acceptor_atom': acceptor_atom,
        'residue1': res1,
        'residue2': res2,
        'distance': 2.8,
        'angle': 150.0
    }


def test_hbond_subtype_classification_basic():
    # Construct interactions dict analogous to result['interactions']
    interactions = {
        'hydrogen_bond': [
            # backbone-backbone: both donor/acceptor heavy atoms are backbone (N/O)
            _make_hbond('N', 'O', 'ALA12', 'VAL15'),
            # backbone-sidechain: one backbone (N) one sidechain (OG)
            _make_hbond('N', 'OG', 'GLY5', 'SER8'),
            # sidechain-sidechain: both sidechain atoms
            _make_hbond('OG', 'ND2', 'SER20', 'ASN25'),
            # ligand: non-canonical residue LIG
            _make_hbond('N', 'O1', 'LIG301', 'TYR55'),
            # water mediated: water residue
            _make_hbond('N', 'O', 'HOH401', 'GLU60'),
        ]
    }
    fake_result = {'interactions': interactions}

    out = compute_hbond_subtypes(fake_result)
    counts = out['counts']
    assert counts['backbone_backbone'] == 1
    assert counts['backbone_sidechain'] == 1
    assert counts['sidechain_sidechain'] == 1
    assert counts['ligand'] == 1
    assert counts['water_mediated'] == 1
    assert out['total_hbonds'] == 5
    # Fractions should sum (within fp tolerance) to ~1.0
    total_frac = sum(out['fractions'][k] for k in counts)
    assert abs(total_frac - 1.0) < 1e-6
