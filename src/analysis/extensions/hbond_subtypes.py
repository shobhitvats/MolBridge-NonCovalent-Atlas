"""Hydrogen bond subtype classification extension.

Reads already detected hydrogen bonds (interaction key 'hydrogen_bond') and classifies them into subtypes:
  - backbone_backbone (both donor and acceptor heavy atoms in backbone atoms: N, O of main chain)
  - backbone_sidechain (one in backbone, one in side chain)
  - sidechain_sidechain (neither atom in backbone)
  - ligand (either residue marked as HETATM / non-standard three-letter code outside canonical AAs)
  - water_mediated (placeholder: if either residue name is HOH / WAT)

We do NOT alter the original hydrogen bond list—this module produces a summary dict plus optional annotated list.

Output structure placed under extensions['hbond_subtypes']:
{
  'counts': {subtype: count, ...},
  'total_hbonds': int,
  'fractions': {subtype: fraction},
  'annotated': [ { original hbond fields + 'subtype': ... }, ... ]
}

Fails gracefully if hydrogen bonds are absent.
"""
from typing import Dict, Any, List

CANONICAL_AA = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
}
BACKBONE_ATOMS = {'N','O','C','CA'}  # heavy backbone atoms; using N/O for classification of donor/acceptor, but broaden set to detect context

class HydrogenBondSubtypeExtension:
    def __init__(self, config=None):
        self.config = config
    
    def _is_backbone_atom(self, atom_name: str) -> bool:
        return atom_name in BACKBONE_ATOMS
    
    def _is_canonical(self, residue_name: str) -> bool:
        return residue_name.upper() in CANONICAL_AA
    
    def classify(self, interactions: Dict[str, Any]) -> Dict[str, Any]:
        hbonds = interactions.get('hydrogen_bond') or interactions.get('hydrogen_bonds') or []
        if not hbonds:
            return {
                'counts': {},
                'total_hbonds': 0,
                'fractions': {},
                'annotated': []
            }
        counts = {
            'backbone_backbone': 0,
            'backbone_sidechain': 0,
            'sidechain_sidechain': 0,
            'ligand': 0,
            'water_mediated': 0
        }
        annotated: List[Dict[str, Any]] = []  # type: ignore
        total = 0
        for hb in hbonds:
            # Expect keys: donor_atom, acceptor_atom, residue1, residue2, chain1, chain2 (based on existing detector pattern)
            donor_atom = hb.get('donor_atom') or hb.get('donor') or ''
            acceptor_atom = hb.get('acceptor_atom') or hb.get('acceptor') or ''
            res1 = hb.get('residue1','')
            res2 = hb.get('residue2','')
            resname1 = ''.join([c for c in res1 if c.isalpha()])[:3].upper()
            resname2 = ''.join([c for c in res2 if c.isalpha()])[:3].upper()
            subtype = None
            # water mediated simplistic detection
            if resname1 in {'HOH','WAT'} or resname2 in {'HOH','WAT'}:
                subtype = 'water_mediated'
            # ligand involvement
            elif not self._is_canonical(resname1) or not self._is_canonical(resname2):
                subtype = 'ligand'
            else:
                donor_backbone = self._is_backbone_atom(donor_atom)
                acceptor_backbone = self._is_backbone_atom(acceptor_atom)
                if donor_backbone and acceptor_backbone:
                    subtype = 'backbone_backbone'
                elif donor_backbone != acceptor_backbone:
                    subtype = 'backbone_sidechain'
                else:
                    subtype = 'sidechain_sidechain'
            counts[subtype] += 1
            total += 1
            annotated_entry = dict(hb)
            annotated_entry['subtype'] = subtype
            annotated.append(annotated_entry)
        fractions = {k: (counts[k]/total if total else 0.0) for k in counts}
        return {
            'counts': counts,
            'total_hbonds': total,
            'fractions': fractions,
            'annotated': annotated
        }

# Helper function so extension loader pattern is similar to others

def run_hbond_subtype_extension(interactions: Dict[str, Any], config=None) -> Dict[str, Any]:
    """Backward compatible helper (legacy name)."""
    ext = HydrogenBondSubtypeExtension(config)
    return ext.classify(interactions)


def compute(result: Dict[str, Any], config=None) -> Dict[str, Any]:
    """Entry point matching other extension modules.

    Expects a full analysis result dict with key 'interactions'. Returns the
    hydrogen bond subtype classification structure (counts, fractions, etc.).
    """
    interactions = result.get('interactions', {}) if isinstance(result, dict) else {}
    return run_hbond_subtype_extension(interactions, config)
