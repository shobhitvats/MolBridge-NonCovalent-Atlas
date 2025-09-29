"""Extension modules for advanced protein interaction analytics.

Each module exposes a `compute(result, config)` function returning a dictionary
that can be inserted under `analysis_results[pdb_id]['extensions'][<name>]`.
Lazy, non-destructive, and optional based on AppConfig feature toggles.
"""

from .residue_profiles import compute as compute_residue_profiles  # noqa: F401
from .interfaces import compute as compute_interface_analysis  # noqa: F401
from .outliers import compute as compute_outliers  # noqa: F401
from .provenance import compute as compute_provenance  # noqa: F401
from .motifs import compute as compute_motifs  # noqa: F401
from .secondary_structure import compute as compute_secondary_structure  # noqa: F401
from .sasa_bsa import compute as compute_sasa_bsa  # noqa: F401
from .geometry_quality import compute as compute_geometry_quality  # noqa: F401
from .disulfides import compute as compute_disulfides  # noqa: F401
from .pockets import compute as compute_pockets  # noqa: F401
from .conservation_stub import compute as compute_conservation  # noqa: F401
from .pi_pi_refinement import compute as compute_pi_pi_refinement  # noqa: F401
from .hbond_subtypes import compute as compute_hbond_subtypes  # noqa: F401

__all__ = [
    'compute_residue_profiles',
    'compute_interface_analysis',
    'compute_outliers',
    'compute_provenance',
    'compute_motifs',
    'compute_secondary_structure',
    'compute_sasa_bsa',
    'compute_geometry_quality',
    'compute_disulfides',
    'compute_pockets',
    'compute_conservation',
    'compute_pi_pi_refinement'
    ,'compute_hbond_subtypes'
]
