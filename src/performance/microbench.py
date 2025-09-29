"""MIT License

Copyright (c) 2025 MolBridge Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Microbenchmark harness for interaction detectors.

Runs selected detectors over representative small/medium/large PDB structures
(or synthetic coordinate expansions) and reports timing + funnel metrics.
Intended for CI p95 regression gating.

Environment variables:
  MOLBRIDGE_BENCH_PDB_SMALL  - path to small PDB (optional)
  MOLBRIDGE_BENCH_PDB_MEDIUM - path to medium PDB (optional)
  MOLBRIDGE_BENCH_PDB_LARGE  - path to large PDB (optional)

If files not provided, attempts to synthesize dummy coordinate-only structures
using Bio.PDB internal classes (no chemistry validity guaranteed) to exercise
geometry pipelines.
"""
from __future__ import annotations
import os, time, json, statistics
try:
    import psutil  # for memory profiling
    _HAVE_PSUTIL = True
except Exception:  # pragma: no cover
    _HAVE_PSUTIL = False
from typing import List, Dict, Any, Tuple
from pathlib import Path

try:
    from Bio.PDB import PDBParser, StructureBuilder
except Exception:  # pragma: no cover
    PDBParser = None
    StructureBuilder = None

from performance.parallel_processor import HighPerformanceProcessor
from utils.config import AppConfig  # assumes existing config
from analysis.registry import DETECTOR_REGISTRY

# Detector subset for microbenchmark (extend as desired)
# Use registry keys (see analysis.registry) and allow legacy underscores -> normalized by stripping them.
BENCH_DETECTORS = [
    'hydrogenbond', 'pipi', 'ionicinteraction', 'hydrophobiccontact',
    'halogenbond', 'chpi', 'chalcogenbond', 'tetrelbond', 'pnictogenbond',
    'anionpi', 'npistar'
]

REPEAT = int(os.getenv('MOLBRIDGE_BENCH_REPEAT', '2'))


def _load_structure(path: str):
    if PDBParser is None:
        raise RuntimeError('Bio.PDB not available; cannot load PDB file')
    parser = PDBParser(QUIET=True)
    return parser.get_structure(Path(path).stem, path)


def _synthesize_structure(name: str, residues: int = 100):
    # Minimal synthetic structure: chain A, residues with two atoms placed on grid
    if StructureBuilder is None:
        raise RuntimeError('StructureBuilder not available for synthesis')
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure(name)
    sb.init_model(0)
    sb.init_chain('A')
    import numpy as np
    for i in range(residues):
        sb.init_seg('    ')
        sb.init_residue('GLY', ' ', i+1, ' ')
        # two pseudo atoms (N, O) for donor/acceptor exercise
        coord = np.array([i*1.5, 0.0, 0.0])
        for atom_name, offset in [('N', 0.0), ('O', 0.6)]:
            sb.init_atom(atom_name, coord + np.array([0.0, offset, 0.0]), 1.0, 1.0, ' ', atom_name, i*2 + (1 if atom_name=='N' else 2), 'N' if atom_name=='N' else 'O')
    return sb.get_structure()


def _get_structures() -> List[Tuple[str, Any]]:
    structs = []
    small = os.getenv('MOLBRIDGE_BENCH_PDB_SMALL')
    medium = os.getenv('MOLBRIDGE_BENCH_PDB_MEDIUM')
    large = os.getenv('MOLBRIDGE_BENCH_PDB_LARGE')
    try:
        if small and Path(small).is_file():
            structs.append(('small', _load_structure(small)))
        else:
            structs.append(('small_synth', _synthesize_structure('small_synth', residues=80)))
    except Exception as e:
        print(f"Skipping small: {e}")
    try:
        if medium and Path(medium).is_file():
            structs.append(('medium', _load_structure(medium)))
        else:
            structs.append(('medium_synth', _synthesize_structure('medium_synth', residues=200)))
    except Exception as e:
        print(f"Skipping medium: {e}")
    try:
        if large and Path(large).is_file():
            structs.append(('large', _load_structure(large)))
        else:
            structs.append(('large_synth', _synthesize_structure('large_synth', residues=400)))
    except Exception as e:
        print(f"Skipping large: {e}")
    return structs


def run_microbench():
    config = AppConfig()
    # Seed for reproducibility (structure synthesis ordering)
    seed_env = os.getenv('MOLBRIDGE_BENCH_SEED')
    if seed_env:
        try:
            import random, numpy as np  # type: ignore
            random.seed(int(seed_env))
            np.random.seed(int(seed_env))
        except Exception:
            pass

    # Map registry names to detector classes (instantiate later via processor path)
    detector_specs = []  # list[(class, method_name, key)]
    for name in BENCH_DETECTORS:
        key = name.replace('_', '')
        entry = DETECTOR_REGISTRY.get(key)
        if entry:
            cls, method_name = entry
            detector_specs.append((cls, method_name, key))
    if not detector_specs:
        raise RuntimeError('No detectors available for microbenchmark (check BENCH_DETECTORS keys)')
    structures = _get_structures()
    processor = HighPerformanceProcessor(config=config)
    results: Dict[str, Any] = {'runs': []}
    # Build list of classes only for processor API compatibility
    detector_classes_only = [cls for cls, _, _ in detector_specs]
    for s_name, s_obj in structures:
        for r in range(REPEAT):
            t0 = time.time()
            start_mem = psutil.Process(os.getpid()).memory_info().rss if _HAVE_PSUTIL else None
            _ = processor.process_interactions_parallel(s_obj, detector_classes_only, config=config)
            dur = time.time() - t0
            funnel = getattr(processor, '_detector_instrumentation', {})
            end_mem = psutil.Process(os.getpid()).memory_info().rss if _HAVE_PSUTIL else None
            mem_delta = (end_mem - start_mem) if (_HAVE_PSUTIL and start_mem is not None and end_mem is not None) else None
            results['runs'].append({
                'structure': s_name,
                'repeat': r,
                'duration_sec': dur,
                'funnel': funnel,
                'rss_delta_bytes': mem_delta
            })
    # Aggregate simple stats
    agg: Dict[str, Any] = {}
    for dname in BENCH_DETECTORS:
        key = dname.replace('_', '')
        detector_funnels = [r['funnel'][key] for r in results['runs'] if key in r.get('funnel', {}) and isinstance(r['funnel'][key], dict)]
        times = [r['duration_sec'] for r in results['runs'] if key in r.get('funnel', {})]
        if times and detector_funnels:
            times_sorted = sorted(times)
            p95 = times_sorted[int(len(times_sorted)*0.95)-1] if len(times_sorted) > 1 else times_sorted[-1]
            # Acceptance ratio distribution
            acc_ratios = [f.get('acceptance_ratio', 0.0) for f in detector_funnels]
            mean_acc = statistics.mean(acc_ratios) if acc_ratios else 0.0
            # Simple sanity flags (could evolve to env-configured bounds)
            warnings = []
            if mean_acc > 0.9:
                warnings.append('high_accept_ratio')
            if mean_acc < 0.0001 and any(f.get('candidate_pairs', 0) > 100 for f in detector_funnels):
                warnings.append('near_zero_accept_ratio')
            # Aggregate funnel totals
            total_raw = sum(f.get('raw_pairs', 0) for f in detector_funnels)
            total_candidate = sum(f.get('candidate_pairs', 0) for f in detector_funnels)
            total_accepted = sum(f.get('accepted_pairs', 0) for f in detector_funnels)
            agg[key] = {
                'mean': statistics.mean(times),
                'p95': p95,
                'n': len(times),
                'mean_acceptance_ratio': mean_acc,
                'total_raw_pairs': total_raw,
                'total_candidate_pairs': total_candidate,
                'total_accepted_pairs': total_accepted,
                'warnings': warnings
            }
    results['aggregate'] = agg
    # Optional OpenTelemetry stub (no-op if env not set) for future integration
    if os.getenv('MOLBRIDGE_OTEL_EXPORT', '0') in {'1','true','True'}:
        try:  # pragma: no cover
            # Placeholder: user can plug real OTLP exporter; we just emit a simple line.
            print(json.dumps({'event': 'otel_stub', 'detectors': list(agg.keys())}))
        except Exception:
            pass
    print(json.dumps(results, indent=2, default=lambda o: '<nonserializable>'))
    return results

if __name__ == '__main__':  # pragma: no cover
    run_microbench()
