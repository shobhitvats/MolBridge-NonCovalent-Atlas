"""Pytest-benchmark hook script for batch processing performance.

Usage (after installing pytest-benchmark extra):
  pytest -k benchmark_batch --benchmark-only
"""
from analysis.high_performance_batch import HighPerformanceBatchProcessor
from utils.config import AppConfig
import os

TEST_PDBS = os.getenv("MOLBRIDGE_BENCHMARK_PDBS", "1CRN,4HHB").split(",")

def test_benchmark_single_structure(benchmark):
    proc = HighPerformanceBatchProcessor(use_parallel=True)
    pdb_id = TEST_PDBS[0]
    def run():
        res = proc.process_single_protein(pdb_id)
        assert res['success']
    benchmark(run)

def test_benchmark_multi_structure(benchmark):
    proc = HighPerformanceBatchProcessor(use_parallel=True)
    def run():
        for pdb_id in TEST_PDBS:
            res = proc.process_single_protein(pdb_id)
            assert res['success']
    benchmark(run)