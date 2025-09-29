"""Lightweight settings layer wrapping environment variables with validation.

Does not replace AppConfig; augments it. Use get_settings() where env-driven behavior
is needed (e.g., toggling JSON logging, cache hot size override already handled separately).
"""
from __future__ import annotations
from functools import lru_cache

# Pydantic v2 moved BaseSettings into separate package `pydantic-settings`.
# Provide a backward-compatible import so environments with pydantic v1 still work.
try:  # Prefer Pydantic v2
    from pydantic_settings import BaseSettings  # type: ignore
    from pydantic import Field, ConfigDict
    _IS_PYDANTIC_V2 = True
except ImportError:  # Fallback for pydantic v1.x
    from pydantic import BaseSettings, Field  # type: ignore
    _IS_PYDANTIC_V2 = False

class Settings(BaseSettings):
    """Environment settings with backward compatibility for Pydantic v1 and v2.

    Uses ConfigDict when on v2; retains inner Config class for v1.
    """
    if _IS_PYDANTIC_V2:  # type: ignore
        model_config = ConfigDict(case_sensitive=False, extra='ignore')  # type: ignore

    log_level: str = Field("INFO", validation_alias=None, json_schema_extra={'env': 'MOLBRIDGE_LOG_LEVEL'})
    json_logging: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_JSON_LOGS'})
    enable_metrics_endpoint: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_METRICS_ENDPOINT'})
    cli_default_preset: str = Field("literature_default", json_schema_extra={'env': 'MOLBRIDGE_CLI_PRESET'})
    # Feature flags
    enable_unified_processor: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_ENABLE_UNIFIED_PROCESSOR'})
    enable_normalization: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_ENABLE_NORMALIZATION'})
    enable_vector_geom: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_ENABLE_VECTOR_GEOM'})
    enable_provenance: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_ENABLE_PROVENANCE'})
    performance_mode: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_PERFORMANCE_MODE'})
    enable_columnar: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_ENABLE_COLUMNAR'})
    direct_columnar_serialization: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_DIRECT_COLUMNAR_JSON'})
    adaptive_parallel_threshold_ms: int = Field(5, json_schema_extra={'env': 'MOLBRIDGE_ADAPTIVE_PAR_MS'})
    partial_normalization: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_PARTIAL_NORMALIZATION'})
    rolling_acceptance_window: int = Field(10, json_schema_extra={'env': 'MOLBRIDGE_ROLLING_ACCEPT_WINDOW'})
    force_columnar_in_performance: bool = Field(True, json_schema_extra={'env': 'MOLBRIDGE_FORCE_COLUMNAR_PERF'})
    metrics_history_size: int = Field(500, json_schema_extra={'env': 'MOLBRIDGE_METRICS_HISTORY_SIZE'})
    # Logging verbosity gating (demote detector-level info logs in performance mode unless explicitly allowed)
    verbose_detector_logs: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_VERBOSE_DETECTOR_LOGS'})
    # Detector/object pooling & spatial index reuse
    enable_detector_pool: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_ENABLE_DETECTOR_POOL'})
    detector_pool_max_size: int = Field(64, json_schema_extra={'env': 'MOLBRIDGE_DETECTOR_POOL_MAX'})
    # API compression & payload tuning
    api_enable_gzip: bool = Field(False, json_schema_extra={'env': 'MOLBRIDGE_API_GZIP'})
    api_gzip_min_size: int = Field(1024, json_schema_extra={'env': 'MOLBRIDGE_API_GZIP_MIN'})
    # Precision control for serialization (fallback when user does not specify query param)
    default_float_precision: int = Field(3, json_schema_extra={'env': 'MOLBRIDGE_DEFAULT_PRECISION'})

    if not _IS_PYDANTIC_V2:  # Keep legacy Config for v1
        class Config:
            case_sensitive = False

@lru_cache(maxsize=1)
def get_settings() -> Settings:
    return Settings()  # type: ignore
