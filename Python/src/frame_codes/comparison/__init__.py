"""Code-comparison configuration, results, and performance statistics."""

from frame_codes.comparison.config import ComparisonConfig
from frame_codes.comparison.results import ComparisonResult, write_comparison_results
from frame_codes.comparison.statistics import (
    CodePerformanceStatistics,
    create_dual_frame_matrix,
    sample_code_statistics,
    sample_condition_numbers,
    select_surviving_worker_ids,
)

__all__ = [
    "CodePerformanceStatistics",
    "ComparisonConfig",
    "ComparisonResult",
    "create_dual_frame_matrix",
    "sample_code_statistics",
    "sample_condition_numbers",
    "select_surviving_worker_ids",
    "write_comparison_results",
]
