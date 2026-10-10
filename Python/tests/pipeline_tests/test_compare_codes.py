"""End-to-end tests for the code-comparison workflow."""

import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from frame_codes.coding_scheme.codes.code_type import CodeType
from frame_codes.coding_scheme.comparison.config import ComparisonConfig
from frame_codes.coding_scheme.comparison.results import write_comparison_results
from frame_codes.coding_scheme.comparison.run import run_comparison


class TestCompareCodesPipeline(unittest.TestCase):
    """Test encoding, straggler decoding, metrics, and result output together."""

    def test_comparison_writes_results_for_active_code_families(self) -> None:
        """Both active code families produce one row per straggler count."""
        config = ComparisonConfig(
            m_data_sets=(2,),
            n_nodes=(4,),
            snr_db=(80,),
            num_trials=1,
            inner_dimension=5,
            code_types=(
                CodeType.NON_CONSECUTIVE_POWERS,
                CodeType.CIRCULANT_PERMUTATION,
            ),
        )

        results = run_comparison(config, rng=np.random.default_rng(41))

        self.assertEqual(len(results), 6)
        self.assertEqual([row.straggler_count for row in results[:3]], [0, 1, 2])
        self.assertEqual([row.straggler_count for row in results[3:]], [0, 1, 2])
        with tempfile.TemporaryDirectory() as temporary_directory:
            result_path = write_comparison_results(results, Path(temporary_directory))
            with result_path.open(newline="", encoding="utf-8") as result_file:
                rows = list(csv.DictReader(result_file))
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[0]["code_type"], CodeType.NON_CONSECUTIVE_POWERS.value)


if __name__ == "__main__":
    unittest.main()
