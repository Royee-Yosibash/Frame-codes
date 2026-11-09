"""Compare worker precision in one configured matrix multiplication run."""

from pathlib import Path
import unittest

from hydra import compose, initialize_config_dir
import numpy as np

from frame_codes.coding_scheme.hydra_pipeline.runner import run_experiment
from frame_codes.comparison.metrics import relative_frobenius_error


class TestHydraMatrixMultiplicationPrecision(unittest.TestCase):
    """Verify the YAML run compares float64 and float16 worker results."""

    def test_precision_comparison_uses_one_yaml_run(self) -> None:
        """Compare both worker precisions against the exact product."""
        config_directory = Path(__file__).resolve().parents[2] / "experiments"
        with initialize_config_dir(version_base=None, config_dir=str(config_directory)):
            config = compose(config_name="matrix_multiplication_precision")

        workspace = run_experiment(config)

        exact_result = np.eye(3)
        np.testing.assert_allclose(
            np.asarray(config.init.lhs) @ np.asarray(config.init.rhs),
            exact_result,
            atol=1e-15,
        )
        results_by_precision = {
            trial["worker_calculation_type"]: trial["result"]
            for trial in workspace["sweep_results"]
        }
        float64_relative_error = relative_frobenius_error(
            results_by_precision["float64"],
            exact_result,
        )
        float16_relative_error = relative_frobenius_error(
            results_by_precision["float16"],
            exact_result,
        )
        precision_difference = relative_frobenius_error(
            results_by_precision["float16"],
            results_by_precision["float64"],
        )

        self.assertEqual(len(workspace["sweep_results"]), 2)
        self.assertEqual(set(results_by_precision), {"float64", "float16"})
        self.assertLess(float64_relative_error, 1e-12)
        self.assertGreater(float16_relative_error, 1e-5)
        self.assertGreater(float16_relative_error, float64_relative_error)
        self.assertGreater(precision_difference, 1e-5)
