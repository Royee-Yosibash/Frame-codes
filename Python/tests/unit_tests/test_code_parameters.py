"""Unit tests for shared code parameters and performance statistics."""

import unittest

import numpy as np

from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.factory import create_code_family
from frame_codes.coding_scheme.worker_runtime.statistics import (
    CodePerformanceStatistics,
    create_dual_frame_matrix,
    sample_code_statistics,
)


class TestCodeParameters(unittest.TestCase):
    """Test shared dimensions and normalization settings."""

    def test_gamma_is_derived_from_dimensions(self) -> None:
        """The aspect ratio equals m divided by n."""
        parameters = CodeParameters(m=2, n=4)

        self.assertEqual(parameters.gamma, 0.5)

    def test_code_family_uses_the_same_parameters(self) -> None:
        """The parameters object is passed directly to code construction."""
        parameters = CodeParameters(m=2, n=4, norm_dim="column")
        family = create_code_family("LPF")

        family.generate_code(parameters)
        generator = family.get_code()

        self.assertEqual(generator.shape, (4, 2))
        np.testing.assert_allclose(np.linalg.norm(generator, axis=1), 1)

    def test_none_normalization_is_default(self) -> None:
        """No axis normalization is selected by default."""
        self.assertIsNone(CodeParameters(m=2, n=4).norm_dim)

    def test_normalization_setting_is_validated(self) -> None:
        """Unsupported normalization axes fail at parameter construction."""
        with self.assertRaisesRegex(ValueError, "norm_dim"):
            CodeParameters(m=2, n=4, norm_dim="diagonal")

    def test_rejects_invalid_dimensions(self) -> None:
        """Code parameters require positive dimensions and m no greater than n."""
        with self.assertRaisesRegex(ValueError, "m must not exceed n"):
            CodeParameters(m=5, n=4)


class TestCodePerformanceStatistics(unittest.TestCase):
    """Test dual-frame construction and erasure statistics."""

    def test_dual_frame_reconstructs_the_identity(self) -> None:
        """The dual frame is a left inverse of its generator matrix."""
        parameters = CodeParameters(m=2, n=4)
        family = create_code_family("LPF")
        family.generate_code(parameters)
        generator = family.get_code()

        dual_frame = create_dual_frame_matrix("LPF", parameters)

        np.testing.assert_allclose(dual_frame @ generator, np.eye(2), atol=1e-12)

    def test_statistics_return_one_value_per_trial(self) -> None:
        """Statistics contain one noise and condition value per erasure trial."""
        parameters = CodeParameters(m=2, n=4)

        statistics = sample_code_statistics(
            "LPF",
            parameters,
            retention_probability=1.0,
            num_trials=4,
            rng=np.random.default_rng(4),
        )

        self.assertIsInstance(statistics, CodePerformanceStatistics)
        self.assertEqual(statistics.noise_amplification.shape, (4,))
        self.assertEqual(statistics.condition_numbers.shape, (4,))
        self.assertTrue(np.all(np.isfinite(statistics.noise_amplification)))
        self.assertTrue(np.all(np.isfinite(statistics.condition_numbers)))

    def test_rejects_retention_that_cannot_span_code_dimension(self) -> None:
        """Too few retained vectors cannot recover all message dimensions."""
        parameters = CodeParameters(m=2, n=4)

        with self.assertRaisesRegex(ValueError, "between m and n"):
            sample_code_statistics("LPF", parameters, 0.25, num_trials=2)


if __name__ == "__main__":
    unittest.main()
