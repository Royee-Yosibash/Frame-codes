"""Unit tests for probability distribution utilities."""

import unittest

import numpy as np
from scipy.integrate import trapezoid

from frame_codes.numerics.distributions import (
    histogram_bin_edges,
    histogram_pdf,
    marchenko_pastur_pdf,
    normalize_pdf,
    sample_from_pdf,
)


class TestNormalizePDF(unittest.TestCase):
    """Test numerical density normalization."""

    def test_normalizes_nonuniform_sample_grid(self) -> None:
        """The trapezoidal integral is one on a nonuniform grid."""
        points = np.array([0.0, 0.25, 1.0, 2.0])
        density = np.array([0.0, 1.0, 0.5, 0.0])

        normalized = normalize_pdf(points, density)

        self.assertAlmostEqual(trapezoid(normalized, points), 1.0)

    def test_rejects_zero_integral(self) -> None:
        """An all-zero density cannot be normalized."""
        with self.assertRaisesRegex(ValueError, "positive integral"):
            normalize_pdf([0, 1], [0, 0])


class TestTheoreticalDistributions(unittest.TestCase):
    """Test theoretical density support and probability mass."""

    def test_marchenko_pastur_density_has_unit_mass_and_expected_support(self) -> None:
        """The beta-below-one MP density integrates to one on its support."""
        beta = 0.25
        points = np.linspace(0, 4, 20_001)

        density = marchenko_pastur_pdf(points, beta)

        self.assertAlmostEqual(trapezoid(density, points), 1.0, places=4)
        self.assertEqual(density[0], 0.0)
        self.assertEqual(density[-1], 0.0)

    def test_marchenko_pastur_rejects_beta_above_one(self) -> None:
        """The continuous-only implementation restricts beta to at most one."""
        with self.assertRaisesRegex(ValueError, "0 < beta <= 1"):
            marchenko_pastur_pdf([0.5, 1.0], 1.2)


class TestEmpiricalDistributionUtilities(unittest.TestCase):
    """Test histogram conversion and sampling from a supplied PDF."""

    def test_histogram_pdf_has_unit_mass(self) -> None:
        """Histogram density values integrate to one over their bins."""
        samples = np.array([0.1, 0.2, 0.3, 0.8, 0.9])
        edges = histogram_bin_edges(samples, bins=4)

        density = histogram_pdf(samples, edges)

        self.assertAlmostEqual(np.sum(density * np.diff(edges)), 1.0)

    def test_sampling_returns_values_within_pdf_support(self) -> None:
        """Sampling stays within the input support and follows a flat PDF."""
        rng = np.random.default_rng(5)
        points = np.linspace(0, 1, 101)

        samples = sample_from_pdf(points, np.ones(points.shape), 5_000, rng)

        self.assertTrue(np.all((samples >= 0) & (samples <= 1)))
        self.assertAlmostEqual(float(np.mean(samples)), 0.5, delta=0.03)


if __name__ == "__main__":
    unittest.main()
