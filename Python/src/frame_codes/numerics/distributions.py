"""Probability densities and empirical distribution utilities."""

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid


def normalize_pdf(sample_points: ArrayLike, density_values: ArrayLike) -> NDArray[np.float64]:
    """Normalize sampled density values to have unit integral.

    Args:
        sample_points: Sample coordinates.
        density_values: Density at each coordinate.

    Returns:
        A normalized density array.

    Raises:
        ValueError: If the sampled density has a non-positive integral.
    """
    points = np.asarray(sample_points, dtype=np.float64)
    density = np.asarray(density_values, dtype=np.float64)
    area = trapezoid(density, points)
    if area <= 0:
        raise ValueError("density must have a positive integral")
    return density / area


def marchenko_pastur_pdf(sample_points: ArrayLike, beta: float) -> NDArray[np.float64]:
    """Evaluate the continuous Marchenko-Pastur density for `0 < beta <= 1`.

    Args:
        sample_points: Points at which to evaluate the density.
        beta: Aspect ratio in `(0, 1]`.

    Returns:
        Density values, zero outside the Marchenko-Pastur support.

    Raises:
        ValueError: If `beta` is outside `(0, 1]`.
    """
    if not 0 < beta <= 1:
        raise ValueError("beta must satisfy 0 < beta <= 1")
    points = np.asarray(sample_points, dtype=np.float64)
    lower = (1 - np.sqrt(beta)) ** 2
    upper = (1 + np.sqrt(beta)) ** 2
    support = (points >= lower) & (points <= upper) & (points > 0)
    density = np.zeros(points.shape, dtype=np.float64)
    density[support] = np.sqrt(
        (points[support] - lower) * (upper - points[support])
    ) / (2 * np.pi * beta * points[support])
    return density


def manova_pdf(
    sample_points: ArrayLike,
    beta: float,
    gamma: float,
) -> NDArray[np.float64]:
    """Evaluate the continuous MANOVA density.

    Args:
        sample_points: Points at which to evaluate the density.
        beta: Sub-frame aspect ratio.
        gamma: Frame aspect ratio.

    Returns:
        MANOVA density values at the points.

    Raises:
        ValueError: If the aspect ratios do not satisfy
            `0 < gamma <= beta <= 1`.
    """
    if not 0 < gamma <= beta <= 1:
        raise ValueError("MANOVA parameters must satisfy 0 < gamma <= beta <= 1")
    points = np.asarray(sample_points, dtype=np.float64)
    density = np.zeros(points.shape, dtype=np.float64)

    left = np.sqrt((1 - gamma) / beta)
    right = np.sqrt(1 - gamma / beta)
    r_minus = (left - right) ** 2
    r_plus = (left + right) ** 2
    support = (points >= r_minus) & (points <= r_plus)
    above_fraction = beta * np.sqrt(
        (points[support] - r_minus) * (r_plus - points[support])
    )
    below_fraction = (2 * np.pi * points[support]) * (1 - gamma * points[support])
    density[support] = above_fraction / below_fraction
    density = normalize_pdf(points, density)

    atom_mass = (gamma * (1 + 1 / beta) - 1) / min(gamma, gamma / beta)
    if atom_mass > 0:
        atom_index = np.abs(points - 1 / gamma).argmin()
        density = density * (1 - atom_mass)
        density.flat[atom_index] += atom_mass
    return density


def histogram_bin_edges(eigenvalues: ArrayLike, bins: int) -> NDArray[np.float64]:
    """Calculate equal-width histogram edges for eigenvalue samples.

    Args:
        eigenvalues: Eigenvalue samples.
        bins: Number of bins.

    Returns:
        Histogram bin edges.
    """
    values = np.asarray(eigenvalues, dtype=np.float64)
    return np.histogram_bin_edges(values, bins=bins)


def histogram_pdf(eigenvalues: ArrayLike, bin_edges: ArrayLike) -> NDArray[np.float64]:
    """Calculate empirical density values for histogram bins.

    Args:
        eigenvalues: Eigenvalue samples.
        bin_edges: Histogram bin edges.

    Returns:
        Density per bin, normalized by sample count and bin width.
    """
    values = np.asarray(eigenvalues, dtype=np.float64)
    edges = np.asarray(bin_edges, dtype=np.float64)
    counts, _ = np.histogram(values, bins=edges)
    return counts / (values.size * np.diff(edges))


def sample_from_pdf(
    sample_points: ArrayLike,
    density_values: ArrayLike,
    sample_count: int,
    rng: np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Draw samples from a PDF represented by points and density values.

    Args:
        sample_points: Sample coordinates.
        density_values: Density at each coordinate.
        sample_count: Number of samples to draw.
        rng: Optional random generator.

    Returns:
        A one-dimensional array of samples.
    """
    points = np.asarray(sample_points, dtype=np.float64)
    density = normalize_pdf(points, density_values)
    segment_mass = (density[:-1] + density[1:]) * np.diff(points) / 2
    cumulative_mass = np.concatenate(([0.0], np.cumsum(segment_mass)))
    generator = np.random.default_rng() if rng is None else rng
    return np.interp(generator.random(sample_count), cumulative_mass, points)
