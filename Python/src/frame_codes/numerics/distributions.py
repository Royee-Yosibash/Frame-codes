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

    The returned density has total mass `1 - manova_atom_mass(beta, gamma)`.

    Args:
        sample_points: Points at which to evaluate the density.
        beta: Sub-frame aspect ratio.
        gamma: Frame aspect ratio.

    Returns:
        Continuous density values with total mass `1 - manova_atom_mass`,
        zero outside the MANOVA support.

    Raises:
        ValueError: If the aspect ratios do not satisfy
            `0 < gamma <= beta <= 1`.
    """
    if not 0 < gamma <= beta <= 1:
        raise ValueError("MANOVA parameters must satisfy 0 < gamma <= beta <= 1")
    points = np.asarray(sample_points, dtype=np.float64)
    left = np.sqrt((1 - gamma) / beta)
    right = np.sqrt(1 - gamma / beta)
    lower = (left - right) ** 2
    upper = (left + right) ** 2
    continuous_mass = 1 - manova_atom_mass(beta, gamma)
    density = np.zeros(points.shape, dtype=np.float64)
    if continuous_mass == 0:
        return density
    raw_mass = _manova_raw_continuous_mass(beta, gamma, lower, upper)
    support = (
        (points >= lower)
        & (points <= upper)
        & (points > 0)
        & (1 - gamma * points > 0)
    )
    density[support] = (
        beta
        * np.sqrt((points[support] - lower) * (upper - points[support]))
        / (2 * np.pi * points[support] * (1 - gamma * points[support]))
        * (continuous_mass / raw_mass)
    )
    return density


def manova_atom_mass(beta: float, gamma: float) -> float:
    """Calculate the MANOVA discrete mass at `1 / gamma`.

    Args:
        beta: Sub-frame aspect ratio.
        gamma: Frame aspect ratio.

    Returns:
        Discrete probability mass, or zero if there is no atom.

    Raises:
        ValueError: If aspect ratios do not satisfy `0 < gamma <= beta <= 1`.
    """
    if not 0 < gamma <= beta <= 1:
        raise ValueError("MANOVA parameters must satisfy 0 < gamma <= beta <= 1")
    excess = gamma * (1 + 1 / beta) - 1
    if excess <= 0:
        return 0.0
    return float(excess / min(gamma, gamma / beta))


def _manova_raw_continuous_mass(
    beta: float,
    gamma: float,
    lower: float,
    upper: float,
) -> float:
    """Calculate the unscaled integral of the continuous MANOVA density.

    Args:
        beta: Sub-frame aspect ratio.
        gamma: Frame aspect ratio.
        lower: Lower support endpoint.
        upper: Upper support endpoint.

    Returns:
        The total mass of the unscaled continuous density.
    """
    center = (lower + upper) / 2
    half_width = (upper - lower) / 2
    reciprocal_gamma = 1 / gamma
    root_at_zero = np.sqrt(lower * upper)
    root_at_reciprocal = np.sqrt(
        max((reciprocal_gamma - center) ** 2 - half_width**2, 0)
    )
    return float(
        beta
        / 2
        * (reciprocal_gamma - root_at_zero - root_at_reciprocal)
    )


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
    cumulative_mass /= cumulative_mass[-1]
    generator = np.random.default_rng() if rng is None else rng
    return np.interp(generator.random(sample_count), cumulative_mass, points)
