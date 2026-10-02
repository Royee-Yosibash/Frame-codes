"""Probability density functions and empirical distribution utilities."""

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid

from frame_codes.utils.numpy_types import as_inexact_array


def normalize_pdf(sample_points: ArrayLike, density_values: ArrayLike) -> NDArray[np.float64]:
    """Normalize sampled density values to have unit integral.

    Args:
        sample_points: At least two strictly increasing real sample points.
        density_values: Non-negative real density values at each sample point.

    Returns:
        A normalized density array with the same shape as `density_values`.

    Raises:
        ValueError: If the arrays have incompatible shapes, invalid sample
            points, negative density values, or a non-positive integral.
        TypeError: If either input is complex or has a non-numeric dtype.
    """
    points = _as_real_vector(sample_points, "sample_points")
    density = _as_real_vector(density_values, "density_values")
    if points.size < 2 or points.shape != density.shape:
        raise ValueError("sample_points and density_values must have the same length >= 2")
    if not np.all(np.isfinite(points)) or np.any(np.diff(points) <= 0):
        raise ValueError("sample_points must be finite and strictly increasing")
    if not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ValueError("density_values must be finite and non-negative")

    area = float(trapezoid(density, points))
    if area <= 0:
        raise ValueError("density_values must have a positive integral")
    return density / area


def marchenko_pastur_pdf(sample_points: ArrayLike, beta: float) -> NDArray[np.float64]:
    """Evaluate the continuous Marchenko-Pastur probability density.

    Args:
        sample_points: Real points at which to evaluate the density.
        beta: Aspect ratio satisfying `0 < beta <= 1`.

    Returns:
        The density values, zero outside the Marchenko-Pastur support.

    Raises:
        ValueError: If `beta` is outside `(0, 1]` or sample points are not
            finite real values.
        TypeError: If sample points are complex or non-numeric.
    """
    points = _as_real_vector(sample_points, "sample_points")
    _validate_beta(beta)
    _validate_finite_points(points, "sample_points")

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
    """Evaluate the continuous MANOVA probability density.

    The density integrates to one minus `manova_atom_mass(beta, gamma)`.

    Args:
        sample_points: Real points at which to evaluate the density.
        beta: Sub-frame aspect ratio satisfying `0 < gamma <= beta <= 1`.
        gamma: Frame aspect ratio satisfying `0 < gamma <= beta`.

    Returns:
        The continuous density values, zero outside the MANOVA support.

    Raises:
        ValueError: If the aspect ratios are outside their supported range or
            sample points are not finite real values.
        TypeError: If sample points are complex or non-numeric.
    """
    points = _as_real_vector(sample_points, "sample_points")
    _validate_manova_parameters(beta, gamma)
    left = np.sqrt((1 - gamma) / beta)
    right = np.sqrt(1 - gamma / beta)
    lower = (left - right) ** 2
    upper = (left + right) ** 2
    _validate_finite_points(points, "sample_points")
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
    """Return the MANOVA point-mass probability at `1 / gamma`.

    Args:
        beta: Sub-frame aspect ratio satisfying `0 < gamma <= beta <= 1`.
        gamma: Frame aspect ratio satisfying `0 < gamma <= beta`.

    Returns:
        The probability mass at `1 / gamma`, or zero when there is no atom.

    Raises:
        ValueError: If the aspect ratios are outside their supported range.
    """
    _validate_manova_parameters(beta, gamma)
    excess = gamma * (1 + 1 / beta) - 1
    if excess <= 0:
        return 0.0
    return float(excess / min(gamma, gamma / beta))


def histogram_bin_edges(eigenvalues: ArrayLike, bins: int) -> NDArray[np.float64]:
    """Calculate equal-width histogram edges covering the eigenvalues.

    Args:
        eigenvalues: Non-empty finite real eigenvalue samples.
        bins: Positive integer number of bins.

    Returns:
        The histogram bin edges.

    Raises:
        ValueError: If samples are empty or non-finite, or if `bins` is not a
            positive integer.
        TypeError: If samples are complex or non-numeric.
    """
    values = _as_real_vector(eigenvalues, "eigenvalues")
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("eigenvalues must be non-empty and finite")
    if isinstance(bins, bool) or not isinstance(bins, (int, np.integer)) or bins <= 0:
        raise ValueError("bins must be a positive integer")
    return np.histogram_bin_edges(values, bins=int(bins))


def histogram_pdf(eigenvalues: ArrayLike, bin_edges: ArrayLike) -> NDArray[np.float64]:
    """Calculate a normalized empirical density for histogram bins.

    Args:
        eigenvalues: Non-empty finite real eigenvalue samples.
        bin_edges: At least two finite, strictly increasing real bin edges.

    Returns:
        One density value per bin, normalized so bin mass sums to one.

    Raises:
        ValueError: If the samples or bin edges are invalid.
        TypeError: If samples or bin edges are complex or non-numeric.
    """
    values = _as_real_vector(eigenvalues, "eigenvalues")
    edges = _as_real_vector(bin_edges, "bin_edges")
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("eigenvalues must be non-empty and finite")
    if edges.size < 2 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError("bin_edges must be finite and strictly increasing")
    if np.any(values < edges[0]) or np.any(values > edges[-1]):
        raise ValueError("bin_edges must cover all eigenvalue samples")

    counts, _ = np.histogram(values, bins=edges)
    widths = np.diff(edges)
    return counts / (values.size * widths)


def sample_from_pdf(
    sample_points: ArrayLike,
    density_values: ArrayLike,
    sample_count: int,
    rng: np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Draw samples from a PDF represented by points and density values.

    Args:
        sample_points: At least two strictly increasing real sample points.
        density_values: Non-negative density values at each sample point.
        sample_count: Non-negative number of samples to draw.
        rng: Optional NumPy random generator. A default generator is created
            when omitted.

    Returns:
        A one-dimensional array of sampled values.

    Raises:
        ValueError: If the PDF inputs or sample count are invalid.
        TypeError: If PDF inputs are complex or non-numeric, or if
            `sample_count` is not an integer.
    """
    if isinstance(sample_count, bool) or not isinstance(sample_count, (int, np.integer)):
        raise TypeError("sample_count must be an integer")
    if sample_count < 0:
        raise ValueError("sample_count must be a non-negative integer")
    points = _as_real_vector(sample_points, "sample_points")
    density = normalize_pdf(sample_points, density_values)
    if points.size < 2:
        raise ValueError("sample_points and density_values must have the same length >= 2")

    segment_mass = (density[:-1] + density[1:]) * np.diff(points) / 2
    cumulative_mass = np.concatenate(([0.0], np.cumsum(segment_mass)))
    cumulative_mass /= cumulative_mass[-1]
    generator = np.random.default_rng() if rng is None else rng
    uniform_samples = generator.random(sample_count)
    return np.interp(uniform_samples, cumulative_mass, points)


def _as_real_vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    """Convert a real numeric input to a one-dimensional float array.

    Args:
        values: NumPy-compatible real numeric values.
        name: Argument name used in validation messages.

    Returns:
        A one-dimensional `float64` array.

    Raises:
        ValueError: If `values` is not one-dimensional.
        TypeError: If values are complex or have a non-numeric dtype.
    """
    array = as_inexact_array(values)
    if np.iscomplexobj(array):
        raise TypeError(f"{name} must contain real values")
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return array


def _validate_finite_points(points: NDArray[np.float64], name: str) -> None:
    """Require finite real sample points.

    Args:
        points: One-dimensional real sample points.
        name: Argument name used in validation messages.

    Raises:
        ValueError: If any sample point is not finite.
    """
    if not np.all(np.isfinite(points)):
        raise ValueError(f"{name} must contain only finite values")


def _validate_beta(beta: float) -> None:
    """Validate the supported Marchenko-Pastur aspect ratio.

    Args:
        beta: Aspect ratio to validate.

    Raises:
        ValueError: If `beta` is not finite or outside `(0, 1]`.
    """
    if not np.isfinite(beta) or not 0 < beta <= 1:
        raise ValueError("beta must satisfy 0 < beta <= 1")


def _validate_manova_parameters(beta: float, gamma: float) -> None:
    """Validate MANOVA aspect ratios used by the frame workflows.

    Args:
        beta: Sub-frame aspect ratio.
        gamma: Frame aspect ratio.

    Raises:
        ValueError: If either ratio is non-finite or they do not satisfy
            `0 < gamma <= beta <= 1`.
    """
    if not np.isfinite(beta) or not np.isfinite(gamma) or not 0 < gamma <= beta <= 1:
        raise ValueError("MANOVA parameters must satisfy 0 < gamma <= beta <= 1")


def _manova_raw_continuous_mass(
    beta: float,
    gamma: float,
    lower: float,
    upper: float,
) -> float:
    """Calculate the integral of the unscaled continuous MANOVA density.

    Args:
        beta: Sub-frame aspect ratio.
        gamma: Frame aspect ratio.
        lower: Lower support endpoint.
        upper: Upper support endpoint.

    Returns:
        The continuous mass before accounting for any point mass.
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
