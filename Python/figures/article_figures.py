"""Create theoretical and empirical article distribution figures."""

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.code_type import CodeType
from frame_codes.coding_scheme.codes.factory import create_code_family
from frame_codes.numerics.distributions import (
    histogram_pdf,
    manova_pdf,
    marchenko_pastur_pdf,
)
from frame_codes.numerics.linear_algebra import gram_matrix_eigenvalues
from numpy.typing import NDArray


@dataclass(frozen=True)
class ArticleFigureConfig:
    """Dimensions and trial counts for article distribution figures."""

    n: int = 400
    gamma: float = 0.5
    beta: float = 0.8
    coded_trials: int = 1_000
    wishart_trials: int = 10_000
    pdf_points: int = 100_001

    def __post_init__(self) -> None:
        """Validate article figure dimensions and aspect ratios."""
        if isinstance(self.n, bool) or not isinstance(self.n, (int, np.integer)):
            raise TypeError("n must be an integer")
        if self.n <= 0:
            raise ValueError("n must be positive")
        if not 0 < self.gamma <= self.beta <= 1:
            raise ValueError("aspect ratios must satisfy 0 < gamma <= beta <= 1")
        if self.coded_trials <= 0 or self.wishart_trials <= 0:
            raise ValueError("trial counts must be positive")
        if self.pdf_points < 2:
            raise ValueError("pdf_points must be at least two")


def create_article_figures(
    config: ArticleFigureConfig | None = None,
    output_directory: Path | None = None,
    rng: np.random.Generator | None = None,
) -> list[Path]:
    """Create theoretical and empirical distribution figures.

    Args:
        config: Figure dimensions and simulation counts.
        output_directory: Figure destination. Defaults to
            `Python/outputs/article_figures`.
        rng: Optional generator controlling stochastic code and matrix draws.

    Returns:
        Paths to the four generated PNG figures.

    Raises:
        TypeError: If `rng` is not a NumPy random generator.
        ValueError: If configuration values are invalid.
    """
    config = ArticleFigureConfig() if config is None else config
    if rng is not None and not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")
    generator = np.random.default_rng() if rng is None else rng
    destination = (
        Path(__file__).resolve().parents[1] / "outputs" / "article_figures"
        if output_directory is None
        else output_directory
    )
    destination.mkdir(parents=True, exist_ok=True)

    m = int(np.floor(config.n * config.gamma + 0.5))
    wishart_m = int(np.floor(config.n * config.beta + 0.5))
    sample_points = np.linspace(0, 10, config.pdf_points)
    manova_density = manova_pdf(sample_points, config.beta, config.gamma)
    mp_density = marchenko_pastur_pdf(sample_points, config.beta)
    manova_path = destination / "manova_theoretical.png"
    _save_distribution_figure(
        sample_points,
        manova_density,
        manova_path,
        "MANOVA probability density",
    )

    coded_eigenvalues = _collect_coded_frame_eigenvalues(config, m, generator)
    coded_path = destination / "manova_coded_frame.png"
    _save_distribution_figure(
        sample_points,
        manova_density,
        coded_path,
        "Coded-frame eigenvalue distribution",
        eigenvalues=coded_eigenvalues,
    )

    mp_path = destination / "marchenko_pastur_theoretical.png"
    _save_distribution_figure(
        sample_points,
        mp_density,
        mp_path,
        "Marchenko-Pastur probability density",
    )

    wishart_eigenvalues = _collect_wishart_eigenvalues(
        config.n,
        wishart_m,
        config.wishart_trials,
        generator,
    )
    wishart_path = destination / "marchenko_pastur_wishart.png"
    _save_distribution_figure(
        sample_points,
        mp_density,
        wishart_path,
        "Wishart eigenvalue distribution",
        eigenvalues=wishart_eigenvalues,
    )
    return [manova_path, coded_path, mp_path, wishart_path]


def _collect_coded_frame_eigenvalues(
    config: ArticleFigureConfig,
    m: int,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Collect eigenvalues from random surviving coded-frame columns.

    Args:
        config: Article frame dimensions and aspect ratios.
        m: Rounded frame dimension derived from gamma and n.
        rng: Random generator for code construction and node selection.

    Returns:
        A flat array of Gram-matrix eigenvalues from all coded trials.
    """
    frame = create_code_family(
        CodeType.OMITTED_VANDERMONDE,
        rng,
    ).generate_code(CodeParameters(m=m, n=config.n, norm_dim="column")).T
    retained_count = int(np.floor(m / config.beta + 0.5))
    eigenvalues = []
    for _ in range(config.coded_trials):
        indices = np.sort(rng.choice(config.n, size=retained_count, replace=False))
        eigenvalues.append(gram_matrix_eigenvalues(frame[:, indices]))
    return np.concatenate(eigenvalues)


def _collect_wishart_eigenvalues(
    n: int,
    m: int,
    num_trials: int,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Collect eigenvalues from independent normalized Gaussian matrices.

    Args:
        n: Number of columns in each Wishart sample.
        m: Number of rows in each Wishart sample.
        num_trials: Number of Gaussian matrices to sample.
        rng: Random generator for matrix entries.

    Returns:
        A flat array of Gram-matrix eigenvalues from all trials.
    """
    code_parameters = CodeParameters(m=m, n=n)
    eigenvalues = []
    for _ in range(num_trials):
        matrix = create_code_family(CodeType.WISHART, rng).generate_code(code_parameters).T
        eigenvalues.append(gram_matrix_eigenvalues(matrix))
    return np.concatenate(eigenvalues)


def _save_distribution_figure(
    sample_points: NDArray[np.float64],
    density: NDArray[np.float64],
    output_path: Path,
    title: str,
    eigenvalues: NDArray[np.float64] | None = None,
) -> None:
    """Plot a theoretical density and optional empirical eigenvalue histogram.

    Args:
        sample_points: Points used to evaluate the theoretical density.
        density: Theoretical continuous density values.
        output_path: Destination image path.
        title: Figure title.
        eigenvalues: Optional empirical eigenvalue samples.
    """
    figure, axes = plt.subplots()
    axes.plot(sample_points, density, color="tab:red", label="Theoretical PDF")
    if eigenvalues is not None:
        _add_empirical_histogram(axes, eigenvalues)
    axes.set_xlabel("Eigenvalue")
    axes.set_ylabel("Probability density")
    axes.set_title(title)
    axes.grid(True)
    axes.legend()
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)


def _add_empirical_histogram(
    axes: plt.Axes,
    eigenvalues: NDArray[np.float64],
) -> None:
    """Add an empirical density histogram to a plot.

    Args:
        axes: Matplotlib axes to draw on.
        eigenvalues: Empirical eigenvalue samples.
    """
    maximum = max(2 * float(np.max(np.abs(eigenvalues))), 1.0)
    edges = np.linspace(0, maximum, 401)
    density = histogram_pdf(eigenvalues, edges)
    centers = (edges[:-1] + edges[1:]) / 2
    axes.plot(
        centers,
        density,
        color="black",
        marker=".",
        linestyle="",
        label="Empirical density",
    )
