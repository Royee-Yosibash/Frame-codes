"""Per-trial comparison metrics and experiment-level aggregates."""

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike

from frame_codes.utils.numpy_types import as_inexact_array


@dataclass(init=False)
class TrialMetrics:
    """Scalar error metrics and metadata for one comparison trial."""

    trial_index: int
    straggler_worker_ids: tuple[int, ...]
    responding_worker_ids: tuple[int, ...]
    worker_calculation_type: str
    mse: float
    relative_frobenius_error: float

    def __init__(
        self,
        expected_result: ArrayLike,
        actual_result: ArrayLike,
        trial_index: int,
        straggler_worker_ids: Sequence[int],
        responding_worker_ids: Sequence[int],
        worker_calculation_type: str,
    ) -> None:
        """Calculate trial metrics without retaining the matrix results.

        Args:
            expected_result: Reference matrix product.
            actual_result: Decoded worker result for this trial.
            trial_index: Zero-based index of this trial within its experiment.
            straggler_worker_ids: Workers omitted in this trial.
            responding_worker_ids: Workers that returned results.
            worker_calculation_type: Calculation dtype used by the workers.

        Raises:
            ValueError: If the result shapes differ or the reference has zero
                squared magnitude.
            TypeError: If either result has an unsupported numeric dtype.
        """
        expected = as_inexact_array(expected_result)
        actual = as_inexact_array(actual_result)
        if expected.shape != actual.shape:
            raise ValueError("actual and expected must have matching shapes")
        expected_energy = np.sum(np.abs(expected) ** 2)
        if expected_energy == 0:
            raise ValueError("expected result must have nonzero Frobenius norm")
        error_energy = np.sum(np.abs(actual - expected) ** 2)

        self.trial_index = trial_index
        self.straggler_worker_ids = tuple(straggler_worker_ids)
        self.responding_worker_ids = tuple(responding_worker_ids)
        self.worker_calculation_type = worker_calculation_type
        self.mse = float(error_energy / expected_energy)
        self.relative_frobenius_error = float(np.sqrt(error_energy / expected_energy))


@dataclass
class ExperimentResult:
    """Trial metrics and aggregates for one configured code-family experiment."""

    name: str
    code_type: str
    n_nodes: int
    m_data_sets: int
    worker_calculation_type: str
    trials: Sequence[TrialMetrics]
    mean_mse: float = field(init=False)
    mean_relative_frobenius_error: float = field(init=False)
    min_relative_frobenius_error: float = field(init=False)
    max_relative_frobenius_error: float = field(init=False)

    def __post_init__(self) -> None:
        """Normalize trials and calculate aggregate error metrics.

        Raises:
            ValueError: If no trial results were supplied.
        """
        trial_results = tuple(self.trials)
        if not trial_results:
            raise ValueError("an experiment requires at least one trial result")
        relative_errors = np.asarray(
            [trial.relative_frobenius_error for trial in trial_results],
            dtype=np.float64,
        )
        self.trials = trial_results
        self.mean_mse = float(np.mean([trial.mse for trial in trial_results]))
        self.mean_relative_frobenius_error = float(np.mean(relative_errors))
        self.min_relative_frobenius_error = float(np.min(relative_errors))
        self.max_relative_frobenius_error = float(np.max(relative_errors))
