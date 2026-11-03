"""Simple band-circulant code family with discrete worker counts."""

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class SimpleBandCirculantCode(CodeFamily):
    """Build identity and cyclic band rows from the requested worker count."""

    @staticmethod
    def _validate_parameters(parameters: CodeParameters) -> None:
        """Validate dimensions supported by the banded row construction.

        Args:
            parameters: Information-set and worker counts.

        Raises:
            ValueError: If fewer than two encoded sets are requested or if the
                worker count does not correspond to a supported band prefix.
        """
        if parameters.m < 2:
            raise ValueError("SimpleBandCirculantCode requires m >= 2")

        regular_band_worker_counts = {
            parameters.m * band_count for band_count in range(1, parameters.m)
        }
        terminal_band_worker_counts = {
            parameters.m * band_count + 1
            for band_count in range(1, parameters.m)
        }
        supported_worker_counts = regular_band_worker_counts | terminal_band_worker_counts
        if parameters.n not in supported_worker_counts:
            raise ValueError(
                "n must match a supported SimpleBandCirculantCode band count"
            )

    def _generate_new_code(
        self,
        parameters: CodeParameters,
    ) -> NDArray[np.float64]:
        """Construct the configured identity and cyclic band rows.

        Args:
            parameters: Number of encoded sets and requested workers.

        Returns:
            A real generator matrix shaped as workers by encoded sets.
        """
        m = parameters.m
        has_terminal_band = parameters.n % m == 1
        regular_band_count = (
            (parameters.n - 1) // m if has_terminal_band else parameters.n // m
        )
        rows = [row.copy() for row in np.eye(m, dtype=np.float64)]

        for band_index in range(1, regular_band_count):
            pattern = np.concatenate(
                (
                    np.ones(band_index + 1, dtype=np.float64),
                    np.zeros(m - band_index - 1, dtype=np.float64),
                )
            )
            rows.extend(np.roll(pattern, shift) for shift in range(m))

        if has_terminal_band:
            rows.append(np.ones(m, dtype=np.float64))
        return np.asarray(rows, dtype=np.float64)
