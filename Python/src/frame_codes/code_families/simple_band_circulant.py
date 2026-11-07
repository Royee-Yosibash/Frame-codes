"""Simple band-circulant code family implementation.

The generator matrix has shape (n, m) - each row gives the encoding coefficients
assigned per 'm' items of the encoded set.

The code is constructed by a series of 'bands' each containing a unique block matrix.
Every band 'b' contains all m cyclic shifts of a binary row containing b consecutive '1's
and the rest of the row is '0's.
For example, Band 0 shifts a single '1' giving the identity matrix.
Band 1 shifts two adjacent ones, band 2 shifts three, and so on.
The last band, or the terminal band, is the only band with one row of all '1's.

For example, m = 4 and n = 9 produce two complete blocks and a final row::

    [[1, 0, 0, 0],
     [0, 1, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 1],
     [1, 1, 0, 0],
     [0, 1, 1, 0],
     [0, 0, 1, 1],
     [1, 0, 0, 1],
     [1, 1, 1, 1]]

Notable properties:
1. Every item in the encoded set appears equally often in the rows.
2. Every band adds one more '1' per row -> each encoding element participates (b_i + 1)
times in construction of the band rows (not including the terminal band).
3. From the above follows that for 'b_m' the highest complete band level each column
contains (b_m + 1)(b_m + 2) / 2 '1's, and another if the terminal band is present.

"""

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class SimpleBandCirculantCode(CodeFamily):
    """Build cyclic band rows from the requested worker count."""

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
    ) -> NDArray[np.bool_]:
        """Construct the configured cyclic band rows starting at band level 0.

        Args:
            parameters: Number of encoded sets and requested workers.

        Returns:
            A boolean generator matrix shaped as workers by encoded sets.
        """
        m = parameters.m
        has_terminal_band = parameters.n % m == 1
        regular_band_count = (
            (parameters.n - 1) // m if has_terminal_band else parameters.n // m
        )
        rows = []

        for band_index in range(regular_band_count):
            pattern = np.concatenate(
                (
                    np.ones(band_index + 1, dtype=bool),
                    np.zeros(m - band_index - 1, dtype=bool),
                )
            )
            rows.extend(np.roll(pattern, shift) for shift in range(m))

        if has_terminal_band:
            rows.append(np.ones(m, dtype=bool))
        return np.asarray(rows, dtype=bool)
