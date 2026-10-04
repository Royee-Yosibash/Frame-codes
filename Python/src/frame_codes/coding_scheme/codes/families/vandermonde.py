"""Random row-selected root-of-unity Vandermonde code families."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import RandomCodeFamily


class VandermondeCodeFamily(RandomCodeFamily):
    """Root-of-unity Vandermonde families with random row selection."""

    @staticmethod
    @cache
    def _create_base_code(n: int) -> NDArray[np.complex128]:
        """Create and cache the full root-of-unity Vandermonde matrix.

        Args:
            n: Number of available roots of unity and powers.

        Returns:
            A node-by-power Vandermonde matrix.
        """
        roots = np.exp(-2j * np.pi * np.arange(n) / n)
        powers = np.arange(n)
        matrix = roots[:, np.newaxis] ** powers[np.newaxis, :]
        return matrix

    def _create_code(
        self,
        n: int,
        m: int,
    ) -> NDArray[np.complex128]:
        """Select nonconsecutive powers from a cached Vandermonde matrix.

        Args:
            n: Number of available roots of unity and powers.
            m: Number of powers to select.

        Returns:
            An `n`-by-`m` generator matrix with selected power columns.
        """
        power_indices = np.sort(
            self._rng.choice(n, size=m, replace=False)
        )
        return self._create_base_code(n)[:, power_indices]
