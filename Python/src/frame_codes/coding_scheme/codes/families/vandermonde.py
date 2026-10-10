"""Random row-selected root-of-unity Vandermonde code families."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import RandomCodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


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

    def _generate_new_code(
        self,
        parameters: CodeParameters,
    ) -> NDArray[np.complex128]:
        """Select nonconsecutive powers from a cached Vandermonde matrix.

        Args:
            parameters: Shared dimensions and normalization setting.

        Returns:
            An `n`-by-`m` generator matrix with selected power columns.
        """
        power_indices = np.sort(
            self._rng.choice(parameters.n, size=parameters.m, replace=False)
        )
        return self._create_base_code(parameters.n)[:, power_indices]
