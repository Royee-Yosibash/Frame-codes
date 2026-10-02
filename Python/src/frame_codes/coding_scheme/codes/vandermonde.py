"""Random row-selected root-of-unity Vandermonde code family."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.code_family import RandomCodeFamily


class VandermondeCodeFamily(RandomCodeFamily):
    """Root-of-unity Vandermonde family with random row selection."""

    @staticmethod
    @cache
    def _create_base_code(
        n_nodes: int,
        message_dimension: int,
    ) -> NDArray[np.complex128]:
        """Create and cache the full root-of-unity Vandermonde matrix.

        Args:
            n_nodes: Number of available roots of unity.
            message_dimension: Number of polynomial terms.

        Returns:
            A node-by-power Vandermonde matrix.
        """
        roots = np.exp(-2j * np.pi * np.arange(n_nodes) / n_nodes)
        powers = np.arange(message_dimension)
        matrix = roots[:, np.newaxis] ** powers[np.newaxis, :]
        return matrix

    def _create_code(
        self,
        n_nodes: int,
        message_dimension: int,
    ) -> NDArray[np.complex128]:
        """Select message-dimension rows from a cached Vandermonde matrix.

        Args:
            n_nodes: Number of available roots of unity.
            message_dimension: Number of selected nodes and polynomial terms.

        Returns:
            A complex generator matrix with selected rows.
        """
        row_indices = np.sort(
            self._rng.choice(n_nodes, size=message_dimension, replace=False)
        )
        return self._create_base_code(n_nodes, message_dimension)[row_indices]
