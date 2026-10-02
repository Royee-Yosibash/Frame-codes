"""Chebyshev polynomial code family used by OrthoMatDot."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.code_family import CodeFamily


class OrthoMatDotCodeFamily(CodeFamily):
    """Chebyshev polynomial evaluation code family."""

    @staticmethod
    @cache
    def _create_code(
        n_nodes: int,
        message_dimension: int,
    ) -> NDArray[np.float64]:
        """Return the cached Chebyshev evaluation matrix.

        Args:
            n_nodes: Number of Chebyshev evaluation nodes.
            message_dimension: Number of polynomial terms.

        Returns:
            A real node-by-polynomial matrix.
        """
        nodes = np.cos((2 * np.arange(n_nodes) + 1) * np.pi / (2 * n_nodes))
        matrix = np.empty((n_nodes, message_dimension), dtype=np.float64)
        matrix[:, 0] = 1
        if message_dimension > 1:
            matrix[:, 1] = nodes
        for degree in range(2, message_dimension):
            matrix[:, degree] = 2 * nodes * matrix[:, degree - 1] - matrix[:, degree - 2]
        matrix[:, 0] /= np.sqrt(2)
        return matrix
