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
        n: int,
        m: int,
    ) -> NDArray[np.float64]:
        """Return the cached Chebyshev evaluation matrix.

        Args:
            n: Number of Chebyshev evaluation nodes.
            m: Number of polynomial terms.

        Returns:
            A real node-by-polynomial matrix.
        """
        nodes = np.cos((2 * np.arange(n) + 1) * np.pi / (2 * n))
        matrix = np.empty((n, m), dtype=np.float64)
        matrix[:, 0] = 1
        if m > 1:
            matrix[:, 1] = nodes
        for degree in range(2, m):
            matrix[:, degree] = 2 * nodes * matrix[:, degree - 1] - matrix[:, degree - 2]
        matrix[:, 0] /= np.sqrt(2)
        return matrix
