"""Chebyshev polynomial code families used by OrthoMatDot."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class OrthoMatDotCodeFamily(CodeFamily):
    """Chebyshev polynomial evaluation code families."""

    @staticmethod
    @cache
    def _generate_new_code(parameters: CodeParameters) -> NDArray[np.float64]:
        """Return the cached Chebyshev evaluation matrix.

        Args:
            parameters: Shared code dimensions and normalization setting.

        Returns:
            A real node-by-polynomial matrix.
        """
        n = parameters.n
        m = parameters.m
        nodes = np.cos((2 * np.arange(n) + 1) * np.pi / (2 * n))
        matrix = np.empty((n, m), dtype=np.float64)
        matrix[:, 0] = 1
        if m > 1:
            matrix[:, 1] = nodes
        for degree in range(2, m):
            matrix[:, degree] = 2 * nodes * matrix[:, degree - 1] - matrix[:, degree - 2]
        matrix[:, 0] /= np.sqrt(2)
        return matrix
