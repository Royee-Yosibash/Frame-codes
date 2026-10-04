"""Block rotation code families for circulant permutation coding."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class CirculantPermutationCodeFamily(CodeFamily):
    """Block rotation code families."""

    def _validate_parameters(self, parameters: CodeParameters) -> None:
        """Require even dimensions for paired rotations.

        Args:
            parameters: Code dimensions and normalization setting.

        Raises:
            ValueError: If either dimension is odd.
        """
        super()._validate_parameters(parameters)
        if parameters.n % 2 or parameters.m % 2:
            raise ValueError("Circulant Permutation requires even code dimensions")

    @staticmethod
    @cache
    def _create_code(parameters: CodeParameters) -> NDArray[np.float64]:
        """Return the cached block rotation matrix.

        Args:
            parameters: Shared code dimensions and normalization setting.

        Returns:
            A real block rotation matrix.

        """
        n = parameters.n
        m = parameters.m
        matrix = np.zeros((n, m), dtype=np.float64)
        angle_step = 4 * np.pi / n
        for row_block in range(n // 2):
            for column_block in range(m // 2):
                angle = angle_step * row_block * column_block
                cosine = np.cos(angle)
                sine = np.sin(angle)
                row = 2 * row_block
                column = 2 * column_block
                matrix[row : row + 2, column : column + 2] = (
                    (cosine, -sine),
                    (sine, cosine),
                )
        return matrix
