"""Block rotation code family for circulant permutation coding."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.code_family import CodeFamily


class CirculantPermutationCodeFamily(CodeFamily):
    """Block rotation code family."""

    def _validate_dimensions(self, n_nodes: int, message_dimension: int) -> None:
        """Require even dimensions for paired rotations.

        Args:
            n_nodes: Number of encoded outputs.
            message_dimension: Number of uncoded message elements.

        Raises:
            ValueError: If either dimension is odd.
        """
        super()._validate_dimensions(n_nodes, message_dimension)
        if n_nodes % 2 or message_dimension % 2:
            raise ValueError("Circulant Permutation requires even code dimensions")

    @staticmethod
    @cache
    def _create_code(
        n: int,
        m: int,
    ) -> NDArray[np.float64]:
        """Return the cached block rotation matrix.

        Args:
            n_nodes: Number of encoded outputs; must be even.
            message_dimension: Number of message elements; must be even.

        Returns:
            A real block rotation matrix.

        """
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
