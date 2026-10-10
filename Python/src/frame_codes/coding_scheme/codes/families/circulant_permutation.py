"""Block rotation code families for circulant permutation coding."""

from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class CirculantPermutationCodeFamily(CodeFamily):
    """Block rotation code families."""

    @staticmethod
    def _validate_code(code: NDArray, parameters: CodeParameters) -> None:
        """Validate the expanded code shape for paired rows and columns.

        Args:
            code: Generated code matrix.
            parameters: Information-set and worker counts.

        Raises:
            ValueError: If the matrix does not have two rows per worker and
                two columns per information set.
        """
        if code.shape != (2 * parameters.n, 2 * parameters.m):
            raise ValueError(
                "Code shape must have two rows per worker and two columns per information set"
            )

    def number_of_workers_required(self) -> int:
        """Return the number of workers represented by the code.

        Returns:
            Half the number of code rows, since each worker owns a row pair.
        """
        return self.get_code().shape[0] // 2

    def _get_coefficient_per_set(
        self,
        worker_id: int,
        set_index: int,
    ) -> NDArray[np.float64]:
        """Return the coefficient pair for one set and worker.

        Args:
            worker_id: Zero-based worker identifier.
            set_index: Column index of the encoded set.

        Returns:
            Coefficients from the worker's adjacent code rows.
        """
        code = self.get_code()
        first_row = 2 * worker_id
        return code[first_row : first_row + 2, set_index]

    @staticmethod
    @cache
    def _generate_new_code(parameters: CodeParameters) -> NDArray[np.float64]:
        """Return the block rotation matrix with two rows per worker.

        Args:
            parameters: Number of information sets and workers.

        Returns:
            A real block rotation matrix with shape `(2 * n, 2 * m)`.

        """
        row_count = 2 * parameters.n
        column_count = 2 * parameters.m
        matrix = np.zeros((row_count, column_count), dtype=np.float64)
        angle_step = 4 * np.pi / row_count
        for row_block in range(parameters.n):
            for column_block in range(parameters.m):
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
