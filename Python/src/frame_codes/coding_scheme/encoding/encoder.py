"""Encoding of one information set for distributed workers."""

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from numpy.typing import ArrayLike, NDArray

from frame_codes.utils.numpy_types import as_inexact_array


@dataclass(frozen=True)
class WeightedInformation:
    """One information payload paired with its worker's code coefficient."""

    information: Any
    coefficient: float | complex


class Encoder:
    """Encoder configured with a code family and its parameters."""

    @staticmethod
    def _validate_code_to_information_compatability(information_set: Sequence, code: NDArray):
        if code.ndim != 2 or code.shape[0] == 0:
            raise ValueError("encoding_matrix must have one row per worker")
        if len(information_set) == 0 or code.shape[1] != len(information_set):
            raise ValueError("encoding_matrix columns must match the information item count")

    def encode(
        self,
        information_set: Sequence,
        encoding_matrix: ArrayLike,
    ) -> list[list[WeightedInformation]]:
        """Encode one information set using the configured code matrix.

        Args:
            information_set: Information items in the order represented by the
                code matrix columns. Items may be arbitrary objects.
            encoding_matrix: Code matrix with one row per worker and one column per
                information item.

        Returns:
            Worker-specific weighted information items.
            
        Raises:
            ValueError: If the code matrix is not two-dimensional, has no workers,
                or its columns do not match the number of information items.
            TypeError: If the code matrix does not have a numeric dtype.
        """
        
        code = as_inexact_array(encoding_matrix)
        self._validate_code_to_information_compatability(information_set=information_set, code=code)
        return [
            [
                WeightedInformation(information, coefficient.item())
                for information, coefficient in zip(information_set, code_row)
            ]
            for code_row in code
        ]
