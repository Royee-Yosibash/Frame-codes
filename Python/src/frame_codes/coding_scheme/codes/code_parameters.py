"""Shared code dimensions and normalization settings."""

from dataclasses import dataclass

import numpy as np

from frame_codes.numerics.matrix_shapes import MatrixAxis


@dataclass(frozen=True)
class CodeParameters:
    """Frame dimensions and normalization shared by code and experiment APIs."""

    m: int
    n: int
    norm_dim: MatrixAxis | None = None

    def __post_init__(self) -> None:
        """Validate the normalization setting."""
        if self.n <= 0 or self.m <= 0:
            raise ValueError("Dimensions must be positive")
        if self.m > self.n:
            raise ValueError("m must not exceed n")

        dimensions = (self.n, self.m)
        if any(
            isinstance(value, bool) or not isinstance(value, (int, np.integer))
            for value in dimensions
        ):
            raise TypeError("code dimensions must be integers")
        
        if self.norm_dim not in {None, "row", "column"}:
            raise ValueError("norm_dim must be None, 'row', or 'column'")

    @property
    def gamma(self) -> float:
        """Return the frame aspect ratio derived from `m / n`."""
        return self.m / self.n
