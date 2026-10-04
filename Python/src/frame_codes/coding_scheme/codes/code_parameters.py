"""Shared code dimensions and normalization settings."""

from dataclasses import dataclass

from frame_codes.numerics.matrix_shapes import MatrixAxis


@dataclass(frozen=True)
class CodeParameters:
    """Frame dimensions and normalization shared by code and experiment APIs."""

    m: int
    n: int
    norm_dim: MatrixAxis | None = None

    def __post_init__(self) -> None:
        """Validate the normalization setting."""
        if self.norm_dim not in {None, "row", "column"}:
            raise ValueError("norm_dim must be None, 'row', or 'column'")

    @property
    def gamma(self) -> float:
        """Return the frame aspect ratio derived from `m / n`."""
        return self.m / self.n
