"""Shared code dimensions and normalization settings."""

from dataclasses import dataclass
from typing import Literal

NormDim = Literal["None", "Row", "Column"]


@dataclass(frozen=True)
class CodeParameters:
    """Frame dimensions and normalization shared by code and experiment APIs."""

    m: int
    n: int
    norm_dim: NormDim = "None"

    def __post_init__(self) -> None:
        """Validate the normalization setting."""
        if self.norm_dim not in {"None", "Row", "Column"}:
            raise ValueError("norm_dim must be 'None', 'Row', or 'Column'")

    @property
    def gamma(self) -> float:
        """Return the frame aspect ratio derived from `m / n`."""
        return self.m / self.n
