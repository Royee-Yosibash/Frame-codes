"""Data structures for encoded worker inputs."""

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class WeightedInformation:
    """One information payload paired with its code coefficient."""

    information: Any
    coefficient: float | complex


@dataclass(frozen=True)
class EncodedInformation:
    """Encoded inputs for one logical worker."""

    information_id: int
    information: list[list[WeightedInformation]]