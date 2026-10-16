"""Data structures for encoded worker inputs."""

from dataclasses import dataclass
from typing import Protocol


class Information(Protocol):
    """Payload that supports multiplication by an encoding coefficient."""

    def __mul__(self, coefficient: complex, /):
        """Apply a code coefficient to this information payload.

        Args:
            coefficient: Numeric encoding coefficient.

        Returns:
            The coefficient-weighted information payload.
        """


@dataclass(frozen=True)
class WeightedInformation:
    """One information payload paired with its code coefficient."""

    information: Information
    coefficient: float | complex


@dataclass(frozen=True)
class EncodedInformation:
    """Encoded inputs for one logical worker."""

    information_id: int
    information: list[list[WeightedInformation]]
