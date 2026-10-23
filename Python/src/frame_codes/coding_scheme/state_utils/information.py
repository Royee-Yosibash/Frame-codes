"""Data structures for encoded worker inputs."""

from dataclasses import dataclass
from typing import Protocol, Self


class Information(Protocol):
    """Payload that supports multiplication by an encoding coefficient."""

    def __mul__(self, coefficient: complex, /) -> Self:
        """Apply a code coefficient to this information payload.

        Args:
            coefficient: Numeric encoding coefficient.

        Returns:
            The coefficient-weighted information payload.
        """

    def __add__(self, other: Self, /) -> Self:
        """Add another compatible information payload.

        Args:
            other: Information payload to add.

        Returns:
            The sum of both payloads.
        """

    def __radd__(self, other: int, /) -> Self:
        """Support adding this payload with Python's integer sum start.

        Args:
            other: Integer initial value used by `sum`.

        Returns:
            The sum with the initial value.
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
