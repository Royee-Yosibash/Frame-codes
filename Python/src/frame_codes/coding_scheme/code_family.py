"""Abstract base classes for deterministic and randomized code families."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray


class CodeFamily(ABC):
    """Base class for node-by-message code generators."""

    @staticmethod
    def _validate_dimensions(n_nodes: int, message_dimension: int) -> None:
        """Validate dimensions common to all code families.

        Args:
            n_nodes: Number of encoded outputs.
            message_dimension: Number of uncoded message elements.

        Raises:
            TypeError: If either dimension is not an integer.
            ValueError: If dimensions are not positive or the message
                dimension exceeds the node count.
        """
        dimensions = (n_nodes, message_dimension)
        if any(
            isinstance(value, bool) or not isinstance(value, (int, np.integer))
            for value in dimensions
        ):
            raise TypeError("code dimensions must be integers")
        if n_nodes <= 0 or message_dimension <= 0:
            raise ValueError("code dimensions must be positive")
        if message_dimension > n_nodes:
            raise ValueError("message_dimension must not exceed n_nodes")

    def create_code(self, n_nodes: int, message_dimension: int) -> NDArray:
        """Create a code matrix.

        Args:
            n_nodes: Number of encoded outputs.
            message_dimension: Number of uncoded message elements.

        Returns:
            A code matrix with shape `(n_nodes, message_dimension)`.

        Raises:
            TypeError: If either dimension is not an integer.
            ValueError: If dimensions are not positive, message dimension
                exceeds node count, or the selected family has more specific
                dimension requirements.
        """
        self._validate_dimensions(n_nodes, message_dimension)
        return self._create_code(n_nodes, message_dimension).copy()

    @abstractmethod
    def _create_code(self, n: int, m: int) -> NDArray:
        """Construct the family-specific matrix.

        Args:
            n: Number of encoded outputs.
            m: Number of uncoded message elements.

        Returns:
            A code matrix with shape `(n, m)`.
        """


class RandomCodeFamily(CodeFamily):
    """Base class for families that use randomized construction."""

    def __init__(self, rng: np.random.Generator | None = None) -> None:
        """Initialize a random code family.

        Args:
            rng: Optional NumPy random generator. A fresh unseeded generator
                is created when omitted.

        Raises:
            TypeError: If `rng` is not a NumPy random generator or `None`.
        """
        if rng is not None and not isinstance(rng, np.random.Generator):
            raise TypeError("rng must be a numpy.random.Generator")
        self._rng = np.random.default_rng() if rng is None else rng
