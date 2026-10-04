"""Abstract base classes for deterministic and randomized code families."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.parameters import CodeParameters
from frame_codes.numerics.linear_algebra import normalize_code


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

    def create_code(self, parameters: CodeParameters) -> NDArray:
        """Create a code matrix using shared code parameters.

        Args:
            parameters: Code dimensions and normalization setting.

        Returns:
            A generator matrix with shape `(parameters.n, parameters.m)`.

        Raises:
            TypeError: If either dimension is not an integer.
            ValueError: If dimensions are not positive, message dimension
                exceeds node count, or the selected families has more specific
                dimension requirements.
        """
        self._validate_dimensions(parameters.n, parameters.m)
        code = self._create_code(parameters.n, parameters.m)
        if parameters.norm_dim is None:
            return code.copy()
        frame = code.T
        return normalize_code(frame, parameters.norm_dim).T

    @abstractmethod
    def _create_code(self, n: int, m: int) -> NDArray:
        """Construct the families-specific matrix.

        Args:
            n: Number of encoded outputs.
            m: Number of uncoded message elements.

        Returns:
            A code matrix with shape `(n, m)`.
        """


class RandomCodeFamily(CodeFamily):
    """Base class for families that use randomized construction."""

    def __init__(self, rng: np.random.Generator | None = None) -> None:
        """Initialize a random code families.

        Args:
            rng: Optional NumPy random generator. A fresh unseeded generator
                is created when omitted.

        Raises:
            TypeError: If `rng` is not a NumPy random generator or `None`.
        """
        if rng is not None and not isinstance(rng, np.random.Generator):
            raise TypeError("rng must be a numpy.random.Generator")
        self._rng = np.random.default_rng() if rng is None else rng
