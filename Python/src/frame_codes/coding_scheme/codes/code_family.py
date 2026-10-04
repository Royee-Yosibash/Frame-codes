"""Abstract base classes for deterministic and randomized code families."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.numerics.linear_algebra import normalize_code


class CodeFamily(ABC):
    """Base class for node-by-message code generators."""

    @staticmethod
    def _validate_parameters(parameters: CodeParameters) -> None:
        """Validate dimensions common to all code families.

        Args:
            parameters: Code dimensions and normalization setting.
        """
        return

    def create_code(self, parameters: CodeParameters) -> NDArray:
        """Create a code matrix using shared code parameters.

        Args:
            parameters: Code dimensions and normalization setting.

        Returns:
            A generator matrix with shape `(parameters.n, parameters.m)`.
        """
        self._validate_parameters(parameters)
        code = self._create_code(parameters)
        if parameters.norm_dim is None:
            return code.copy()
        frame = code.T
        return normalize_code(frame, parameters.norm_dim).T

    @abstractmethod
    def _create_code(self, parameters: CodeParameters) -> NDArray:
        """Construct the families-specific matrix.

        Args:
            parameters: Shared dimensions and normalization setting.

        Returns:
            A code matrix with shape `(parameters.n, parameters.m)`.
        """


class RandomCodeFamily(CodeFamily, ABC):
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
