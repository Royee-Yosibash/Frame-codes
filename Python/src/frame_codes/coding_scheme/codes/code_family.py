"""Abstract base classes for deterministic and randomized code families."""

from abc import ABC, abstractmethod
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.state_utils.encoded_information import WeightedInformation
from frame_codes.numerics.linear_algebra import normalize_code


class CodeFamily(ABC):
    """Base class for node-by-message code generators."""
    
    _current_code = None 

    @staticmethod
    def _validate_parameters(parameters: CodeParameters) -> None:
        """Validate dimensions common to all code families.

        Args:
            parameters: Code dimensions and normalization setting.
        """
        return

    def reset_code(self) -> None:
        """Reset the currently saved code"""
        self._current_code = None

    @staticmethod
    def _validate_code(code: NDArray, parameters: CodeParameters) -> None:
        if (code.shape[1] != parameters.m) or (code.shape[0] != parameters.n):
            raise ValueError("Code shape inconsistent with parameters required")

    def _update_code(self, code: NDArray, parameters: CodeParameters) -> None:
        self._validate_code(code, parameters)
        self._current_code = code

    def get_code(self) -> NDArray:
        if self._current_code is None:
            raise AttributeError("No code has been generated.")
        return self._current_code.copy()

    def number_of_encoded_sets(self):
        return self.get_code().shape[1]
        
    def number_of_workers_required(self):
        return self.get_code().shape[0]

    def _get_coefficient_per_set(
        self,
        worker_id: int,
        set_index: int,
    ) -> NDArray:
        """Return this worker's coefficients for one encoded set.

        Args:
            worker_id: Zero-based worker identifier.
            set_index: Column index of the encoded set.

        Returns:
            One coefficient for each code row assigned to the worker.
        """
        return self.get_code()[worker_id : worker_id + 1, set_index]

    def encode_for_worker(
        self,
        worker_id: int,
        information_set: Sequence,
    ) -> list[list[WeightedInformation]]:
        """Encode information items for one worker using this family's code.

        Args:
            worker_id: Zero-based worker identifier.
            information_set: Items represented by the generated code columns.

        Returns:
            The encoded code rows assigned to the requested worker.

        Raises:
            AttributeError: If no code has been generated.
            TypeError: If `worker_id` is not an integer.
            ValueError: If the worker ID, row layout, or information count is
                incompatible with the generated code.
        """
        worker_count = self.number_of_workers_required()
        worker_id = int(worker_id)
        if worker_id < 0 or worker_id >= worker_count:
            raise ValueError("worker_id is outside the generated code's worker range")
        if len(information_set) != self.number_of_encoded_sets():
            raise ValueError("information item count must match the generated code")
        coefficients_by_set = [
            self._get_coefficient_per_set(worker_id, set_index)
            for set_index in range(len(information_set))
        ]

        information = [[] for _ in range(len(coefficients_by_set[0]))]
        for item, coefficients in zip(information_set, coefficients_by_set):
            for row, coefficient in zip(information, coefficients):
                row.append(WeightedInformation(item, coefficient.item()))
        return information

    @abstractmethod
    def _generate_new_code(self, parameters: CodeParameters) -> NDArray:
        """Construct the family-specific code generator matrix.

        Args:
            parameters: CodeParameters.

        Returns:
            A code matrix with shape `(parameters.n, parameters.m)`.
        """

    def generate_code(self, parameters: CodeParameters) -> None:
        """Generate and store a code matrix.

        Args:
            parameters: Code dimensions and normalization setting.
        """
        self._validate_parameters(parameters)
        code = self._generate_new_code(parameters)
        if parameters.norm_dim is not None:
            code = normalize_code(code.T, parameters.norm_dim).T
        self._update_code(code, parameters)


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
