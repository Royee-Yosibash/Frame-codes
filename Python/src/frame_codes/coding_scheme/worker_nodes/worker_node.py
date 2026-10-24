"""Worker logic units addressed by a stable worker identifier."""

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, Protocol

import numpy as np
from numpy.typing import DTypeLike

from frame_codes.coding_scheme.state_utils.information import Information


class WorkerOperation(Protocol):
    """Callable worker logic receiving encoded information by name."""

    def __call__(
        self,
        *,
        information: Sequence[Information],
        calculation_type: np.dtype,
        **kwargs: Any,
    ) -> Any:
        """Process an encoded information sequence.

        Args:
            information: Coefficient-weighted information items for one code
                row.
            **kwargs: Operation-specific keyword arguments.

        Returns:
            The operation's result.
        """


@dataclass(frozen=True)
class WorkerNodeConfig:
    """Configuration for one worker node."""

    operation: WorkerOperation
    calculation_type: DTypeLike

    def __post_init__(self) -> None:
        """Normalize the worker's calculation type to a NumPy dtype.

        The dataclass is frozen, so normalization replaces the initially
        supplied dtype-like value during construction.
        """
        object.__setattr__(self, "calculation_type", np.dtype(self.calculation_type))


class WorkerNode:
    """Represent one persistent worker and its operation."""

    def __init__(
        self,
        worker_id: int,
        configuration: WorkerNodeConfig,
    ) -> None:
        """Create a worker node from its single configuration source.

        Args:
            worker_id: Identifier used to route encoded inputs and results.
            configuration: Operation and calculation dtype for this worker.
        """
        self.worker_id = worker_id
        self._configuration = configuration
        self._validate()

    @property
    def calculation_type(self) -> np.dtype:
        """Return the calculation dtype assigned to this worker.

        Returns:
            This worker's calculation dtype.
        """
        return self._configuration.calculation_type

    def _validate(self):
        """Validate the worker operation and its calculation dtype.

        Raises:
            TypeError: If the operation is not callable or the dtype is not
                numeric.
        """
        if not isinstance(self._configuration.operation, Callable):
            raise TypeError("operation must be callable")
        if not np.issubdtype(self.calculation_type, np.number):
            raise TypeError("calculation_type must be a numeric dtype")

    def execute(
        self,
        *,
        information: Sequence[Information],
        **kwargs: Any,
    ) -> Any:
        """Execute this worker's operation on its encoded information.

        Args:
            information: Coefficient-weighted information items for one code
                row.
            **kwargs: Additional arguments forwarded to the worker operation,
                such as the deterministic straggler flag.

        Returns:
            The output produced by this worker's operation.
        """
        calculation_type = self.calculation_type
        calculation_information = []
        for value in information:
            numeric_information = np.asarray(value)
            if np.issubdtype(numeric_information.dtype, np.number):
                numeric_information = numeric_information.astype(
                    calculation_type
                )[()]
            calculation_information.append(numeric_information)
        return self._configuration.operation(
            information=calculation_information,
            calculation_type=calculation_type,
            **kwargs,
        )
