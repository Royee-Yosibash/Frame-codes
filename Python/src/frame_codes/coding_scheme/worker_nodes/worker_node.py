"""Worker logic units addressed by a stable worker identifier."""

from collections.abc import Callable, Sequence
from typing import Any, Protocol

from frame_codes.coding_scheme.state_utils.encoded_information import (
    WeightedInformation,
)


class WorkerOperation(Protocol):
    """Callable worker logic receiving encoded information by name."""

    def __call__(
        self,
        *,
        information: Sequence[WeightedInformation],
        **kwargs: Any,
    ) -> Any:
        """Process an encoded information sequence.

        Args:
            information: Weighted inputs for this worker operation.
            **kwargs: Operation-specific keyword arguments.

        Returns:
            The operation's result.
        """


class WorkerNode:
    """Represent one persistent worker and its operation."""

    def __init__(self, worker_id: int, operation: WorkerOperation) -> None:
        """Create a worker node with a stable ID and operation.

        Args:
            worker_id: Identifier used to route encoded inputs and results.
            operation: Function that performs this worker's computation.
        """
        self.worker_id = worker_id
        self._operation = operation
        self._validate()

    def _validate(self):
        if not isinstance(self._operation, Callable):
            raise TypeError('Operation must be a callable.')

    def execute(
        self,
        *,
        information: Sequence[WeightedInformation],
        **kwargs: Any,
    ) -> Any:
        """Execute this worker's operation on its encoded information.

        Args:
            information: Weighted inputs for this worker operation.
            **kwargs: Additional arguments forwarded to the worker operation,
                such as the deterministic straggler flag.

        Returns:
            The output produced by this worker's operation.
        """
        return self._operation(information=information, **kwargs)
