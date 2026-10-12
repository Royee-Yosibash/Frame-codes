"""Worker logic units addressed by a stable worker identifier."""

from collections.abc import Callable
from typing import Any

WorkerOperation = Callable[..., Any]


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