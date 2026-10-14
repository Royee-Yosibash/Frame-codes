"""Dispatch encoded information to workers and model their responses."""

from collections.abc import Callable, Mapping, Sequence
from typing import Any

from frame_codes.coding_scheme.state_utils.encoded_information import EncodedInformation
from frame_codes.coding_scheme.worker_nodes.worker_node import (
    WorkerNode,
    WorkerOperation,
)
from frame_codes.coding_scheme.worker_runtime.runtime import WorkerResult

ResponsePolicy = Callable[[Sequence[int]], Sequence[int]]
OutputAlteration = Callable[[int, Any], Any]


class WorkerDispatcher:
    """Coordinate worker execution and simulate which results are returned."""

    def __init__(self):
        """Store the persistent worker pool and response models.
        """
        self.workers: Mapping[int, WorkerNode] = {}

    def setup(
        self,
        number_of_workers_required: int,
        parameters_per_worker: list[WorkerOperation] | WorkerOperation,
    ) -> None:
        """Set up the requested number of workers from per-worker parameters.

        Args:
            number_of_workers_required: Number of worker nodes to configure.
            parameters_per_worker: Worker operations, ordered by worker ID.

        Raises:
            ValueError: If there is not one operation per requested worker.
        """

        if not isinstance(parameters_per_worker, list):
            parameters_per_worker = [
                parameters_per_worker
            ] * number_of_workers_required

        if len(parameters_per_worker) != number_of_workers_required:
            raise ValueError("one worker operation per worker")

        self.workers = {
            worker_id: WorkerNode(worker_id, operation)
            for worker_id, operation in enumerate(parameters_per_worker)
        }

    def dispatch(
        self,
        encoded_information_sets: Sequence[Sequence[EncodedInformation]],
    ) -> list[WorkerResult]:
        """Run workers and omit any worker that returns a failed answer.

        Args:
            encoded_information_sets: Encoded information grouped by source
                set, with each group containing worker assignments.

        Returns:
            Results from workers whose operations returned answers.
        """
        raise NotImplementedError


    def _run_workers(
        self,
    ) -> list[WorkerResult]:
        """Execute the configured operation for every encoded worker row.
        
        Returns:
            One result containing ordered row outputs per worker.
        """
        raise NotImplementedError
