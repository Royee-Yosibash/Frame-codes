"""Dispatch encoded information to workers and model their responses."""

from collections.abc import Mapping, Sequence

from frame_codes.coding_scheme.state_utils.information import (
    EncodedInformation,
    Information,
)
from frame_codes.coding_scheme.worker_nodes.worker_node import (
    WorkerNode,
    WorkerNodeConfig,
)
from frame_codes.coding_scheme.worker_nodes.worker_result import WorkerResult


class WorkerDispatcher:
    """Coordinate worker execution and simulate which results are returned."""

    def __init__(self):
        """Store the persistent worker pool and response models.
        """
        self.workers: Mapping[int, WorkerNode] = {}

    def setup(
        self,
        number_of_workers_required: int,
        parameters_per_worker: Sequence[WorkerNodeConfig],
    ) -> None:
        """Set up the requested number of workers from per-worker parameters.

        Args:
            number_of_workers_required: Number of worker nodes to configure.
            parameters_per_worker: Named operation and calculation-type fields
                for each worker, ordered by worker ID.

        Raises:
            ValueError: If there is not one operation per requested worker.
        """

        if len(parameters_per_worker) != number_of_workers_required:
            raise ValueError("one worker configuration is required per worker")

        self.workers = {
            worker_id: WorkerNode(worker_id, parameters)
            for worker_id, parameters in enumerate(parameters_per_worker)
        }

    def dispatch(
        self,
        encoded_information: Sequence[EncodedInformation],
    ) -> list[WorkerResult]:
        """Run workers and omit any worker that returns a failed answer.

        Args:
            encoded_information: One encoded-information record per worker.

        Returns:
            Results from workers whose operations returned answers.
        """
        worker_information = self._consolidate_information(encoded_information)
        return self._run_workers(worker_information)

    @staticmethod
    def _consolidate_information(
        encoded_information: Sequence[EncodedInformation],
    ) -> dict[int, list[list[Information]]]:
        """Apply coefficients and group information by worker and code row.

        Args:
            encoded_information: One encoded-information record per worker.

        Returns:
            Coefficient-weighted information grouped by worker and code row.
        """
        return {
            encoded.information_id: [
                [item.information * item.coefficient for item in row]
                for row in encoded.information
            ]
            for encoded in encoded_information
        }

    def _run_workers(
        self,
        worker_information: Mapping[int, Sequence[Sequence[Information]]],
    ) -> list[WorkerResult]:
        """Execute all workers from their consolidated information inputs.

        Args:
            worker_information: Inputs grouped by worker ID and code row.

        Returns:
            Results from workers that returned answers for all assigned rows.
        """
        worker_results = []
        for worker_id, information_rows in worker_information.items():
            worker = self.workers[worker_id]
            values = []
            for information in information_rows:
                value = worker.execute(information=information)
                if value is None:
                    break
                values.append(value)
            else:
                worker_results.append(WorkerResult(worker_id, values))
        return worker_results
