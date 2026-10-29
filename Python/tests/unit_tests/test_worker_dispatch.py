"""Tests for worker dispatch and worker-level straggler handling."""

import unittest
from collections.abc import Sequence
from functools import partial
from typing import Any, cast

import numpy as np
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.factory import create_code_family
from frame_codes.coding_scheme.encoding.encoder import Encoder
from frame_codes.coding_scheme.state_utils.information import Information
from frame_codes.coding_scheme.worker_nodes.operation_wrappers import (
    with_deterministic_straggler,
)
from frame_codes.coding_scheme.worker_nodes.worker_dispatcher import WorkerDispatcher
from frame_codes.coding_scheme.worker_nodes.worker_node import (
    WorkerNodeConfig,
    WorkerOperation,
)


class TestWorkerDispatch(unittest.TestCase):
    """Verify worker results retain Circulant's paired-row ownership."""

    @staticmethod
    def _sum_information(
        information: Sequence[Information],
        calculation_type: np.dtype,
        **kwargs: Any,
    ) -> np.ndarray:
        """Sum a worker's coefficient-weighted information values.

        Args:
            information: Coefficient-weighted values for one code row.
            calculation_type: Dtype configured for this worker.

        Returns:
            The worker's result for the code row.
        """
        return np.asarray(sum(information))

    def test_straggling_worker_drops_all_of_its_row_pair(self) -> None:
        """A straggling worker contributes none of its paired row results."""
        parameters = CodeParameters(m=4, n=8)
        code_family = create_code_family("Circulant Permutation")
        code_family.generate_code(parameters)
        matrix = np.arange(80, dtype=np.float64).reshape(16, 5)
        information_set = list(
            np.split(
                matrix,
                code_family.number_of_encoded_sets(),
                axis=0,
            )
        )
        encoded_workers = Encoder.encode(
            information_set,
            code_family,
        )
        worker_count = code_family.number_of_workers_required()
        operation = with_deterministic_straggler(self._sum_information)
        worker_parameters = [
            WorkerNodeConfig(
                operation=cast(
                    WorkerOperation,
                    partial(operation, straggle=(worker_id == 0)),
                ),
                calculation_type=code_family.get_code().dtype,
            )
            for worker_id in range(worker_count)
        ]
        dispatcher = WorkerDispatcher()
        dispatcher.setup(worker_count, worker_parameters)

        results = dispatcher.dispatch(encoded_workers)

        self.assertEqual(len(results), worker_count - 1)
        self.assertEqual(
            [result.worker_id for result in results],
            list(range(1, worker_count)),
        )

    def test_dispatch_returns_expected_rows_for_each_configured_worker(self) -> None:
        """Dispatch executes each worker and preserves its calculated rows."""
        parameters = CodeParameters(m=4, n=8)
        code_family = create_code_family("Circulant Permutation")
        code_family.generate_code(parameters)
        code = code_family.get_code()
        matrix = np.arange(80, dtype=np.float64).reshape(16, 5)
        information_set = list(
            np.split(
                matrix,
                code_family.number_of_encoded_sets(),
                axis=0,
            )
        )
        encoded_workers = Encoder.encode(information_set, code_family)
        worker_count = code_family.number_of_workers_required()
        dispatcher = WorkerDispatcher()
        dispatcher.setup(
            worker_count,
            [
                WorkerNodeConfig(
                    operation=self._sum_information,
                    calculation_type=code.dtype,
                )
                for _ in range(worker_count)
            ],
        )

        results = dispatcher.dispatch(encoded_workers)

        self.assertEqual(
            [result.worker_id for result in results], list(range(worker_count))
        )
        for result in results:
            expected_rows = [
                np.asarray(
                    sum(
                        information * code[code_row, set_index]
                        for set_index, information in enumerate(information_set)
                    )
                )
                for code_row in range(2 * result.worker_id, 2 * result.worker_id + 2)
            ]
            self.assertEqual(len(result.values), 2)
            for actual, expected in zip(result.values, expected_rows):
                np.testing.assert_allclose(actual, expected, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
