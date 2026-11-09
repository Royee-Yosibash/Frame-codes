"""End-to-end pipeline test for encoding, worker execution, and decoding."""

import unittest
from collections.abc import Sequence
from functools import partial
from typing import Any, cast

import numpy as np

from frame_codes.code_families.factory import create_code_family
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.decoding.decoders.pseudo_inverse import (
    PseudoInverseDecoder,
)
from frame_codes.coding_scheme.encoding.encoder import Encoder
from frame_codes.coding_scheme.state_utils.information import Information
from frame_codes.coding_scheme.worker_nodes.operation_wrappers import (
    with_deterministic_straggler,
)
from frame_codes.coding_scheme.worker_nodes.worker_dispatcher import (
    WorkerDispatcher,
)
from frame_codes.coding_scheme.worker_nodes.worker_node import (
    WorkerNodeConfig,
    WorkerOperation,
)
from frame_codes.comparison.metrics import relative_frobenius_error
from frame_codes.distributed_matrix_multiplication.operations import matrix_product


def apply_worker_function(
    information: Sequence[Information],
    calculation_type: np.dtype,
    **kwargs: Any,
) -> np.ndarray:
    """Sum coefficient-weighted information for one worker row.

    Args:
        information: Coefficient-weighted information values.
        calculation_type: Dtype configured for the worker.

    Returns:
        The worker's encoded answer.
    """
    return np.asarray(sum(information))


class TestWorkerCodingPipeline(unittest.TestCase):
    """Check recovery across the actual encode-to-worker-to-decode path."""

    def setUp(self):
        """
        Set up the pipeline for testing.

        It is important that the test sets this up once because it is important to check sanitation is handled correctly
        between usages since it is crucial for allowing a clean and robust comparison suite
        """

        self.encoder = Encoder()
        self.dispatcher = WorkerDispatcher()
        self.decoder = PseudoInverseDecoder()

    def test_regular_code_recovers_information_after_worker_erasure(self) -> None:
        """The full pipeline recovers original items after one worker is lost."""
        code_family = create_code_family("LPF")
        code_family.generate_code(CodeParameters(m=2, n=4))
        information_set = [
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            np.array([[-2.0, 4.0], [3.0, -1.0]]),
        ]


        encoded_workers = self.encoder.encode(information_set, code_family)
        
        worker_count = code_family.number_of_workers_required()
        self.dispatcher.setup(
            number_of_workers_required=worker_count,
            parameters_per_worker=
                [
                WorkerNodeConfig(
                    operation=cast(
                        WorkerOperation,
                        partial(
                            with_deterministic_straggler(apply_worker_function),
                            straggle=(worker_id == 0),
                        ),
                    ),
                    calculation_type=code_family.get_code().dtype,
                )
                    for worker_id in range(worker_count)
                ],
        )
        worker_results = self.dispatcher.dispatch(encoded_workers)
        decoded_information = self.decoder.decode(worker_results, code_family)

        self.assertEqual(len(worker_results), worker_count - 1)
        self.assertNotIn(0, [result.worker_id for result in worker_results])
        self.assertEqual(len(decoded_information), len(information_set))
        for decoded, original in zip(decoded_information, information_set):
            np.testing.assert_allclose(decoded, original, atol=1e-12)

    def test_complex_code_recovers_complex_information(self) -> None:
        """The full pipeline preserves complex values through worker execution."""
        code_family = create_code_family("LPF")
        code_family.generate_code(CodeParameters(m=2, n=4))
        information_set = [
            np.array([[1.0 + 2.0j, 3.0 - 1.0j]]),
            np.array([[-2.0 + 0.5j, 4.0 + 3.0j]]),
        ]
        encoded_workers = self.encoder.encode(information_set, code_family)
        worker_count = code_family.number_of_workers_required()
        calculation_type = code_family.get_code().dtype
        self.dispatcher.setup(
            worker_count,
            [
                WorkerNodeConfig(
                    operation=apply_worker_function,
                    calculation_type=calculation_type,
                )
                for _ in range(worker_count)
            ],
        )

        worker_results = self.dispatcher.dispatch(encoded_workers)
        decoded_information = self.decoder.decode(worker_results, code_family)

        self.assertEqual(len(decoded_information), len(information_set))
        for decoded, original in zip(decoded_information, information_set):
            np.testing.assert_allclose(decoded, original, atol=1e-12)

    def test_mat_vect_multiplication_in_distributed_framework(self) -> None:
        """The pipeline computes `A @ v` and recovers the expected vector b."""
        code_family = create_code_family("LPF")
        code_family.generate_code(CodeParameters(m=3, n=4))
        matrix_a = np.array(
            [
                [1.0, 2.0, 0.0],
                [0.0, -1.0, 3.0],
                [4.0, 1.0, 2.0],
            ]
        )
        vector_v = np.array([2.0, -1.0, 3.0])
        expected_b = np.array([0.0, 10.0, 13.0])
        np.testing.assert_allclose(matrix_a @ vector_v, expected_b)
        information_set = list(
            np.split(matrix_a, code_family.number_of_encoded_sets(), axis=0)
        )

        encoded_workers = self.encoder.encode(information_set, code_family)
        worker_count = code_family.number_of_workers_required()
        self.dispatcher.setup(
            worker_count,
            [
                WorkerNodeConfig(
                    operation=partial(matrix_product, rhs=vector_v),
                    calculation_type=code_family.get_code().dtype,
                )
                for _ in range(worker_count)
            ],
        )
        worker_results = self.dispatcher.dispatch(encoded_workers)
        decoded_rows = self.decoder.decode(worker_results, code_family)
        result = np.concatenate(decoded_rows, axis=0)

        np.testing.assert_allclose(result, expected_b, atol=1e-12)

    def test_happy_path_matrix_vector_pipeline_has_low_relative_error(self) -> None:
        """With no stragglers, the pipeline computes `A @ v` accurately."""
        code_family = create_code_family("LPF")
        code_family.generate_code(CodeParameters(m=3, n=4))
        matrix_a = np.eye(3)
        vector_v = np.ones((3, 1))
        expected_b = matrix_a @ vector_v
        information_set = list(
            np.split(matrix_a, code_family.number_of_encoded_sets(), axis=0)
        )

        encoded_workers = self.encoder.encode(information_set, code_family)
        worker_count = code_family.number_of_workers_required()
        self.dispatcher.setup(
            worker_count,
            [
                WorkerNodeConfig(
                    operation=partial(matrix_product, rhs=vector_v),
                    calculation_type=code_family.get_code().dtype,
                )
                for _ in range(worker_count)
            ],
        )
        worker_results = self.dispatcher.dispatch(encoded_workers)
        decoded_rows = self.decoder.decode(worker_results, code_family)
        actual_b = np.concatenate(decoded_rows, axis=0)
        self.assertLessEqual(
            relative_frobenius_error(actual_b, expected_b),
            1e-12,
        )


if __name__ == "__main__":
    unittest.main()
