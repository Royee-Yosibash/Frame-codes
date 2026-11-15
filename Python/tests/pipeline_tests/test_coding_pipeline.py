"""End-to-end pipeline test for encoding, worker execution, and decoding."""

import unittest
from functools import partial

import numpy as np
from codedcompute.codes.code_parameters import CodeParameters
from codedcompute.decoding.decoders.pseudo_inverse import (
    PseudoInverseDecoder,
)
from codedcompute.encoding.encoder import Encoder
from codedcompute.worker_nodes.worker_dispatcher import (
    WorkerDispatcher,
)
from codedcompute.worker_nodes.worker_node import (
    WorkerNodeConfig,
)

from frame_codes.code_families.factory import create_code_family
from frame_codes.comparison.metrics import relative_frobenius_error
from frame_codes.distributed_matrix_multiplication.operations import matrix_product


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
