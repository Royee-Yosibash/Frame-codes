"""End-to-end pipeline test for encoding, worker execution, and decoding."""

import unittest
from collections.abc import Sequence
from functools import partial

import numpy as np
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.factory import create_code_family
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


def apply_worker_function(
    information: Sequence[Information],
    calculation_type: np.dtype,
) -> np.ndarray:
    """Sum coefficient-weighted information for one worker row.

    Args:
        information: Coefficient-weighted information values.
        calculation_type: Dtype configured for the worker.

    Returns:
        The worker's encoded answer.
    """
    return sum(information)


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
            worker_count,
            [
                (
                    partial(
                        with_deterministic_straggler(apply_worker_function),
                        straggle=(worker_id == 0),
                    ),
                    code_family.get_code().dtype,
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


if __name__ == "__main__":
    unittest.main()
