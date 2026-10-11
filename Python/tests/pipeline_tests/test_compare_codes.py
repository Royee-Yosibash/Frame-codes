"""Pipeline test for information encoding, worker execution, and decoding."""

import unittest

import numpy as np
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.factory import create_code_family
from frame_codes.coding_scheme.decoding.decoders.pseudo_inverse import (
    PseudoInverseDecoder,
)
from frame_codes.coding_scheme.encoding.encoder import Encoder
from frame_codes.coding_scheme.state_utils.encoded_information import (
    WeightedInformation,
)
from frame_codes.coding_scheme.worker_runtime.runtime import (
    execute_workers,
    sample_worker_results,
)


class TestWorkerCodingPipeline(unittest.TestCase):
    """Check recovery across the actual encode-to-worker-to-decode path."""

    @staticmethod
    def _apply_worker_function(
        encoded_information: list[WeightedInformation],
    ) -> np.ndarray:
        """Apply the worker's weighted information sum.

        Args:
            encoded_information: Weighted information items for one worker.

        Returns:
            The worker's encoded answer.
        """
        return sum(
            item.coefficient * item.information for item in encoded_information
        )

    def test_regular_code_recovers_information_after_worker_erasure(self) -> None:
        """The full pipeline recovers original items after one worker is lost."""
        code_family = create_code_family("LPF")
        code_family.generate_code(CodeParameters(m=2, n=4))
        information_set = [
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            np.array([[-2.0, 4.0], [3.0, -1.0]]),
        ]

        encoded_workers = Encoder.encode(information_set, code_family)
        worker_results = execute_workers(
            [encoded_workers], self._apply_worker_function
        )
        surviving_results = sample_worker_results(
            worker_results,
            straggler_count=1,
            snr_db=300,
            rng=np.random.default_rng(53),
        )
        decoded_information = PseudoInverseDecoder().decode(
            surviving_results, code_family
        )

        self.assertEqual(len(worker_results), 4)
        self.assertEqual(len(surviving_results), 3)
        self.assertEqual(len(decoded_information), len(information_set))
        for decoded, original in zip(decoded_information, information_set):
            np.testing.assert_allclose(decoded, original, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
