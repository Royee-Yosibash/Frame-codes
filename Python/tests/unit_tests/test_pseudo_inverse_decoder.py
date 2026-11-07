"""Tests for pseudoinverse decoding across code-family row layouts."""

import unittest

import numpy as np

from frame_codes.code_families.factory import create_code_family
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.decoding.decoders.pseudo_inverse import (
    PseudoInverseDecoder,
)
from frame_codes.coding_scheme.worker_nodes.worker_result import WorkerResult


class TestPseudoInverseDecoder(unittest.TestCase):
    """Verify recovery from surviving workers in regular and paired codes."""

    def test_decodes_arbitrary_worker_output_shapes(self):
        """The pseudoinverse recovers one function output per encoded set."""
        family = create_code_family("LPF")
        family.generate_code(CodeParameters(m=2, n=4))
        code = family.get_code()
        expected_outputs = [
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            np.array([[-2.0, 4.0], [3.0, -1.0]]),
        ]
        worker_ids = [3, 0, 2]
        worker_results = [
            WorkerResult(
                worker_id=worker_id,
                values=[
                    sum(
                        code[worker_id, set_index] * output
                        for set_index, output in enumerate(expected_outputs)
                    )
                ],
            )
            for worker_id in worker_ids
        ]

        decoded = PseudoInverseDecoder().decode(worker_results, family)

        self.assertEqual(len(decoded), len(expected_outputs))
        for actual, expected in zip(decoded, expected_outputs):
            self.assertEqual(actual.shape, expected.shape)
            np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_rejects_survivors_that_cannot_recover_all_sets(self):
        """A rank-deficient survivor code cannot decode all encoded sets."""
        family = create_code_family("LPF")
        family.generate_code(CodeParameters(m=2, n=4))
        worker_results = [WorkerResult(worker_id=0, values=[np.array([1.0])])]

        with self.assertRaises(ValueError):
            PseudoInverseDecoder().decode(worker_results, family)

        with self.assertRaises(ValueError):
            PseudoInverseDecoder().decode(worker_results, family)

    def test_decodes_circulant_workers_with_adjacent_row_pairs(self):
        """The decoder uses both code rows assigned to each Circulant worker."""
        family = create_code_family("Circulant Permutation")
        family.generate_code(CodeParameters(m=2, n=3))
        expected_outputs = [
            np.array([[1.0, 2.0], [3.0, 4.0]]),
            np.array([[-2.0, 4.0], [3.0, -1.0]]),
            np.array([[0.5, -1.5], [2.5, 3.5]]),
            np.array([[5.0, 2.0], [-1.0, 0.0]]),
        ]
        worker_results = []
        for worker_id in (1, 0):
            code_rows = family.get_code_rows_for_worker(worker_id)
            worker_results.append(
                WorkerResult(
                    worker_id=worker_id,
                    values=[
                        sum(
                            coefficient * output
                            for coefficient, output in zip(code_row, expected_outputs)
                        )
                        for code_row in code_rows
                    ],
                )
            )

        decoded = PseudoInverseDecoder().decode(worker_results, family)

        self.assertEqual(len(decoded), len(expected_outputs))
        for actual, expected in zip(decoded, expected_outputs):
            np.testing.assert_allclose(actual, expected, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
