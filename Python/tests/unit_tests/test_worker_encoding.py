"""Tests for family-owned worker-specific encoding."""

import unittest

import numpy as np
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.factory import create_code_family
from frame_codes.coding_scheme.encoding.encoder import Encoder
from frame_codes.coding_scheme.state_utils.information import (
    WeightedInformation,
)


class TestWorkerEncoding(unittest.TestCase):
    """Verify code-family coefficient selection and encoder delegation."""

    def test_standard_family_encodes_one_row_for_requested_worker(self) -> None:
        """A standard worker receives its corresponding code row."""
        family = create_code_family("LPF")
        family.generate_code(CodeParameters(m=2, n=3))
        code = family.get_code()
        information_set = ["first", "second"]

        encoded = family.encode_for_worker(1, information_set)

        self.assertEqual(len(encoded), 1)
        self.assertEqual([item.information for item in encoded[0]], information_set)
        np.testing.assert_allclose(
            [item.coefficient for item in encoded[0]], code[1, :]
        )

    def test_circulant_worker_encodes_its_adjacent_code_rows(self) -> None:
        """A Circulant worker receives coefficients from its row pair."""
        family = create_code_family("Circulant Permutation")
        family.generate_code(CodeParameters(m=2, n=3))
        code = family.get_code()
        information_set = [f"set-{index}" for index in range(code.shape[1])]

        encoded = family.encode_for_worker(1, information_set)

        self.assertEqual(len(encoded), 2)
        for row_index, encoded_row in enumerate(encoded, start=2):
            self.assertEqual(
                [item.information for item in encoded_row], information_set
            )
            np.testing.assert_allclose(
                [item.coefficient for item in encoded_row], code[row_index, :]
            )

    def test_worker_encoded_rows_are_independent_lists(self) -> None:
        """Appending to one encoded row leaves another row unchanged."""
        family = create_code_family("Circulant Permutation")
        family.generate_code(CodeParameters(m=2, n=3))
        information_set = [
            f"set-{index}" for index in range(family.number_of_encoded_sets())
        ]

        encoded = family.encode_for_worker(0, information_set)
        second_row_before = encoded[1].copy()
        encoded[0].append(WeightedInformation("extra", 1.0))

        self.assertEqual(encoded[1], second_row_before)

    def test_encoder_creates_one_record_per_worker(self) -> None:
        """The encoder delegates each worker's payload to the family."""
        family = create_code_family("Circulant Permutation")
        family.generate_code(CodeParameters(m=2, n=3))
        information_set = [
            f"set-{index}" for index in range(family.number_of_encoded_sets())
        ]

        encoded_workers = Encoder.encode(information_set, family)

        self.assertEqual(len(encoded_workers), family.number_of_workers_required())
        self.assertEqual(
            [worker.information_id for worker in encoded_workers], [0, 1, 2]
        )
        self.assertTrue(all(len(worker.information) == 2 for worker in encoded_workers))

    def test_rejects_invalid_worker_and_information_count(self) -> None:
        """The family rejects worker IDs and item counts outside its code."""
        family = create_code_family("LPF")
        family.generate_code(CodeParameters(m=2, n=3))

        with self.assertRaisesRegex(ValueError, "worker range"):
            family.encode_for_worker(3, ["first", "second"])
        with self.assertRaisesRegex(ValueError, "information item count"):
            family.encode_for_worker(0, ["only one"])


if __name__ == "__main__":
    unittest.main()
