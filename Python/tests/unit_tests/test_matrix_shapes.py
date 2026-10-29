"""Unit tests for matrix padding and partitioning utilities."""

import unittest

import numpy as np
from frame_codes.numerics.matrix_shapes import (
    fix_matrix_dimensions,
    pad_to_partition,
    partition_matrix,
)


class TestPadToPartition(unittest.TestCase):
    """Test dtype-preserving matrix padding."""

    def test_pads_rows_with_integer_zeros(self) -> None:
        """Row padding produces a divisible shape and preserves integer dtype."""
        matrix = np.array([[1, 2], [3, 4], [5, 6]], dtype=np.int32)

        padded = pad_to_partition(matrix, partition_count=2, axis="row")

        self.assertEqual(padded.shape, (4, 2))
        self.assertEqual(padded.dtype, np.int32)
        np.testing.assert_array_equal(padded[-1], [0, 0])

    def test_pads_columns_with_complex_zeros(self) -> None:
        """Column padding preserves complex dtype and existing values."""
        matrix = np.array([[1 + 2j], [3 + 4j]], dtype=np.complex64)

        padded = pad_to_partition(matrix, partition_count=3, axis="column")

        self.assertEqual(padded.shape, (2, 3))
        self.assertEqual(padded.dtype, np.complex64)
        np.testing.assert_array_equal(padded[:, 0], matrix[:, 0])
        np.testing.assert_array_equal(
            padded[:, 1:], np.zeros((2, 2), dtype=np.complex64)
        )


class TestPartitionMatrix(unittest.TestCase):
    """Test equal row and column partitioning."""

    def test_partitions_rows_equally(self) -> None:
        """Row partitions preserve their order and dtype."""
        matrix = np.arange(12, dtype=np.int16).reshape(4, 3)

        partitions = partition_matrix(matrix, partition_count=2, axis="row")

        self.assertEqual(len(partitions), 2)
        self.assertTrue(all(part.dtype == np.int16 for part in partitions))
        np.testing.assert_array_equal(np.vstack(partitions), matrix)

    def test_partitions_columns_equally(self) -> None:
        """Column partitions preserve their order and dtype."""
        matrix = np.arange(12, dtype=np.float32).reshape(3, 4)

        partitions = partition_matrix(matrix, partition_count=2, axis="column")

        self.assertEqual(len(partitions), 2)
        self.assertTrue(all(part.dtype == np.float32 for part in partitions))
        np.testing.assert_array_equal(np.hstack(partitions), matrix)

    def test_rejects_uneven_partitioning(self) -> None:
        """The selected dimension must divide evenly."""
        with self.assertRaisesRegex(ValueError, "not divisible"):
            partition_matrix(np.ones((3, 2)), partition_count=2, axis="row")


class TestFixMatrixDimensions(unittest.TestCase):
    """Test code-specific padding conventions."""

    def test_pads_rows_of_a_and_columns_of_b_for_regular_codes(self) -> None:
        """Standard code types pad outer matrix dimensions."""
        matrix_a = np.ones((3, 2), dtype=np.int32)
        matrix_b = np.ones((2, 3), dtype=np.float32)

        padded_a, padded_b = fix_matrix_dimensions(
            matrix_a, matrix_b, "Non consecutive powers", 2, 2
        )

        self.assertEqual(padded_a.shape, (4, 2))
        self.assertEqual(padded_b.shape, (2, 4))
        self.assertEqual(padded_a.dtype, np.int32)
        self.assertEqual(padded_b.dtype, np.float32)

    def test_pads_orthomatdot_shared_dimension_to_common_multiple(self) -> None:
        """OrthoMatDot pads both sides of the product's shared dimension."""
        matrix_a = np.ones((2, 5), dtype=np.int16)
        matrix_b = np.ones((5, 1), dtype=np.complex64)

        padded_a, padded_b = fix_matrix_dimensions(
            matrix_a, matrix_b, "OrthoMatDot", 2, 3
        )

        self.assertEqual(padded_a.shape, (2, 6))
        self.assertEqual(padded_b.shape, (6, 1))
        self.assertEqual(padded_a.dtype, np.int16)
        self.assertEqual(padded_b.dtype, np.complex64)
        np.testing.assert_array_equal(padded_a[:, -1], [0, 0])
        np.testing.assert_array_equal(padded_b[-1, :], [0])
        np.testing.assert_array_equal(padded_a @ padded_b, matrix_a @ matrix_b)

    def test_rejects_incompatible_matrix_product_dimensions(self) -> None:
        """Input matrices must be conformable for multiplication."""
        with self.assertRaisesRegex(ValueError, "columns must match"):
            fix_matrix_dimensions(np.ones((2, 3)), np.ones((4, 2)), "Other", 2, 2)


if __name__ == "__main__":
    unittest.main()
