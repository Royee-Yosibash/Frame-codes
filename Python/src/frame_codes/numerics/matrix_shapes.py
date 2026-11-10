"""Matrix padding and partitioning utilities."""

from math import lcm
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

MatrixAxis = Literal["row", "column"]


def pad_to_partition(
    matrix: ArrayLike,
    partition_count: int,
    axis: MatrixAxis,
) -> NDArray:
    """Pad one matrix axis with zeros until it divides evenly into partitions.

    Args:
        matrix: A two-dimensional array-like matrix.
        partition_count: Positive number of equal partitions required.
        axis: Axis to pad: `row` pads rows and `column` pads columns.

    Returns:
        A matrix padded at the end of the selected axis, preserving its dtype.

    Raises:
        TypeError: If `partition_count` is not an integer.
        ValueError: If the matrix is not two-dimensional, the partition count
            is not positive, or the axis is not supported.
    """
    values = np.asarray(matrix)
    _validate_partition_count(partition_count)
    axis_index = _axis_index(axis)
    remainder = values.shape[axis_index] % partition_count
    if remainder == 0:
        return values

    padding = partition_count - remainder
    return _pad_axis_to_size(values, axis_index, values.shape[axis_index] + padding)


def partition_matrix(
    matrix: ArrayLike,
    partition_count: int,
    axis: MatrixAxis,
) -> list[NDArray]:
    """Split a matrix into equal row or column partitions.

    Args:
        matrix: A two-dimensional array-like matrix.
        partition_count: Positive number of equal partitions to produce.
        axis: Axis to partition: `row` partitions rows and `column` partitions
            columns.

    Returns:
        A list of equal-sized matrix partitions that preserve the input dtype.

    Raises:
        TypeError: If `partition_count` is not an integer.
        ValueError: If the matrix is not two-dimensional, the partition count
            is not positive, the axis is unsupported, or the selected axis
            cannot be divided evenly.
    """
    values = np.asarray(matrix)
    _validate_partition_count(partition_count)
    axis_index = _axis_index(axis)
    if values.shape[axis_index] % partition_count != 0:
        raise ValueError("matrix dimension is not divisible by partition_count")
    return list(np.split(values, partition_count, axis=axis_index))


def fix_matrix_dimensions(
    matrix_a: ArrayLike,
    matrix_b: ArrayLike,
    code_type: str,
    partitions_a: int,
    partitions_b: int,
) -> tuple[NDArray, NDArray]:
    """Pad matrix dimensions to fit the partitioning rules for a code type.

    OrthoMatDot pads the shared inner dimension of both matrices to a common
    multiple of their partition counts. Other code types pad rows of A and
    columns of B.

    Args:
        matrix_a: Left matrix in a matrix product.
        matrix_b: Right matrix in a matrix product.
        code_type: Code type selecting the padding convention.
        partitions_a: Positive partition count for matrix A.
        partitions_b: Positive partition count for matrix B.

    Returns:
        A tuple containing the padded A and B matrices. Both outputs preserve
        the dtype of their corresponding input.

    Raises:
        TypeError: If either partition count is not an integer.
        ValueError: If an input is not two-dimensional, a partition count is
            not positive, or the matrices cannot be multiplied.
    """
    values_a = np.asarray(matrix_a)
    values_b = np.asarray(matrix_b)
    _validate_partition_count(partitions_a)
    _validate_partition_count(partitions_b)
    if values_a.shape[1] != values_b.shape[0]:
        raise ValueError("matrix_a columns must match matrix_b rows")

    if code_type == "OrthoMatDot":
        common_partition = lcm(partitions_a, partitions_b)
        inner_size = _round_up(values_a.shape[1], common_partition)
        padded_a = _pad_axis_to_size(values_a, 1, inner_size)
        padded_b = _pad_axis_to_size(values_b, 0, inner_size)
        return padded_a, padded_b

    padded_a = pad_to_partition(values_a, partitions_a, "row")
    padded_b = pad_to_partition(values_b, partitions_b, "column")
    return padded_a, padded_b


def _validate_partition_count(partition_count: int) -> None:
    """Validate a partition count.

    Args:
        partition_count: Requested number of partitions.

    Raises:
        TypeError: If the partition count is not an integer.
        ValueError: If the partition count is not positive.
    """
    if isinstance(partition_count, bool) or not isinstance(
        partition_count, (int, np.integer)
    ):
        raise TypeError("partition_count must be an integer")
    if partition_count <= 0:
        raise ValueError("partition_count must be positive")


def _axis_index(axis: MatrixAxis) -> int:
    """Convert an axis name to its NumPy dimension index.

    Args:
        axis: Matrix axis name.

    Returns:
        Zero for rows or one for columns.

    Raises:
        ValueError: If the axis name is not `row` or `column`.
    """
    if axis == "row":
        return 0
    if axis == "column":
        return 1
    raise ValueError("axis must be 'row' or 'column'")


def _pad_axis_to_size(matrix: NDArray, axis: int, target_size: int) -> NDArray:
    """Pad a matrix axis to a requested size with dtype-compatible zeros.

    Args:
        matrix: Input matrix.
        axis: Axis to pad.
        target_size: Target length of the selected axis.

    Returns:
        The padded matrix, preserving the input dtype.
    """
    padding = target_size - matrix.shape[axis]
    if padding <= 0:
        return matrix
    pad_width = [(0, 0), (0, 0)]
    pad_width[axis] = (0, padding)
    return np.pad(matrix, pad_width, mode="constant", constant_values=0)


def _round_up(value: int, multiple: int) -> int:
    """Round a non-negative value up to the nearest multiple.

    Args:
        value: Value to round.
        multiple: Positive multiple to round to.

    Returns:
        The smallest multiple of `multiple` greater than or equal to `value`.
    """
    return ((value + multiple - 1) // multiple) * multiple
