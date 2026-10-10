"""NumPy dtype helpers shared by numerical modules."""

import numpy as np
from numpy.typing import ArrayLike, NDArray


def as_inexact_array(values: ArrayLike) -> NDArray[np.float64] | NDArray[np.complex128]:
    """Convert supported numeric inputs to a safe floating-point dtype.

    Args:
        values: Real or complex NumPy-compatible numeric values.

    Returns:
        A real array converted to `float64`, or a complex array converted to
        `complex128`.

    Raises:
        TypeError: If the input dtype is not an integer, unsigned integer,
            floating-point, or complex numeric type.
    """
    array = np.asarray(values)
    if array.dtype.kind not in "iufc":
        raise TypeError("Values must have an integer, floating, or complex dtype")
    dtype = np.complex128 if array.dtype.kind == "c" else np.float64
    return np.asarray(array, dtype=dtype)
