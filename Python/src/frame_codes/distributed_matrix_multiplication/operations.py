"""Worker operations for distributed matrix multiplication."""

from collections.abc import Sequence
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from frame_codes.coding_scheme.state_utils.information import Information


def matrix_product(
    *,
    information: Sequence[Information],
    calculation_type: np.dtype,
    rhs: ArrayLike | None = None,
    **kwargs: Any,
) -> NDArray:
    """Multiply an encoded left-hand matrix block by the supplied RHS.

    Args:
        information: Coefficient-weighted matrix blocks for one code row.
        calculation_type: Dtype configured for the worker.
    rhs: Uncoded vector or matrix shared by workers for this dispatch.
        **kwargs: Additional worker options, which are ignored.

    Returns:
        The worker's encoded matrix product.

    Raises:
        TypeError: If the right operand cannot be cast to the worker dtype
            using NumPy's same-kind casting rule.
    """
    calculation_rhs = np.asarray(rhs)
    if not np.can_cast(
        calculation_rhs.dtype,
        calculation_type,
        casting="same_kind",
    ):
        raise TypeError(
            "right operand dtype cannot be cast to worker using same-kind casting"
        )
    encoded_matrix = np.asarray(sum(information))
    calculation_rhs = calculation_rhs.astype(calculation_type, copy=False)
    return encoded_matrix @ calculation_rhs
