"""Wrappers for modeling worker straggle outcomes."""

from collections.abc import Sequence
from functools import wraps
from typing import Any

import numpy as np

from frame_codes.coding_scheme.state_utils.information import Information
from frame_codes.coding_scheme.worker_nodes.worker_node import WorkerOperation


def with_probabilistic_straggler(
    operation: WorkerOperation,
    straggler_probability: float,
    rng: np.random.Generator | None = None,
) -> WorkerOperation:
    """Wrap an operation with a probabilistic straggler outcome.

    Args:
        operation: Worker operation to invoke when the worker responds.
        straggler_probability: Probability that the worker returns no answer,
            in the inclusive range from zero to one.
        rng: Optional generator controlling the straggler draws.

    Returns:
        A callable that returns `None` when it straggles, otherwise the
        operation's result.

    Raises:
        ValueError: If `straggler_probability` is outside `[0, 1]`.
        TypeError: If `rng` is not a NumPy random generator.
    """
    if not 0 <= straggler_probability <= 1:
        raise ValueError("straggler_probability must be between 0 and 1")
    if rng is not None and not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")
    random_generator = np.random.default_rng() if rng is None else rng

    @wraps(operation)
    def wrapped_operation(
        *,
        information: Sequence[Information],
        calculation_type: np.dtype,
        **kwargs: Any,
    ) -> Any:
        """Return no answer on a straggler draw or run the operation.

        Args:
            information: Information supplied to the worker operation.
            calculation_type: Dtype configured for the worker.
            **kwargs: Keyword arguments forwarded to the operation.

        Returns:
            `None` for a straggler, otherwise the operation's result.
        """
        if random_generator.random() < straggler_probability:
            return None
        return operation(
            information=information,
            calculation_type=calculation_type,
            **kwargs,
        )

    return wrapped_operation


def with_deterministic_straggler(
    operation: WorkerOperation,
) -> WorkerOperation:
    """Wrap an operation with an explicit per-call straggler flag.

    Args:
        operation: Worker operation to invoke when `straggle` is false.

    Returns:
        A callable that accepts a keyword-only `straggle` flag and returns
        `None` when it is true.
    """

    @wraps(operation)
    def wrapped_operation(
        *,
        information: Sequence[Information],
        calculation_type: np.dtype,
        straggle: bool = False,
        **kwargs: Any,
    ) -> Any:
        """Return no answer when explicitly marked as a straggler.

        Args:
            information: Information supplied to the worker operation.
            calculation_type: Dtype configured for the worker.
            straggle: Whether this worker should return no answer.
            **kwargs: Keyword arguments forwarded to the operation.

        Returns:
            `None` when `straggle` is true, otherwise the operation's result.

        Raises:
            TypeError: If `straggle` is not a boolean.
        """
        if not isinstance(straggle, bool):
            raise TypeError("straggle must be a bool")
        if straggle:
            return None
        return operation(
            information=information,
            calculation_type=calculation_type,
            **kwargs,
        )

    return wrapped_operation
