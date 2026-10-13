"""Wrappers for modeling worker straggle outcomes."""

from collections.abc import Callable
from functools import wraps
from typing import Any

import numpy as np


def with_probabilistic_straggler(
    operation: Callable[..., Any],
    straggler_probability: float,
    rng: np.random.Generator | None = None,
) -> Callable[..., Any | None]:
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
    def wrapped_operation(*args: Any, **kwargs: Any) -> Any | None:
        """Return no answer on a straggler draw or run the operation.

        Args:
            *args: Positional arguments forwarded to the operation.
            **kwargs: Keyword arguments forwarded to the operation.

        Returns:
            `None` for a straggler, otherwise the operation's result.
        """
        if random_generator.random() < straggler_probability:
            return None
        return operation(*args, **kwargs)

    return wrapped_operation


def with_deterministic_straggler(
    operation: Callable[..., Any],
) -> Callable[..., Any | None]:
    """Wrap an operation with an explicit per-call straggler flag.

    Args:
        operation: Worker operation to invoke when `straggle` is false.

    Returns:
        A callable that accepts a keyword-only `straggle` flag and returns
        `None` when it is true.
    """

    @wraps(operation)
    def wrapped_operation(
        *args: Any,
        straggle: bool = False,
        **kwargs: Any,
    ) -> Any | None:
        """Return no answer when explicitly marked as a straggler.

        Args:
            *args: Positional arguments forwarded to the operation.
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
        return operation(*args, **kwargs)

    return wrapped_operation
