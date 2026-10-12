"""Abstract interface for decoding worker results."""

from abc import ABC, abstractmethod
from collections.abc import Sequence
from typing import Any

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.worker_runtime.runtime import WorkerResult


class AbstractDecoder(ABC):
    """Base interface for decoding results returned by workers."""

    @abstractmethod
    def decode(
        self,
        worker_results: Sequence[WorkerResult],
        code_family: CodeFamily,
    ) -> Any:
        """Decode the available worker results using a generated code family.

        Args:
            worker_results: Results from workers that returned an answer.
            code_family: Family containing the generated code.

        Returns:
            The output defined by the concrete decoding strategy.
        """
