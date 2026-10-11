"""Pseudoinverse decoding for regular one-row-per-worker codes."""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.decoding.abstract_decoder import AbstractDecoder
from frame_codes.coding_scheme.worker_runtime.runtime import WorkerResult
from frame_codes.utils.numpy_types import as_inexact_array


class PseudoInverseDecoder(AbstractDecoder):
    """Decode regular code outputs with the pseudoinverse of surviving rows."""

    def decode(
        self,
        worker_results: Sequence[WorkerResult],
        code_family: CodeFamily,
    ) -> list[NDArray]:
        """Decode one output per encoded set from surviving regular workers.

        Args:
            worker_results: One answer from each worker that returned a result.
            code_family: Generated regular code with one row per worker.

        Returns:
            Decoded worker-function outputs, one for each encoded set.

        Raises:
            ValueError: If worker IDs or outputs do not match a regular code,
                or the surviving code rows cannot recover every encoded set.
        """
        code = as_inexact_array(code_family.get_code())
        worker_ids = [result.worker_id for result in worker_results]
        surviving_code = code[worker_ids]
        if np.linalg.matrix_rank(surviving_code) < code_family.number_of_encoded_sets():
            raise ValueError("Surviving workers cannot recover every encoded set")

        decoder_matrix = np.linalg.pinv(surviving_code)

        worker_outputs = [
            as_inexact_array(result.values[0]) for result in worker_results
        ]

        output_stack = np.stack(worker_outputs, axis=0)
        decoded_outputs = np.tensordot(decoder_matrix, output_stack, axes=(1, 0))
        return [np.asarray(output) for output in decoded_outputs]
