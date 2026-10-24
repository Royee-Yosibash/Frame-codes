"""Pseudoinverse decoding for code-family worker outputs."""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily
from frame_codes.coding_scheme.decoding.abstract_decoder import AbstractDecoder
from frame_codes.coding_scheme.worker_nodes.worker_result import WorkerResult
from frame_codes.utils.numpy_types import as_inexact_array


class PseudoInverseDecoder(AbstractDecoder):
    """Decode code outputs with the pseudoinverse of surviving worker rows."""

    def decode(
        self,
        worker_results: Sequence[WorkerResult],
        code_family: CodeFamily,
    ) -> list[NDArray]:
        """Decode one output per encoded set from surviving workers.

        Args:
            worker_results: Answers from workers that returned results.
            code_family: Generated family.

        Returns:
            Decoded worker-function outputs, one for each encoded set.

        Raises:
            ValueError: If worker IDs or outputs do not match a regular code,
                or the surviving code rows cannot recover every encoded set.
        """
        worker_code_rows = []
        worker_outputs = []
        for result in worker_results:
            worker_code_rows.extend(
                code_family.get_code_rows_for_worker(result.worker_id)
            )
            worker_outputs.extend(
                as_inexact_array(value) for value in result.values
            )
        surviving_code = as_inexact_array(np.asarray(worker_code_rows))
        if np.linalg.matrix_rank(surviving_code) < code_family.number_of_encoded_sets():
            raise ValueError("Surviving workers cannot recover every encoded set")

        decoder_matrix = np.linalg.pinv(surviving_code)
        output_stack = np.stack(worker_outputs, axis=0)
        decoded_outputs = np.tensordot(decoder_matrix, output_stack, axes=(1, 0))
        return [np.asarray(output) for output in decoded_outputs]
