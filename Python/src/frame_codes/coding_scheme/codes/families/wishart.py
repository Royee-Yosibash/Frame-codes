"""Wishart-style Gaussian code families."""

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import RandomCodeFamily
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class WishartCodeFamily(RandomCodeFamily):
    """Real Gaussian code families with node-count normalization."""

    def _create_code(
        self,
        parameters: CodeParameters,
    ) -> NDArray[np.float64]:
        """Draw and normalize a real Gaussian code matrix.

        Args:
            parameters: Shared code dimensions and normalization setting.

        Returns:
            A real Gaussian generator matrix.
        """
        return self._rng.standard_normal((parameters.n, parameters.m)) / np.sqrt(
            parameters.n
        )
