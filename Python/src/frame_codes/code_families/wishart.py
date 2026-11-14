"""Wishart-style Gaussian code family implementation."""

import numpy as np
from codedcompute.codes.code_family import RandomCodeFamily
from codedcompute.codes.code_parameters import CodeParameters
from numpy.typing import NDArray


class WishartCodeFamily(RandomCodeFamily):
    """Real Gaussian code families with node-count normalization."""

    def _generate_new_code(
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
