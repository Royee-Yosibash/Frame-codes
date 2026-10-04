"""Wishart-style Gaussian code families."""

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import RandomCodeFamily


class WishartCodeFamily(RandomCodeFamily):
    """Real Gaussian code families with node-count normalization."""

    def _create_code(
        self,
        n: int,
        m: int,
    ) -> NDArray[np.float64]:
        """Draw and normalize a real Gaussian code matrix.

        Args:
            n: Number of encoded outputs.
            m: Number of message elements.

        Returns:
            A real Gaussian generator matrix.
        """
        return self._rng.standard_normal((n, m)) / np.sqrt(n)
