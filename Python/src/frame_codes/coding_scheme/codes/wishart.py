"""Wishart-style Gaussian code family."""

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.code_family import RandomCodeFamily


class WishartCodeFamily(RandomCodeFamily):
    """Real Gaussian code family with node-count normalization."""

    def _create_code(
        self,
        n_nodes: int,
        message_dimension: int,
    ) -> NDArray[np.float64]:
        """Draw and normalize a real Gaussian code matrix.

        Args:
            n_nodes: Number of encoded outputs.
            message_dimension: Number of message elements.

        Returns:
            A real Gaussian generator matrix.
        """
        return self._rng.standard_normal((n_nodes, message_dimension)) / np.sqrt(n_nodes)
