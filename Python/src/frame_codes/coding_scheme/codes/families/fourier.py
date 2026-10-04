"""Low-pass and band-pass Fourier code families."""
from abc import ABC
from functools import cache

import numpy as np
from numpy.typing import NDArray

from frame_codes.coding_scheme.codes.code_family import CodeFamily, RandomCodeFamily


class FourierCodeFamily(CodeFamily, ABC):
    """Shared normalized inverse-DFT construction for Fourier families."""

    @staticmethod
    @cache
    def _create_base_code(n_nodes: int) -> NDArray[np.complex128]:
        """Create and cache the full inverse-DFT matrix by node count.

        Args:
            n_nodes: Fourier matrix dimension.

        Returns:
            A normalized inverse-DFT matrix.
        """
        indices = np.arange(n_nodes)
        matrix = (
            np.exp(2j * np.pi * np.outer(indices, indices) / n_nodes)
            / np.sqrt(n_nodes)
        )
        return matrix


class LowPassFourierCodeFamily(FourierCodeFamily):
    """Deterministic low-pass Fourier code families."""

    def _create_code(
        self,
        n: int,
        m: int,
    ) -> NDArray[np.complex128]:
        """Select the first inverse-DFT columns.

        Args:
            n: Number of encoded outputs.
            m: Number of message elements.

        Returns:
            A complex generator containing the first Fourier columns.
        """
        return self._create_base_code(n)[:, :m]


class BandPassFourierCodeFamily(FourierCodeFamily, RandomCodeFamily):
    """Fourier code families using a random contiguous frequency band."""

    def _create_code(
        self,
        n: int,
        m: int,
    ) -> NDArray[np.complex128]:
        """Select a random contiguous band from the inverse DFT matrix.

        Args:
            n: Number of encoded outputs.
            m: Number of message elements.

        Returns:
            A complex generator containing the selected Fourier band.
        """
        start = int(self._rng.integers(0, n - m + 1))
        return self._create_base_code(n)[:, start : start + m]
