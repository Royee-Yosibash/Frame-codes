"""Factory for selecting concrete code-families implementations."""

import numpy as np

from frame_codes.code_families.circulant_permutation import (
    CirculantPermutationCodeFamily,
)
from frame_codes.code_families.code_type import CodeType
from frame_codes.code_families.fourier import (
    BandPassFourierCodeFamily,
    LowPassFourierCodeFamily,
)
from frame_codes.code_families.orthomatdot import (
    OrthoMatDotCodeFamily,
)
from frame_codes.code_families.simple_band_circulant import (
    SimpleBandCirculantCode,
)
from frame_codes.code_families.vandermonde import (
    VandermondeCodeFamily,
)
from frame_codes.code_families.wishart import WishartCodeFamily
from frame_codes.coding_scheme.codes.code_family import CodeFamily, RandomCodeFamily

CODE_FAMILY_TYPES: dict[CodeType, type[CodeFamily]] = {
    CodeType.LPF: LowPassFourierCodeFamily,
    CodeType.BPF: BandPassFourierCodeFamily,
    CodeType.WISHART: WishartCodeFamily,
    CodeType.OMITTED_VANDERMONDE: VandermondeCodeFamily,
    CodeType.NON_CONSECUTIVE_POWERS: VandermondeCodeFamily,
    CodeType.ORTHO_MAT_DOT: OrthoMatDotCodeFamily,
    CodeType.CIRCULANT_PERMUTATION: CirculantPermutationCodeFamily,
    CodeType.SIMPLE_BAND_CIRCULANT: SimpleBandCirculantCode,
}


def create_code_family(
    code_type: CodeType | str,
    rng: np.random.Generator | None = None,
) -> CodeFamily:
    """Create the concrete families registered for a code type.

    Args:
        code_type: Code families enum member or its string value.
        rng: Optional random generator for stochastic code families.

    Returns:
        An instance of the registered code-families class.

    Raises:
        NotImplementedError: If a recognized finite-field families has no
            implementation.
        TypeError: If `rng` is not a NumPy random generator.
        ValueError: If the code type is not stochastic but `rng` is supplied.
    """
    code_type = CodeType(code_type)
    family_type = CODE_FAMILY_TYPES.get(code_type)
    if family_type is None:
        raise NotImplementedError(
            f"{code_type} requires finite-field code support that is not implemented"
        )
    if issubclass(family_type, RandomCodeFamily):
        return family_type(rng)
    if rng is not None:
        raise ValueError(f"rng is not used by code type {code_type}")
    return family_type()
