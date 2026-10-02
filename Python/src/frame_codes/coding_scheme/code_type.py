"""Recognized frame-code family names."""

from enum import StrEnum


class CodeType(StrEnum):
    """Recognized code families."""

    LPF = "LPF"
    BPF = "BPF"
    WISHART = "Wishart"
    OMITTED_VANDERMONDE = "Omitted Vandermonde"
    NON_CONSECUTIVE_POWERS = "Non consecutive powers"
    ORTHO_MAT_DOT = "OrthoMatDot"
    CIRCULANT_PERMUTATION = "Circulant Permutation"
    BCH = "BCH"
    REED_SOLOMON = "Reed Solomon"
