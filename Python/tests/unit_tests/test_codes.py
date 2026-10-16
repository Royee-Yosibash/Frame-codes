"""Unit tests for code construction and normalization."""

import unittest

import numpy as np

from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.code_type import CodeType
from frame_codes.coding_scheme.codes.factory import create_code_family
from frame_codes.coding_scheme.codes.families.fourier import (
    LowPassFourierCodeFamily,
)


class TestCodeFamilyConstruction(unittest.TestCase):
    """Test supported code-families construction."""

    def test_code_type_enum_members_are_string_values(self) -> None:
        """Code type enum members retain their external string values."""
        self.assertEqual(CodeType.NON_CONSECUTIVE_POWERS, "Non consecutive powers")

    def test_factory_selects_the_registered_family(self) -> None:
        """The factory delegates code-families selection to its registry."""
        family = create_code_family(CodeType.LPF)

        self.assertIsInstance(family, LowPassFourierCodeFamily)

    def test_cached_deterministic_matrices_can_be_regenerated_after_mutation(self) -> None:
        """Regeneration restores the cached code after current state mutation."""
        family = create_code_family("LPF")
        parameters = CodeParameters(m=3, n=8)
        family.generate_code(parameters)
        first = family.get_code()
        family.generate_code(parameters)
        expected = family.get_code().copy()
        first[0, 0] = 0

        family.generate_code(parameters)
        np.testing.assert_allclose(family.get_code(), expected)

    def test_lpf_columns_are_orthonormal(self) -> None:
        """The low-pass Fourier columns have unitary normalization."""
        family = create_code_family("LPF")
        family.generate_code(CodeParameters(m=3, n=8))
        code = family.get_code()

        np.testing.assert_allclose(code.conj().T @ code, np.eye(3), atol=1e-12)

    def test_randomized_codes_use_explicit_generators(self) -> None:
        """Equal generator states reproduce stochastic code matrices."""
        for code_type in ("BPF", "Wishart", "Non consecutive powers"):
            first_family = create_code_family(code_type, rng=np.random.default_rng(7))
            first_family.generate_code(CodeParameters(m=4, n=8))
            first = first_family.get_code()
            second_family = create_code_family(code_type, rng=np.random.default_rng(7))
            second_family.generate_code(CodeParameters(m=4, n=8))
            second = second_family.get_code()
            np.testing.assert_array_equal(first, second)

    def test_vandermonde_code_uses_powers_as_columns(self) -> None:
        """The generator has one row per node and selected powers as columns."""
        family = create_code_family("Non consecutive powers", rng=np.random.default_rng(2))
        family.generate_code(CodeParameters(m=4, n=9))
        code = family.get_code()

        self.assertEqual(code.shape, (9, 4))
        np.testing.assert_allclose(np.abs(code), 1)

    def test_orthomatdot_uses_chebyshev_recurrence(self) -> None:
        """The code columns follow first-kind Chebyshev polynomials."""
        n_nodes = 8
        family = create_code_family("OrthoMatDot")
        family.generate_code(CodeParameters(m=4, n=n_nodes))
        code = family.get_code()
        nodes = np.cos((2 * np.arange(n_nodes) + 1) * np.pi / (2 * n_nodes))

        np.testing.assert_allclose(code[:, 0], 1 / np.sqrt(2))
        np.testing.assert_allclose(code[:, 1], nodes)
        np.testing.assert_allclose(code[:, 2], 2 * nodes**2 - 1)

    def test_circulant_permutation_expands_worker_and_set_dimensions(self) -> None:
        """Circulant expands each worker and encoded set into a pair."""
        family = create_code_family("Circulant Permutation")
        family.generate_code(CodeParameters(m=4, n=8))
        code = family.get_code()

        self.assertEqual(code.shape, (16, 8))
        self.assertTrue(np.isrealobj(code))
        self.assertEqual(family.number_of_workers_required(), 8)
        self.assertEqual(family.number_of_encoded_sets(), 8)

        family.generate_code(CodeParameters(m=3, n=7))
        code = family.get_code()

        self.assertEqual(code.shape, (14, 6))
        self.assertEqual(family.number_of_workers_required(), 7)
        self.assertEqual(family.number_of_encoded_sets(), 6)

    def test_finite_field_codes_are_explicitly_unimplemented(self) -> None:
        """Finite-field families are not silently approximated by another code."""
        for code_type in ("BCH", "Reed Solomon"):
            with self.assertRaisesRegex(NotImplementedError, "finite-field"):
                create_code_family(code_type).generate_code(CodeParameters(m=26, n=31))

    def test_rejects_rng_for_deterministic_codes(self) -> None:
        """Random generators are rejected when a code ignores randomness."""
        with self.assertRaisesRegex(ValueError, "rng is not used"):
            create_code_family("LPF", rng=np.random.default_rng(1))


if __name__ == "__main__":
    unittest.main()
