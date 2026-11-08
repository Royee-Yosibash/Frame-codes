"""Unit tests for shared code parameters and performance statistics."""

import unittest

import numpy as np
from frame_codes.code_families.factory import create_code_family
from frame_codes.coding_scheme.codes.code_parameters import CodeParameters


class TestCodeParameters(unittest.TestCase):
    """Test shared dimensions and normalization settings."""

    def test_gamma_is_derived_from_dimensions(self) -> None:
        """The aspect ratio equals m divided by n."""
        parameters = CodeParameters(m=2, n=4)

        self.assertEqual(parameters.gamma, 0.5)

    def test_code_family_uses_the_same_parameters(self) -> None:
        """The parameters object is passed directly to code construction."""
        parameters = CodeParameters(m=2, n=4, norm_dim="column")
        family = create_code_family("LPF")

        family.generate_code(parameters)
        generator = family.get_code()

        self.assertEqual(generator.shape, (4, 2))
        np.testing.assert_allclose(np.linalg.norm(generator, axis=1), 1)

    def test_none_normalization_is_default(self) -> None:
        """No axis normalization is selected by default."""
        self.assertIsNone(CodeParameters(m=2, n=4).norm_dim)

    def test_normalization_setting_is_validated(self) -> None:
        """Unsupported normalization axes fail at parameter construction."""
        with self.assertRaisesRegex(ValueError, "norm_dim"):
            CodeParameters(m=2, n=4, norm_dim="diagonal")

    def test_rejects_invalid_dimensions(self) -> None:
        """Code parameters require positive dimensions and m no greater than n."""
        with self.assertRaisesRegex(ValueError, "m must not exceed n"):
            CodeParameters(m=5, n=4)

if __name__ == "__main__":
    unittest.main()
