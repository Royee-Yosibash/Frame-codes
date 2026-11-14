"""Unit tests for shared code parameters and performance statistics."""

import unittest

import numpy as np
from codedcompute.codes.code_parameters import CodeParameters

from frame_codes.code_families.factory import create_code_family


class TestCodeParameters(unittest.TestCase):
    """Test shared dimensions and normalization settings."""

    def test_code_family_uses_the_same_parameters(self) -> None:
        """The parameters object is passed directly to code construction."""
        parameters = CodeParameters(m=2, n=4, norm_dim="column")
        family = create_code_family("LPF")

        family.generate_code(parameters)
        generator = family.get_code()

        self.assertEqual(generator.shape, (4, 2))
        np.testing.assert_allclose(np.linalg.norm(generator, axis=1), 1)


if __name__ == "__main__":
    unittest.main()
