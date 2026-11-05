"""Tests for Hydra construction of the experiment workspace."""

import unittest
from pathlib import Path

from hydra import compose, initialize_config_dir
from hydra.utils import instantiate

from frame_codes.coding_scheme.codes.code_parameters import CodeParameters
from frame_codes.coding_scheme.codes.families.simple_band_circulant import (
    SimpleBandCirculantCode,
)
from frame_codes.coding_scheme.decoding.decoders.pseudo_inverse import (
    PseudoInverseDecoder,
)
from frame_codes.coding_scheme.encoding.encoder import Encoder
from frame_codes.coding_scheme.worker_nodes.worker_dispatcher import WorkerDispatcher


class TestHydraExperimentInit(unittest.TestCase):
    """Verify configured components populate the workspace."""

    def test_init_config_constructs_shared_components(self) -> None:
        """Hydra builds the code, encoder, dispatcher, and decoder objects."""
        config_directory = Path(__file__).resolve().parents[2] / "experiments"
        with initialize_config_dir(version_base=None, config_dir=str(config_directory)):
            config = compose(config_name="experiment_init")

        workspace = instantiate(config.init)

        self.assertIsInstance(workspace, dict)
        self.assertIsInstance(workspace["code_family"], SimpleBandCirculantCode)
        self.assertIsInstance(workspace["code_parameters"], CodeParameters)
        self.assertEqual(workspace["code_parameters"].m, 3)
        self.assertEqual(workspace["code_parameters"].n, 6)
        self.assertIsInstance(workspace["encoder"], Encoder)
        self.assertIsInstance(workspace["dispatcher"], WorkerDispatcher)
        self.assertIsInstance(workspace["decoder"], PseudoInverseDecoder)
        self.assertTrue(callable(workspace["worker_operation"]))


if __name__ == "__main__":
    unittest.main()
