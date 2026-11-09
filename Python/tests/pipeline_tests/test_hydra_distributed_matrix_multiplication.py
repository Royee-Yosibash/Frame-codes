"""End-to-end Hydra pipeline tests for distributed matrix multiplication."""

import unittest
from itertools import combinations
from pathlib import Path

import numpy as np
from hydra import compose, initialize_config_dir
from omegaconf import DictConfig, OmegaConf

from frame_codes.coding_scheme.hydra_pipeline.runner import run_experiment


class TestHydraDistributedMatrixMultiplication(unittest.TestCase):
    """Verify configured pipeline execution for matrix and vector RHS inputs."""

    @staticmethod
    def _compose_config() -> DictConfig:
        """Compose the distributed matrix multiplication Hydra configuration.

        Returns:
            The composed experiment configuration.
        """
        config_directory = Path(__file__).resolve().parents[2] / "experiments"
        with initialize_config_dir(version_base=None, config_dir=str(config_directory)):
            return compose(config_name="simple_band_circulant_matrix_multiplication")

    def test_matrix_rhs_runs_full_distributed_pipeline(self) -> None:
        """Every two-worker straggler pair decodes the matrix-RHS product."""
        config = self._compose_config()

        workspace = run_experiment(config)

        expected = np.asarray(config.init.lhs) @ np.asarray(config.init.rhs)
        worker_count = config.init.code_parameters.n
        expected_stragglers = list(combinations(range(worker_count), 2))
        trials = workspace["sweep_results"]

        self.assertEqual(len(trials), len(expected_stragglers))
        self.assertEqual(
            [trial["straggler_worker_ids"] for trial in trials],
            expected_stragglers,
        )
        for trial in trials:
            with self.subTest(stragglers=trial["straggler_worker_ids"]):
                expected_responders = tuple(
                    worker_id
                    for worker_id in range(worker_count)
                    if worker_id not in trial["straggler_worker_ids"]
                )
                self.assertEqual(trial["responding_worker_ids"], expected_responders)
                np.testing.assert_allclose(trial["result"], expected, atol=1e-12)

    def test_vector_rhs_runs_full_distributed_pipeline(self) -> None:
        """Every two-worker straggler pair decodes the vector-RHS product."""
        config = self._compose_config()
        OmegaConf.update(config, "init.rhs", [2.0, -1.0, 3.0])

        workspace = run_experiment(config)

        expected = np.asarray(config.init.lhs) @ np.asarray(config.init.rhs)
        worker_count = config.init.code_parameters.n
        self.assertEqual(
            len(workspace["sweep_results"]),
            len(list(combinations(range(worker_count), 2))),
        )
        for trial in workspace["sweep_results"]:
            with self.subTest(stragglers=trial["straggler_worker_ids"]):
                np.testing.assert_allclose(trial["result"], expected, atol=1e-12)
                self.assertEqual(trial["result"].ndim, 1)


if __name__ == "__main__":
    unittest.main()
