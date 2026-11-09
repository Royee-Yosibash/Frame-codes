"""Run a workspace pipeline described by a Hydra configuration."""

from typing import Any

from hydra.utils import instantiate
from omegaconf import DictConfig

from frame_codes.coding_scheme.hydra_pipeline.workspace_pipeline import (
    WorkspacePipeline,
)


def run_experiment(config: DictConfig) -> dict[str, Any]:
    """Instantiate the workspace and execute its configured pipeline.

    Args:
        config: Hydra configuration with an `init` workspace and pipeline.

    Returns:
        The populated workspace after all configured pipeline steps complete.
    """
    workspace = instantiate(config.init)
    workspace["config"] = config
    WorkspacePipeline(workspace).run(config.pipeline)
    return workspace
