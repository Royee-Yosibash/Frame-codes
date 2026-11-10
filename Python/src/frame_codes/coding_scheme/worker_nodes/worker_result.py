"""Result returned by a worker node."""

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class WorkerResult:
    """Store a worker's outputs under its stable identifier."""

    worker_id: int
    values: list[Any]
