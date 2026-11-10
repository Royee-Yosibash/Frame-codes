"""Execute configured calls and workspace-bound loops."""

from collections.abc import Mapping, Sequence
from typing import Any

from hydra.utils import get_method
from omegaconf import DictConfig


class WorkspacePipeline:
    """Execute configured calls against a shared workspace dictionary."""

    def __init__(self, workspace: dict[str, Any]) -> None:
        """Store the workspace shared by all configured steps.

        Args:
            workspace: Mutable state containing experiment inputs and objects.
        """
        self.workspace = workspace

    def run(self, steps: Sequence[DictConfig]) -> None:
        """Execute call steps and workspace loops in order.

        Args:
            steps: Ordered pipeline call descriptions.
        """
        for step in steps:
            if "for_each" in step:
                loop = step.for_each
                for value in self._resolve_value(loop["values"]):
                    if "as" in loop:
                        self.workspace[loop["as"]] = value
                    self.run(loop.steps)
            elif "run" in step:
                self.run(self._resolve_value(step["run"]))
            else:
                self._run_step(step)

    def _run_step(self, step: DictConfig) -> None:
        """Resolve and call one operation, optionally naming its result.

        Args:
            step: A call description, optionally nested under the workspace key
                that should hold its result.
        """
        if "call" in step:
            output_name = None
            call = step
        else:
            output_name, call = next(iter(step.items()))

        operation = self._resolve_operation(call.call)
        arguments = [self._resolve_value(value) for value in call.get("args", [])]
        keyword_arguments = {
            name: self._resolve_value(value)
            for name, value in call.get("kwargs", {}).items()
        }
        result = operation(*arguments, **keyword_arguments)
        if output_name is not None:
            self.workspace[output_name] = result

    def _resolve_operation(self, call: str) -> Any:
        """Resolve a workspace method or importable callable.

        Args:
            call: Dotted workspace path or Python import path.

        Returns:
            The callable referenced by the pipeline step.
        """
        first, *remaining = call.split(".")
        if first not in self.workspace:
            return get_method(call)
        operation = self.workspace[first]
        for attribute in remaining:
            operation = getattr(operation, attribute)
        return operation

    def _resolve_value(self, value: Any) -> Any:
        """Resolve workspace references and leave literal values unchanged.

        Args:
            value: Literal input, a string prefixed by `@` referencing the
                workspace, or a sequence of values to resolve recursively.

        Returns:
            The resolved literal, sequence, or workspace value.
        """
        if isinstance(value, str) and value.startswith("@"):
            return self._resolve_workspace_path(value[1:])
        if isinstance(value, Mapping):
            return {key: self._resolve_value(item) for key, item in value.items()}
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
            return [self._resolve_value(item) for item in value]
        return value

    def _resolve_workspace_path(self, path: str) -> Any:
        """Follow a dotted key/attribute path from the shared workspace.

        Args:
            path: Workspace key followed by optional nested keys or attributes.
                The `workspace` path refers to the full shared dictionary.

        Returns:
            The value found at the requested workspace path.
        """
        if path == "workspace":
            return self.workspace
        first, *remaining = path.split(".")
        value = self.workspace[first]
        for part in remaining:
            value = value[part] if isinstance(value, Mapping) else getattr(value, part)
        return value
