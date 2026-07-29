"""Command-line interface for SpectroPy's existing calculation stages."""

from __future__ import annotations

import argparse
import importlib
import os
from collections.abc import Callable


def _run_in_directory(workdir: str, operation: Callable[[], None]) -> None:
    """Run an existing stage in *workdir*, restoring the caller's CWD."""
    previous = os.getcwd()
    try:
        os.chdir(workdir)
        operation()
    finally:
        os.chdir(previous)


def _commands(mode: str) -> dict[str, tuple[str, Callable[[], None]]]:
    """Return the available CLI stages without importing their dependencies."""
    def stage(module: str, function: str) -> Callable[[], None]:
        def run() -> None:
            getattr(importlib.import_module(module), function)()
        return run

    def stages(*stage_specs: tuple[str, str]) -> Callable[[], None]:
        def run() -> None:
            for module, function in stage_specs:
                getattr(importlib.import_module(module), function)()
        return run

    def derivatives() -> None:
        getattr(importlib.import_module("process_symmetry"), "run_mapping")()
        getattr(
            importlib.import_module("calculate_dielectric_derivatives"),
            "run_generate_derivatives_for_input",
        )()

    def spectrum() -> None:
        getattr(importlib.import_module("calculate_spectrum"), "run_raman_tensors_for_input")()

    if mode == "full":
        prepare = stages(
            ("create_displacements", "run_displacements"),
            ("prepare_vasp_inputs", "run_generate_displacements"),
        )
        prepare_description = (
            "Generate the full +/- Cartesian displacement set and VASP-ready directories."
        )
    elif mode == "atoms":
        prepare = stages(
            ("generate_atom_displacements", "run_generate_atoms"),
            ("prepare_vasp_inputs", "run_generate_displacements"),
        )
        prepare_description = (
            "Generate all +/- Cartesian displacements for symmetry-inequivalent atoms."
        )
    else:
        prepare = stage("generate_minimal_displacements", "run_generate")
        prepare_description = "Generate Phonopy symmetry-reduced displaced POSCAR directories."

    return {
        "displacements": (
            prepare_description,
            prepare,
        ),
        "derivatives": (
            "Process symmetry and evaluate derivatives at every laser energy in input.",
            derivatives,
        ),
        "spectrum": (
            "Build Raman tensors and intensities for every laser energy in input.",
            spectrum,
        ),
        "plot": (
            "Interactively broaden and plot generated Raman intensities.",
            stage("generate_raman_plots", "run_automation"),
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="spectropy",
        description=(
            "Prepare finite-displacement Raman calculations and analyze their "
            "dielectric-tensor results. DFT calculations are run externally."
        ),
    )
    parser.add_argument("command", nargs="?", help="SpectroPy stage to run")
    parser.add_argument(
        "--workdir", "-C", default=".", metavar="PATH",
        help="Directory containing the stage's input files (default: current directory).",
    )
    parser.add_argument(
        "--mode", choices=("full", "atoms", "minimal"), default="full",
        help="Displacement strategy for `displacements` (default: full).",
    )
    parser.add_argument("--version", action="version", version="spectropy 0.1.0")
    args = parser.parse_args(argv)

    commands = _commands(args.mode)
    if args.command is None:
        parser.print_help()
        print("\nCommands:")
        for name, (description, _) in commands.items():
            print(f"  {name:19} {description}")
        return 0
    if args.command not in commands:
        parser.error(f"unknown command {args.command!r}; choose from {', '.join(commands)}")
    if not os.path.isdir(args.workdir):
        parser.error(f"work directory does not exist: {args.workdir}")

    _run_in_directory(args.workdir, commands[args.command][1])
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
