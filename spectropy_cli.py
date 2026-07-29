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


def _commands() -> dict[str, tuple[str, Callable[[], None]]]:
    """Return the available CLI stages without importing their dependencies."""
    def stage(module: str, function: str) -> Callable[[], None]:
        def run() -> None:
            getattr(importlib.import_module(module), function)()
        return run

    return {
        "displacements": (
            "Create the legacy +/- Cartesian displacement list from CONTCAR.",
            stage("create_displacements", "run_displacements"),
        ),
        "prepare-inputs": (
            "Create displaced POSCAR directories from displacements.dat.",
            stage("prepare_vasp_inputs", "run_generate_displacements"),
        ),
        "minimal-displacements": (
            "Create phonopy symmetry-reduced displaced POSCAR directories.",
            stage("generate_minimal_displacements", "run_generate"),
        ),
        "symmetry": (
            "Process phonopy's symmetry file into atom mapping matrices.",
            stage("process_symmetry", "run_mapping"),
        ),
        "derivatives": (
            "Read dielectric outputs and write frequency-dependent derivatives.",
            stage("calculate_dielectric_derivatives", "run_generate_derivatives"),
        ),
        "static-derivatives": (
            "Read static dielectric outputs and write static derivatives.",
            stage("calculate_dielectric_derivatives_static", "run_calculate_dielectric_derivatives_static"),
        ),
        "spectrum": (
            "Build Placzek Raman tensors and intensities from band.yaml.",
            stage("calculate_spectrum", "run_raman_tensor"),
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
    parser.add_argument("--version", action="version", version="spectropy 0.1.0")
    args = parser.parse_args(argv)

    commands = _commands()
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
