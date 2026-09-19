from __future__ import annotations

import argparse
import importlib
import os
from collections.abc import Callable


def _run_in_directory(workdir: str, operation: Callable[[], None]) -> None:
    previous = os.getcwd()
    try:
        os.chdir(workdir)
        operation()
    finally:
        os.chdir(previous)


def _commands(mode: str) -> dict[str, tuple[str, Callable[[], None]]]:
    def stage(module: str, function: str) -> Callable[[], None]:
        def run() -> None:
            getattr(importlib.import_module(module), function)()
        return run

    def derivatives() -> None:
        getattr(importlib.import_module("spectropy_derivatives"), "run_derivatives")()

    def spectrum() -> None:
        getattr(importlib.import_module("calculate_spectrum"), "run_raman_tensors_for_input")()

    def displacements() -> None:
        getattr(importlib.import_module("spectropy_displacements"), "run_displacements")(mode)

    descriptions = {
        "full": "Generate the full +/- Cartesian displacement set and VASP-ready directories.",
        "atoms": "Generate all +/- Cartesian displacements for symmetry-inequivalent atoms.",
        "minimal": "Generate the Cartesian site-symmetry-minimal displacement set.",
        "fractional": "Generate Phonopy's fractional-basis (lattice-direction) minimal displacement set.",
    }

    return {
        "displacements": (
            descriptions[mode],
            displacements,
        ),
        "derivatives": (
            "Process symmetry and evaluate derivatives at every laser energy in input.",
            derivatives,
        ),
        "symmetry": (
            "Write the symmetry operations that map equivalent atoms.",
            stage("process_symmetry", "run_mapping"),
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
    parser.add_argument("--version", action="version", version="spectropy 0.1.0")
    subparsers = parser.add_subparsers(dest="command", title="commands", metavar="COMMAND")

    def add_workdir(command_parser: argparse.ArgumentParser) -> None:
        command_parser.add_argument(
            "--workdir", "-C", default=".", metavar="PATH",
            help="Directory containing this stage's input files (default: current directory).",
        )

    displacement_parser = subparsers.add_parser(
        "displacements",
        help="Generate displaced structures and DFT-ready directories.",
        description="Generate displacement structures; this does not run a DFT code.",
    )
    add_workdir(displacement_parser)
    displacement_parser.add_argument(
        "--mode", choices=("full", "atoms", "minimal", "fractional"), default="full",
        help=(
            "full: all atoms and +/- Cartesian axes; atoms: all axes for "
            "symmetry-inequivalent atoms; minimal: Cartesian site-symmetry-minimal "
            "set; fractional: Phonopy's lattice-basis (fractional) minimum set."
        ),
    )
    for name, help_text in (
        ("symmetry", "Write atom-to-atom symmetry-operation matrices."),
        ("derivatives", "Calculate dielectric derivatives."),
        ("spectrum", "Build Raman tensors and intensities from derivatives."),
        ("plot", "Broaden and plot Raman intensities."),
    ):
        command_parser = subparsers.add_parser(name, help=help_text, description=help_text)
        add_workdir(command_parser)

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if not os.path.isdir(args.workdir):
        subparsers.choices[args.command].error(f"work directory does not exist: {args.workdir}")

    commands = _commands(getattr(args, "mode", "full"))
    _run_in_directory(args.workdir, commands[args.command][1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
