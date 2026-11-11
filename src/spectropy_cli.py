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


def _commands(mode: str | None, basis: str | None = None) -> dict[str, tuple[str, Callable[[], None]]]:
    def stage(module: str, function: str) -> Callable[[], None]:
        def run() -> None:
            getattr(importlib.import_module(module), function)()
        return run

    def derivatives() -> None:
        getattr(importlib.import_module("spectropy_derivatives"), "run_derivatives")()

    def spectrum() -> None:
        getattr(importlib.import_module("calculate_spectrum"), "run_raman_tensors_for_input")()

    def displacements() -> None:
        getattr(importlib.import_module("spectropy_displacements"), "run_displacements")(mode, basis=basis)

    descriptions = {
        "full": "Generate all +/- basis-axis displacements on every atom.",
        "atoms": "Generate all +/- basis-axis displacements on atomic representatives.",
        "minimal": "Generate the site-symmetry-reduced displacement set.",
    }

    return {
        "displacements": (
            descriptions.get(mode, "Generate displacements using input settings."),
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
        "--mode", choices=("full", "atoms", "minimal"), default=None,
        help=(
            "Override input displacement_mode (default: full). "
            "full: all atoms and +/- basis axes; atoms: all axes for "
            "symmetry-inequivalent atoms; minimal: site-symmetry-reduced "
            "set including mixed directions in the selected basis."
        ),
    )
    displacement_parser.add_argument("--basis", choices=("cartesian", "fractional"), default=None,
                                     help="Override input displacement_basis (default: cartesian); lengths are in Angstrom.")
    for name, help_text in (
        ("symmetry", "Write atom-to-atom symmetry-operation matrices."),
        ("derivatives", "Calculate dielectric derivatives."),
        ("spectrum", "Build Raman tensors and intensities from derivatives."),
        ("plot", "Broaden and plot Raman intensities."),
    ):
        command_parser = subparsers.add_parser(name, help=help_text, description=help_text)
        add_workdir(command_parser)
        if name == "derivatives":
            command_parser.add_argument("--difference", choices=("central", "forward"), default=None, help="Override input difference_scheme (default: central)")
            command_parser.add_argument("--equilibrium-file", "--equilibrium-vasprun", dest="equilibrium_vasprun", help="Override input equilibrium_file; otherwise use the generated equilibrium XML")
            command_parser.add_argument("--saved-basis", action=argparse.BooleanOptionalAction, default=None,
                                        help="Override input saved_basis (default: false)")

    mode_displacements = subparsers.add_parser("mode-displacements", help="Generate explicit +/- displacements for every Raman-active mode.")
    add_workdir(mode_displacements)
    mode_displacements.add_argument("--modes", nargs="+", type=int, help="Override input mode_indices with 1-based indices (default: automatic activity selection).")
    mode_displacements.add_argument("--q", type=float, default=None, help="Override input mode_q, in Angstrom sqrt(amu) (default: 0.01).")
    mode_displacements.add_argument("--output", default=None, help="Override input mode_displacements_directory (default: mode_displacements).")
    mode_displacements.add_argument("--minimum-frequency", type=float, default=None, help="Override input mode_minimum_frequency; automatic selection requires frequencies above this cm^-1 (default: 1.0).")
    mode_tensors = subparsers.add_parser("mode-tensors", help="Reconstruct mode tensors directly from collective displaced dielectric tensors.")
    add_workdir(mode_tensors)
    mode_tensors.add_argument("--directory", default=None, help="Override input mode_displacements_directory (default: mode_displacements).")
    mode_tensors.add_argument("--output", default=None, help="Override input mode_results_directory (default: mode_results).")
    mode_tensors.add_argument("--source", choices=("vasprun.xml",), default=None, help="Override input dielectric_source (only vasprun.xml).")
    mode_tensors.add_argument("--energies", nargs="+", type=float, help="Override input laser_energies (default: 0.00 eV).")

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if not os.path.isdir(args.workdir):
        subparsers.choices[args.command].error(f"work directory does not exist: {args.workdir}")

    if args.command in {"mode-displacements", "mode-tensors"}:
        def collective():
            module = importlib.import_module("spectropy_modes")
            if args.command == "mode-displacements":
                module.generate(q=args.q, modes=args.modes,
                                output=args.output, minimum_frequency=args.minimum_frequency)
            else:
                module.reconstruct(args.directory, args.source, args.energies, args.output)
        _run_in_directory(args.workdir, collective)
    elif args.command == "derivatives":
        def derive():
            importlib.import_module("spectropy_derivatives").run_derivatives(
                use_saved_basis=args.saved_basis, difference_scheme=args.difference,
                equilibrium_vasprun=args.equilibrium_vasprun)
        _run_in_directory(args.workdir, derive)

    else:
        commands = _commands(getattr(args, "mode", None), getattr(args, "basis", None))
        _run_in_directory(args.workdir, commands[args.command][1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
