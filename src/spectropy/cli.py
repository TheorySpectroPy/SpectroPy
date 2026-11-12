from __future__ import annotations

import argparse
import importlib
import os
from collections.abc import Callable
from . import __version__


def _run_in_directory(workdir: str, operation: Callable[[], None]) -> None:
    previous = os.getcwd()
    try:
        os.chdir(workdir)
        operation()
    finally:
        os.chdir(previous)


def _commands(symmetry: str | None, basis: str | None = None) -> dict[str, tuple[str, Callable[[], None]]]:
    def stage(module: str, function: str) -> Callable[[], None]:
        def run() -> None:
            getattr(importlib.import_module(module), function)()
        return run

    def spectrum() -> None:
        getattr(importlib.import_module('spectropy.spectrum'), "run_spectrum")()

    def displacements() -> None:
        getattr(importlib.import_module('spectropy.displacements'), "run_displacements")(symmetry, basis=basis)

    descriptions = {
        "off": "Generate all +/- basis-axis displacements on every atom.",
        "atoms": "Generate all +/- basis-axis displacements on atomic representatives.",
        "on": "Generate the site-symmetry-reduced displacement set.",
    }

    return {
        "displacements": (
            descriptions.get(symmetry, "Generate displacements using input.yaml settings."),
            displacements,
        ),
        "symmetry": (
            "Write the symmetry operations that map equivalent atoms.",
            stage('spectropy.symmetry', "run_mapping"),
        ),
        "spectrum": (
            "Calculate polarized intensities for every laser energy in input.yaml.",
            spectrum,
        ),
        "plot": (
            "Interactively broaden and plot generated Raman intensities.",
            stage('spectropy.plotting', "run_automation"),
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
    parser.add_argument("--version", action="version", version=f"spectropy {__version__}")
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
        "--symmetry", choices=("on", "off", "atoms"), default=None,
        help=(
            "Override input.yaml symmetry.reduction (default: off). "
            "off: all atoms and +/- basis axes; atoms: all axes for "
            "symmetry-inequivalent atoms; on: site-symmetry-reduced "
            "set including mixed directions in the selected basis."
        ),
    )
    displacement_parser.add_argument("--basis", choices=("cartesian", "fractional"), default=None,
                                     help="Override input.yaml displacements.basis (default: cartesian); lengths are in Angstrom.")
    for name, help_text in (
        ("symmetry", "Write atom-to-atom symmetry-operation matrices."),
        ("tensors", "Calculate final mode Raman tensors from displaced dielectric responses."),
        ("spectrum", "Calculate polarized intensities from final Raman tensors."),
        ("plot", "Broaden and plot Raman intensities."),
    ):
        command_parser = subparsers.add_parser(name, help=help_text, description=help_text)
        add_workdir(command_parser)
        if name == "tensors":
            command_parser.add_argument("--difference", choices=("central", "forward"), default=None, help="Override input.yaml displacements.difference (default: central)")
            command_parser.add_argument("--equilibrium-file", "--equilibrium-vasprun", dest="equilibrium_vasprun", help="Override input.yaml displacements.equilibrium_file; otherwise use the generated equilibrium XML")
            command_parser.add_argument("--saved-basis", action=argparse.BooleanOptionalAction, default=None,
                                        help="Override input.yaml symmetry.saved_basis (default: false)")

    mode_displacements = subparsers.add_parser("mode-displacements", help="Generate explicit +/- displacements for every Raman-active mode.")
    add_workdir(mode_displacements)
    mode_displacements.add_argument("--modes", nargs="+", type=int, help="Override input.yaml modes.indices with 1-based indices (default: automatic activity selection).")
    mode_displacements.add_argument("--q", type=float, default=None, help="Override input.yaml modes.q, in Angstrom sqrt(amu) (default: 0.01).")
    mode_displacements.add_argument("--output", default=None, help="Override input.yaml modes.displacements_directory (default: mode_displacements).")
    mode_displacements.add_argument("--minimum-frequency", type=float, default=None, help="Override input.yaml modes.minimum_frequency; automatic selection requires frequencies above this cm^-1 (default: 1.0).")
    mode_tensors = subparsers.add_parser("mode-tensors", help="Reconstruct mode tensors directly from collective displaced dielectric tensors.")
    add_workdir(mode_tensors)
    mode_tensors.add_argument("--directory", default=None, help="Override input.yaml modes.displacements_directory (default: mode_displacements).")
    mode_tensors.add_argument("--output", default=None, help="Override input.yaml modes.results_directory (default: mode_results).")
    mode_tensors.add_argument("--energies", nargs="+", type=float, help="Override input.yaml spectrum.energies (default: 0.00 eV).")

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if not os.path.isdir(args.workdir):
        subparsers.choices[args.command].error(f"work directory does not exist: {args.workdir}")

    if args.command in {"mode-displacements", "mode-tensors"}:
        def collective():
            module = importlib.import_module('spectropy.modes')
            if args.command == "mode-displacements":
                module.generate(q=args.q, modes=args.modes,
                                output=args.output, minimum_frequency=args.minimum_frequency)
            else:
                module.reconstruct(directory=args.directory, energies=args.energies, output=args.output)
        _run_in_directory(args.workdir, collective)
    elif args.command == "tensors":
        def derive():
            importlib.import_module('spectropy.tensors').run_tensors(
                use_saved_basis=args.saved_basis, difference_scheme=args.difference,
                equilibrium_vasprun=args.equilibrium_vasprun)
        _run_in_directory(args.workdir, derive)

    else:
        commands = _commands(getattr(args, "symmetry", None), getattr(args, "basis", None))
        _run_in_directory(args.workdir, commands[args.command][1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
