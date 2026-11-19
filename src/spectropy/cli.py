from __future__ import annotations

import argparse
import os
from collections.abc import Callable
from . import __version__
from .config import read_settings, resolve_settings


def _run_in_directory(workdir: str, operation: Callable[[], None]) -> None:
    previous = os.getcwd()
    try:
        os.chdir(workdir)
        operation()
    finally:
        os.chdir(previous)


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
    displacement_parser.add_argument("--phonopy", dest="phonopy_file", metavar="YAML", default=None,
                                     help="Import a Phonopy single-atom displacement set instead of selecting directions; overrides displacements.phonopy_file.")
    for name, help_text in (
        ("symmetry", "Write atom-to-atom symmetry-operation matrices."),
        ("tensors", "Calculate final mode Raman tensors from displaced dielectric responses."),
        ("spectrum", "Calculate polarized intensities from final Raman tensors."),
        ("plot", "Broaden and plot Raman intensities."),
    ):
        command_parser = subparsers.add_parser(name, help=help_text, description=help_text)
        add_workdir(command_parser)
        if name == "tensors":
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

    def execute():
        settings = read_settings()
        if args.command == "displacements":
            from .displacements import run_displacements
            resolved = resolve_settings(settings, symmetry=args.symmetry, displacement_basis=args.basis,
                                        phonopy_file=args.phonopy_file)
            run_displacements(settings=resolved)
        elif args.command == "tensors":
            from .tensors import run_tensors
            resolved = resolve_settings(settings, equilibrium_file=args.equilibrium_vasprun,
                                        saved_basis=args.saved_basis)
            run_tensors(settings=resolved)
        elif args.command == "mode-displacements":
            from .modes import generate
            resolved = resolve_settings(settings, mode_q=args.q, mode_indices=args.modes,
                                        mode_displacements_directory=args.output,
                                        mode_minimum_frequency=args.minimum_frequency)
            generate(settings=resolved)
        elif args.command == "mode-tensors":
            from .modes import reconstruct
            resolved = resolve_settings(settings, mode_displacements_directory=args.directory,
                                        mode_results_directory=args.output, laser_energies=args.energies)
            reconstruct(settings=resolved)
        elif args.command == "spectrum":
            from .spectrum import run_spectrum
            run_spectrum(settings=settings)
        elif args.command == "plot":
            from .plotting import run_automation
            run_automation(settings=settings)
        elif args.command == "symmetry":
            from .symmetry import run_mapping
            run_mapping()

    _run_in_directory(args.workdir, execute)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
