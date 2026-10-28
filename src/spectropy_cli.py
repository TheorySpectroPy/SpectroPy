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
        if name == "derivatives":
            command_parser.add_argument("--difference", choices=("central", "forward"), default=None, help="Override input difference_scheme (default: central)")
            command_parser.add_argument("--equilibrium-file", "--equilibrium-vasprun", dest="equilibrium_vasprun", help="Optional external reference matching dielectric_source; defaults to displacements/ra_pos_equilibrium/OUTCAR or vasprun.xml")
            command_parser.add_argument("--saved-basis", action="store_true", help="Optional saved reconstruction lookup; manual reconstruction is the default")

    check_parser = subparsers.add_parser(
        "check", help="Export atomic-projection or direct mode-covariance residuals.",
        description="Check atomic derivative symmetry or mode-specific irrep covariance.",
    )
    add_workdir(check_parser)
    check_parser.add_argument("--kind", choices=("atomic", "mode"), default="atomic",
                              help="atomic: derivative projection; mode: direct irrep covariance (requires irreps.yaml).")
    check_parser.add_argument("--subspace-tolerance", type=float, default=1e-5,
                              help="Maximum mode-subspace leakage for --kind mode (default: 1e-5).")
    check_parser.add_argument("--energies", nargs="+", type=float, help="Laser energies in eV (default: input).")
    check_parser.add_argument("--output", default="check", help="Output directory relative to workdir (default: check).")
    check_parser.add_argument("--epsilon-fraction", type=float, default=1e-12, help="Atomic-projection check only: relative norm regularizer (default: 1e-12).")
    check_parser.add_argument("--plot", action="store_true", help="Also create residual PNG plots; requires Matplotlib.")
    check_parser.add_argument("--style", help="Matplotlib style file for --plot.")

    angular_parser = subparsers.add_parser(
        "angular", help="Export mode-specific parallel and crossed polarization-angle scans.",
        description="Generate angular intensity data from mode tensors independently of symmetry checks.",
    )
    add_workdir(angular_parser)
    angular_parser.add_argument("--modes", nargs="+", type=int, required=True, help="1-based mode numbers to scan.")
    angular_parser.add_argument("--energies", nargs="+", type=float, help="Laser energies in eV (default: input).")
    angular_parser.add_argument("--output", default="angular", help="Output directory relative to workdir (default: angular).")
    angular_parser.add_argument("--angle-step", type=float, default=1.0, help="Angular sampling step in degrees (default: 1).")
    angular_parser.add_argument("--surface-normal", choices=("x", "y", "z"), help="Scan normal (default: input, or z).")
    angular_parser.add_argument("--plot", action="store_true", help="Also create polar PNG plots; requires Matplotlib.")
    angular_parser.add_argument("--style", help="Matplotlib style file for --plot.")

    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    if not os.path.isdir(args.workdir):
        subparsers.choices[args.command].error(f"work directory does not exist: {args.workdir}")

    if args.command == "derivatives":
        def derive():
            importlib.import_module("spectropy_derivatives").run_derivatives(
                use_saved_basis=args.saved_basis, difference_scheme=args.difference,
                equilibrium_vasprun=args.equilibrium_vasprun)
        _run_in_directory(args.workdir, derive)
    elif args.command in ("check", "angular"):
        style = os.path.abspath(args.style) if args.style else None
        def analyze() -> None:
            module = importlib.import_module("spectropy_checks")
            common = dict(energies=args.energies, output=args.output, plot=args.plot, style=style)
            if args.command == "check":
                if args.kind == "mode":
                    module.run_mode_check(subspace_tolerance=args.subspace_tolerance, **common)
                else:
                    module.run_check(epsilon_fraction=args.epsilon_fraction, **common)
            else:
                module.run_angular(args.modes, angle_step=args.angle_step,
                                   surface_normal=args.surface_normal, **common)
        _run_in_directory(args.workdir, analyze)
    else:
        commands = _commands(getattr(args, "mode", "full"))
        _run_in_directory(args.workdir, commands[args.command][1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
