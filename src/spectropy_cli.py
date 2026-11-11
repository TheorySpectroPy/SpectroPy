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


def _commands(mode: str, basis: str = "cartesian") -> dict[str, tuple[str, Callable[[], None]]]:
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
        "--mode", choices=("full", "atoms", "minimal"), default="full",
        help=(
            "full: all atoms and +/- basis axes; atoms: all axes for "
            "symmetry-inequivalent atoms; minimal: site-symmetry-reduced "
            "set including mixed directions in the selected basis."
        ),
    )
    displacement_parser.add_argument("--basis", choices=("cartesian", "fractional"), default="cartesian",
                                     help="Basis for displacement candidates; lengths are always in Angstrom.")
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
            command_parser.add_argument("--equilibrium-file", "--equilibrium-vasprun", dest="equilibrium_vasprun", help="Optional external reference in XML format; defaults to displacements/ra_pos_equilibrium/vasprun.xml")
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

    mode_displacements = subparsers.add_parser("mode-displacements", help="Generate explicit +/- displacements for every Raman-active mode.")
    add_workdir(mode_displacements)
    mode_displacements.add_argument("--modes", nargs="+", type=int, help="Explicit 1-based indices overriding activity/cutoff selection; --reduce expands partner groups.")
    mode_displacements.add_argument("--q", type=float, default=.01, help="Normal coordinate in Angstrom sqrt(amu), matching Liangbo's raman_dis_mode.")
    reduction = mode_displacements.add_mutually_exclusive_group()
    reduction.add_argument("--reduce", action="store_true", help="Experimental opt-in: reduce degenerate partners and reuse symmetry-related signs.")
    reduction.add_argument("--no-reduction", action="store_true", help="Explicit original-mode +/- pairs (the default).")
    mode_displacements.add_argument("--keep-both-signs", action="store_true", help="Reduce partners but compute both signs explicitly.")
    mode_displacements.add_argument("--output", default="mode_displacements")
    mode_displacements.add_argument("--minimum-frequency", type=float, default=1.0, help="Default automatic selection excludes frequencies below this cm^-1.")
    mode_tensors = subparsers.add_parser("mode-tensors", help="Reconstruct mode tensors directly from collective displaced dielectric tensors.")
    add_workdir(mode_tensors)
    mode_tensors.add_argument("--directory", default="mode_displacements")
    mode_tensors.add_argument("--output", default="mode_results")
    mode_tensors.add_argument("--source", choices=("vasprun.xml",), default="vasprun.xml")
    mode_tensors.add_argument("--energies", nargs="+", type=float)
    mode_tensors.add_argument("--denoise", action="store_true", help="Project complete mode partner groups onto Neumann covariance; preserve raw tensors (VASP only).")

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
                                reduce=args.reduce, reuse_sign=not args.keep_both_signs,
                                output=args.output, minimum_frequency=args.minimum_frequency)
            else:
                module.reconstruct(args.directory, args.source, args.energies, args.output, denoise=args.denoise)
        _run_in_directory(args.workdir, collective)
    elif args.command == "derivatives":
        def derive():
            importlib.import_module("spectropy_derivatives").run_derivatives(
                use_saved_basis=args.saved_basis, difference_scheme=args.difference,
                equilibrium_vasprun=args.equilibrium_vasprun)
        _run_in_directory(args.workdir, derive)
    elif args.command == "check":
        style = os.path.abspath(args.style) if args.style else None
        def analyze() -> None:
            module = importlib.import_module("spectropy_checks")
            common = dict(energies=args.energies, output=args.output, plot=args.plot, style=style)
            if args.kind == "mode":
                module.run_mode_check(subspace_tolerance=args.subspace_tolerance, **common)
            else:
                module.run_check(epsilon_fraction=args.epsilon_fraction, **common)
        _run_in_directory(args.workdir, analyze)
    else:
        commands = _commands(getattr(args, "mode", "full"), getattr(args, "basis", "cartesian"))
        _run_in_directory(args.workdir, commands[args.command][1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
