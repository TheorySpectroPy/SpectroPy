from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import math


_SETTING_RE = re.compile(r"^([a-z_]+)\s*(?::|=|\s)\s*(.+)$", re.I)


@dataclass(frozen=True)
class Settings:
    incident_polarization: tuple[float, float, float] | None = None
    scattered_polarization: tuple[float, float, float] | None = None
    surface_normal: str | None = None
    laser_energies: tuple[float, ...] = (0.0,)
    displacement_amplitude: tuple[float, float, float] | None = None
    broadening_hwhm: float = 1.0
    broadening_type: str = "lorentzian"
    polarization: str = "specific"
    frequency_cutoff_cm1: float | None = None
    difference_scheme: str = "central"
    dielectric_source: str = "vasprun.xml"
    site_symmetry_projection: bool = False
    displacement_mode: str = "full"
    displacement_basis: str = "cartesian"
    equilibrium_file: str | None = None
    saved_basis: bool = False
    mode_q: float = 0.01
    mode_indices: tuple[int, ...] | None = None
    mode_minimum_frequency: float = 1.0
    mode_displacements_directory: str = "mode_displacements"
    mode_results_directory: str = "mode_results"

    @property
    def broadening_fwhm(self):
        """Full width for internal line-shape routines; public settings use HWHM."""
        return 2 * self.broadening_hwhm


def _content_lines(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [line.split("!", 1)[0].split("#", 1)[0].strip() for line in path.read_text().splitlines()]


def read_settings(path: str | Path = "input") -> Settings:
    geometry: list[list[str]] = []
    values: dict[str, str] = {}
    for line in _content_lines(Path(path)):
        if not line:
            continue
        match = _SETTING_RE.match(line)
        if match and match.group(1).lower() in {
            "laser_energy", "laser_energies", "displacement_amplitude",
            "broadening_hwhm", "broadening_fwhm", "broadening_type", "polarization", "cutoff",
            "frequency_cutoff", "difference_scheme", "dielectric_source", "site_symmetry_projection",
            "displacement_mode", "displacement_basis", "equilibrium_file", "saved_basis",
            "mode_q", "mode_indices", "mode_minimum_frequency",
            "mode_displacements_directory", "mode_results_directory",
            "incident_polarization", "scattered_polarization", "surface_normal",
        }:
            values[match.group(1).lower()] = match.group(2)
        elif re.match(r"^[a-z_]+\s*[:=]", line, re.I):
            raise ValueError(f"Unknown input setting: {line.split(':', 1)[0].split('=', 1)[0].strip()}")
        else:
            geometry.append(line.split())

    geometry_names = {'incident_polarization', 'scattered_polarization', 'surface_normal'}
    if geometry and geometry_names.intersection(values):
        raise ValueError('Do not mix named geometry settings with positional geometry lines')

    def vector(name: str, index: int) -> tuple[float, float, float] | None:
        parts = values[name].replace(',', ' ').split() if name in values else (geometry[index] if len(geometry) > index else None)
        if parts is None:
            return None
        try:
            result = tuple(float(value) for value in parts)
        except ValueError as error:
            raise ValueError(f'{name} must contain three finite numbers') from error
        if len(result) != 3 or not all(math.isfinite(value) for value in result):
            raise ValueError(f'{name} must contain three finite numbers')
        return result

    surface_normal = values.get('surface_normal', geometry[2][0] if len(geometry) >= 3 else None)
    if surface_normal is not None:
        surface_normal = surface_normal.lower()
        if surface_normal not in {'x', 'y', 'z'}:
            raise ValueError('surface_normal must be x, y, or z')

    energies_text = values.get("laser_energies", values.get("laser_energy", "0.00"))
    try:
        energies = tuple(dict.fromkeys(float(value) for value in energies_text.replace(",", " ").split()))
    except ValueError as error:
        raise ValueError(f"Invalid laser energy setting: {energies_text}") from error
    if not energies:
        raise ValueError("laser_energies must contain at least one value")

    amplitude = values.get("displacement_amplitude")
    if amplitude is not None:
        try:
            components = [float(component) for component in amplitude.replace(",", " ").split()]
        except ValueError as error:
            raise ValueError(f"Invalid displacement_amplitude setting: {amplitude}") from error
        if len(components) == 1:
            components *= 3
        if len(components) != 3:
            raise ValueError("displacement_amplitude takes one value or one per Cartesian axis")
        if min(components) <= 0:
            raise ValueError("displacement_amplitude must be positive")
        displacement_amplitude = tuple(components)
    else:
        displacement_amplitude = None

    if 'broadening_hwhm' in values and 'broadening_fwhm' in values:
        raise ValueError('Specify broadening_hwhm or legacy broadening_fwhm, not both')
    try:
        hwhm = (float(values['broadening_fwhm'].split()[0]) / 2
                if 'broadening_fwhm' in values else
                float(values.get('broadening_hwhm', '1.0').split()[0]))
    except ValueError as error:
        raise ValueError("broadening_hwhm must be numeric") from error
    if not math.isfinite(hwhm) or hwhm <= 0:
        raise ValueError("broadening_hwhm must be finite and positive")
    broadening_type = values.get("broadening_type", "lorentzian").lower()
    if broadening_type in {"l", "lorentzian"}:
        broadening_type = "lorentzian"
    elif broadening_type in {"g", "gaussian"}:
        broadening_type = "gaussian"
    else:
        raise ValueError("broadening_type must be lorentzian or gaussian")
    polarization = values.get("polarization", "specific").lower()
    if polarization not in {"specific", "average"}:
        raise ValueError("polarization must be specific or average")

    cutoff = values.get("cutoff", values.get("frequency_cutoff"))
    if cutoff is not None:
        try:
            frequency_cutoff_cm1 = float(cutoff.split()[0])
        except ValueError as error:
            raise ValueError("cutoff must be numeric") from error
        if frequency_cutoff_cm1 < 0:
            raise ValueError("cutoff must be non-negative")
    else:
        frequency_cutoff_cm1 = None

    difference_scheme = values.get('difference_scheme', 'central').lower()
    if difference_scheme not in {'central', 'forward'}:
        raise ValueError('difference_scheme must be central or forward')
    dielectric_source = values.get('dielectric_source', 'vasprun.xml').lower()
    if dielectric_source != 'vasprun.xml':
        raise ValueError('dielectric_source must be vasprun.xml')
    projection = values.get('site_symmetry_projection', 'false').lower()
    if projection not in {'true', 'false'}:
        raise ValueError('site_symmetry_projection must be true or false')
    displacement_mode = values.get('displacement_mode', 'full').lower()
    if displacement_mode not in {'full', 'atoms', 'minimal'}:
        raise ValueError('displacement_mode must be full, atoms, or minimal')
    displacement_basis = values.get('displacement_basis', 'cartesian').lower()
    if displacement_basis not in {'cartesian', 'fractional'}:
        raise ValueError('displacement_basis must be cartesian or fractional')
    saved_basis = values.get('saved_basis', 'false').lower()
    if saved_basis not in {'true', 'false'}:
        raise ValueError('saved_basis must be true or false')
    try:
        mode_q = float(values.get('mode_q', '0.01'))
        minimum_frequency = float(values.get('mode_minimum_frequency', '1.0'))
    except ValueError as error:
        raise ValueError('mode_q and mode_minimum_frequency must be numeric') from error
    if not math.isfinite(mode_q) or mode_q <= 0:
        raise ValueError('mode_q must be finite and positive (Angstrom sqrt(amu))')
    if not math.isfinite(minimum_frequency) or minimum_frequency < 0:
        raise ValueError('mode_minimum_frequency must be finite and nonnegative (cm^-1)')
    indices_text = values.get('mode_indices', 'auto')
    mode_indices = None
    if indices_text.lower() != 'auto':
        try:
            mode_indices = tuple(dict.fromkeys(int(i) for i in indices_text.replace(',', ' ').split()))
        except ValueError as error:
            raise ValueError('mode_indices must be auto or positive 1-based integers') from error
        if not mode_indices or min(mode_indices) < 1:
            raise ValueError('mode_indices must be auto or positive 1-based integers')
    return Settings(
        vector('incident_polarization', 0), vector('scattered_polarization', 1), surface_normal,
        energies, displacement_amplitude, hwhm, broadening_type, polarization, frequency_cutoff_cm1,
        difference_scheme, dielectric_source, projection == 'true',
        displacement_mode, displacement_basis, values.get('equilibrium_file'), saved_basis == 'true',
        mode_q, mode_indices, minimum_frequency,
        values.get('mode_displacements_directory', 'mode_displacements'),
        values.get('mode_results_directory', 'mode_results'),
    )


def write_settings(path: str | Path, settings: Settings) -> None:
    if settings.incident_polarization is None or settings.scattered_polarization is None or settings.surface_normal is None:
        raise ValueError("Geometry is required when writing input")
    lines = [
        "incident_polarization: " + " ".join(f"{value:.16g}" for value in settings.incident_polarization),
        "scattered_polarization: " + " ".join(f"{value:.16g}" for value in settings.scattered_polarization),
        f"surface_normal: {settings.surface_normal}",
        "laser_energies: " + " ".join(f"{energy:.2f}" for energy in settings.laser_energies),
    ]
    if settings.displacement_amplitude is not None:
        lines.append("displacement_amplitude: " + " ".join(f"{value:g}" for value in settings.displacement_amplitude))
    if settings.frequency_cutoff_cm1 is not None:
        lines.append(f"cutoff: {settings.frequency_cutoff_cm1:g}")
    lines.extend([
        f"polarization: {settings.polarization}",
        f"broadening_hwhm: {settings.broadening_hwhm:g}",
        f"broadening_type: {settings.broadening_type}",
        f"difference_scheme: {settings.difference_scheme}",
        f"dielectric_source: {settings.dielectric_source}",
        f"site_symmetry_projection: {str(settings.site_symmetry_projection).lower()}",
        f"displacement_mode: {settings.displacement_mode}",
        f"displacement_basis: {settings.displacement_basis}",
        f"saved_basis: {str(settings.saved_basis).lower()}",
        f"mode_q: {settings.mode_q:g}",
        "mode_indices: " + ('auto' if settings.mode_indices is None else ' '.join(map(str, settings.mode_indices))),
        f"mode_minimum_frequency: {settings.mode_minimum_frequency:g}",
        f"mode_displacements_directory: {settings.mode_displacements_directory}",
        f"mode_results_directory: {settings.mode_results_directory}",
    ])
    if settings.equilibrium_file is not None:
        lines.append(f"equilibrium_file: {settings.equilibrium_file}")
    Path(path).write_text("\n".join(lines) + "\n")
