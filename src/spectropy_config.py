from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


_SETTING_RE = re.compile(r"^([a-z_]+)\s*(?::|=|\s)\s*(.+)$", re.I)


@dataclass(frozen=True)
class Settings:
    incident_polarization: tuple[float, float, float] | None = None
    scattered_polarization: tuple[float, float, float] | None = None
    surface_normal: str | None = None
    laser_energies: tuple[float, ...] = (0.0,)
    displacement_amplitude: float | None = None
    broadening_fwhm: float = 1.0
    broadening_type: str = "lorentzian"
    polarization: str = "specific"


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
            "broadening_fwhm", "broadening_type", "polarization",
        }:
            values[match.group(1).lower()] = match.group(2)
        else:
            geometry.append(line.split())

    def vector(index: int) -> tuple[float, float, float] | None:
        if len(geometry) <= index:
            return None
        try:
            return tuple(float(value) for value in geometry[index][:3])  # type: ignore[return-value]
        except ValueError as error:
            raise ValueError("The first two input lines must be three-component polarization vectors") from error

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
            displacement_amplitude = float(amplitude.split()[0])
        except ValueError as error:
            raise ValueError(f"Invalid displacement_amplitude setting: {amplitude}") from error
        if displacement_amplitude <= 0:
            raise ValueError("displacement_amplitude must be positive")
    else:
        displacement_amplitude = None

    try:
        fwhm = float(values.get("broadening_fwhm", "1.0").split()[0])
    except ValueError as error:
        raise ValueError("broadening_fwhm must be numeric") from error
    if fwhm <= 0:
        raise ValueError("broadening_fwhm must be positive")
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

    return Settings(
        vector(0), vector(1), geometry[2][0].lower() if len(geometry) >= 3 else None,
        energies, displacement_amplitude, fwhm, broadening_type, polarization,
    )


def write_settings(path: str | Path, settings: Settings) -> None:
    if settings.incident_polarization is None or settings.scattered_polarization is None or settings.surface_normal is None:
        raise ValueError("Geometry is required when writing input")
    lines = [
        " ".join(f"{value:.1f}" for value in settings.incident_polarization),
        " ".join(f"{value:.1f}" for value in settings.scattered_polarization),
        settings.surface_normal,
        "laser_energies: " + " ".join(f"{energy:.2f}" for energy in settings.laser_energies),
    ]
    if settings.displacement_amplitude is not None:
        lines.append(f"displacement_amplitude: {settings.displacement_amplitude:g}")
    lines.extend([
        f"polarization: {settings.polarization}",
        f"broadening_fwhm: {settings.broadening_fwhm:g}",
        f"broadening_type: {settings.broadening_type}",
    ])
    Path(path).write_text("\n".join(lines) + "\n")
