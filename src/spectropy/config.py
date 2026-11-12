from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from fractions import Fraction
from pathlib import Path
import math
import re

import yaml


DEFAULT_INPUT = "input.yaml"


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
    site_symmetry_projection: bool = False
    symmetry: str = "off"
    displacement_basis: str = "cartesian"
    equilibrium_file: str | None = None
    saved_basis: bool = False
    mode_q: float = 0.01
    mode_indices: tuple[int, ...] | None = None
    mode_minimum_frequency: float = 1.0
    mode_displacements_directory: str = "mode_displacements"
    mode_results_directory: str = "mode_results"
    per_atom_output: bool = False
    raw_directional_derivatives: bool = False

    @property
    def broadening_fwhm(self):
        """Full width for internal line-shape routines; public settings use HWHM."""
        return 2 * self.broadening_hwhm


# Public YAML paths map to the existing internal settings used by the solvers.
_SECTIONS = {
    "polarization": {
        "incident": "incident_polarization", "scattered": "scattered_polarization",
        "surface_normal": "surface_normal", "mode": "polarization",
    },
    "spectrum": {"energies": "laser_energies", "cutoff": "frequency_cutoff_cm1"},
    "displacements": {
        "basis": "displacement_basis", "amplitude": "displacement_amplitude",
        "difference": "difference_scheme", "equilibrium_file": "equilibrium_file",
    },
    "symmetry": {
        "reduction": "symmetry", "denoise": "site_symmetry_projection",
        "saved_basis": "saved_basis",
    },
    "broadening": {"hwhm": "broadening_hwhm", "type": "broadening_type"},
    "modes": {
        "q": "mode_q", "indices": "mode_indices",
        "minimum_frequency": "mode_minimum_frequency",
        "displacements_directory": "mode_displacements_directory",
        "results_directory": "mode_results_directory",
    },
    "output": {"per_atom": {
        "enabled": "per_atom_output",
        "raw_directional_derivatives": "raw_directional_derivatives",
    }},
}


def _setting_paths(schema, prefix=""):
    for key, value in schema.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            yield from _setting_paths(value, path)
        else:
            yield value, path


_SETTING_PATHS = dict(_setting_paths(_SECTIONS))


def _flatten_sections(sections, schema=_SECTIONS, prefix=""):
    if not isinstance(sections, dict):
        raise ValueError(f"{prefix or 'Settings'} must be a YAML mapping")
    values = {}
    for key, value in sections.items():
        path = f"{prefix}.{key}" if prefix else key
        if key not in schema:
            raise ValueError(f"Unknown input setting: {path}")
        if isinstance(schema[key], dict):
            values.update(_flatten_sections(value, schema[key], path))
        else:
            values[schema[key]] = value
    return values


class _SettingsLoader(yaml.SafeLoader):
    """Safe YAML with duplicate and non-string setting names rejected."""

    # Resolve booleans as true/false only, so reduction: on/off remain strings.
    # Copy the resolver table: other YAML readers must retain their behavior.
    yaml_implicit_resolvers = {
        key: [(tag, pattern) for tag, pattern in resolvers
              if tag != 'tag:yaml.org,2002:bool']
        for key, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
    }

    def construct_mapping(self, node, deep=False):
        result = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, str):
                raise ValueError("Setting names must be strings")
            if key in result:
                raise ValueError(f"Duplicate input setting: {key}")
            result[key] = self.construct_object(value_node, deep=deep)
        return result


_SettingsLoader.add_implicit_resolver(
    'tag:yaml.org,2002:bool', re.compile(r'^(?:true|True|TRUE|false|False|FALSE)$'), list('tTfF')
)


def _number(value, name):
    name = _SETTING_PATHS.get(name, name)
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"{name} must be numeric")
    try:
        number = float(value)
    except (ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be numeric") from error
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _laser_energies(value):
    """Expand an inclusive decimal range without cumulative float-step error."""
    name = "spectrum.energies"
    if isinstance(value, dict):
        if set(value) != {"start", "stop", "step"}:
            raise ValueError(f"{name} range requires exactly start, stop, and step")
        start, stop, step = (Fraction(str(_number(value[key], f"{name}.{key}")))
                             for key in ("start", "stop", "step"))
        if start < 0 or stop < start or step <= 0:
            raise ValueError(f"{name} range requires 0 <= start <= stop and a positive step")
        intervals = (stop - start) / step
        if intervals.denominator != 1:
            raise ValueError(f"{name}.stop must lie on the start + n*step grid (inclusive)")
        count = intervals.numerator + 1
        if count > 100_000:
            raise ValueError(f"{name} range is too large (maximum 100000 values)")
        energies = tuple(float(start + i * step) for i in range(count))
    else:
        items = value if isinstance(value, list) else [value]
        if not items:
            raise ValueError(f"{name} must contain at least one value")
        energies = tuple(dict.fromkeys(_number(item, name) for item in items))
    if min(energies) < 0:
        raise ValueError(f"{name} must be non-negative")
    if len({f"{energy:.2f}" for energy in energies}) != len(energies):
        raise ValueError(f"{name} values must have distinct filenames at two decimal places")
    return energies


def _settings_from_mapping(values):
    if not isinstance(values, dict):
        raise ValueError("Settings must be a YAML mapping of named keys and values")
    values = values.copy()
    allowed = {field.name for field in fields(Settings)}
    unknown = set(values) - allowed
    if unknown:
        raise ValueError(f"Unknown input setting(s): {', '.join(sorted(map(str, unknown)))}")
    defaults = Settings()

    def vector(name):
        value = values.get(name)
        if value is None:
            return None
        if not isinstance(value, list) or len(value) != 3:
            raise ValueError(f"{_SETTING_PATHS[name]} must be a YAML list of three finite numbers, e.g. [1, 0, 0]")
        return tuple(_number(component, name) for component in value)

    def numbers(name, default):
        value = values.get(name, default)
        items = value if isinstance(value, list) else [value]
        if not items:
            raise ValueError(f"{_SETTING_PATHS[name]} must contain at least one value")
        return tuple(_number(item, name) for item in items)

    def choice(name, options, aliases=None):
        value = values.get(name, getattr(defaults, name))
        if isinstance(value, str):
            value = value.lower()
            value = (aliases or {}).get(value, value)
            if value in options:
                return value
        raise ValueError(f"{_SETTING_PATHS[name]} must be {' or '.join(options)}")

    def boolean(name):
        value = values.get(name, getattr(defaults, name))
        if not isinstance(value, bool):
            raise ValueError(f"{_SETTING_PATHS[name]} must be a YAML boolean (true or false)")
        return value

    def path_value(name):
        value = values.get(name, getattr(defaults, name))
        if value is None and name == "equilibrium_file":
            return None
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{_SETTING_PATHS[name]} must be a nonempty path string")
        return value

    def bounded_number(name, default, *, positive=False):
        value = _number(values.get(name, default), name)
        if value < 0 or (positive and value == 0):
            raise ValueError(f"{_SETTING_PATHS[name]} must be {'positive' if positive else 'non-negative'}")
        return value

    surface_normal = (choice("surface_normal", ("x", "y", "z"))
                      if values.get("surface_normal") is not None else None)
    energies = _laser_energies(values.get("laser_energies", list(defaults.laser_energies)))
    amplitude = None
    if values.get("displacement_amplitude") is not None:
        amplitude = numbers("displacement_amplitude", None)
        if len(amplitude) == 1:
            amplitude *= 3
        if len(amplitude) != 3:
            raise ValueError("displacements.amplitude takes one value or a YAML list of three basis-axis lengths")
        if min(amplitude) <= 0:
            raise ValueError("displacements.amplitude must be positive")

    indices = values.get("mode_indices", "auto")
    if indices is None or indices == "auto":
        indices = None
    elif (isinstance(indices, list) and indices and
          all(type(index) is int and index > 0 for index in indices)):
        indices = tuple(dict.fromkeys(indices))
    else:
        raise ValueError("modes.indices must be auto or a YAML list of positive 1-based integers")

    per_atom_output = boolean("per_atom_output")
    raw_directional_derivatives = boolean("raw_directional_derivatives")
    if raw_directional_derivatives and not per_atom_output:
        raise ValueError("output.per_atom.raw_directional_derivatives requires output.per_atom.enabled: true")

    return Settings(
        incident_polarization=vector("incident_polarization"),
        scattered_polarization=vector("scattered_polarization"),
        surface_normal=surface_normal,
        laser_energies=energies,
        displacement_amplitude=amplitude,
        broadening_hwhm=bounded_number("broadening_hwhm", defaults.broadening_hwhm, positive=True),
        broadening_type=choice("broadening_type", ("lorentzian", "gaussian"),
                               {"l": "lorentzian", "g": "gaussian"}),
        polarization=choice("polarization", ("specific", "average")),
        frequency_cutoff_cm1=(bounded_number("frequency_cutoff_cm1", None)
                              if values.get("frequency_cutoff_cm1") is not None else None),
        difference_scheme=choice("difference_scheme", ("central", "forward")),
        site_symmetry_projection=boolean("site_symmetry_projection"),
        symmetry=choice("symmetry", ("on", "off", "atoms")),
        displacement_basis=choice("displacement_basis", ("cartesian", "fractional")),
        equilibrium_file=path_value("equilibrium_file"),
        saved_basis=boolean("saved_basis"),
        mode_q=bounded_number("mode_q", defaults.mode_q, positive=True),
        mode_indices=indices,
        mode_minimum_frequency=bounded_number("mode_minimum_frequency", defaults.mode_minimum_frequency),
        mode_displacements_directory=path_value("mode_displacements_directory"),
        mode_results_directory=path_value("mode_results_directory"),
        per_atom_output=per_atom_output,
        raw_directional_derivatives=raw_directional_derivatives,
    )


def read_settings(path: str | Path = DEFAULT_INPUT) -> Settings:
    path = Path(path)
    if not path.exists():
        if path.name == DEFAULT_INPUT and path.with_name("input").is_file():
            raise ValueError(
                f"Found legacy {path.with_name('input')} but no {path}. "
                "Convert it to grouped YAML settings; see docs.md."
            )
        return Settings()
    try:
        values = yaml.load(path.read_text(), Loader=_SettingsLoader)
        return _settings_from_mapping(_flatten_sections({} if values is None else values))
    except (yaml.YAMLError, ValueError) as error:
        raise ValueError(f"Invalid settings in {path}: {error}") from error


def write_settings(path: str | Path, settings: Settings) -> None:
    """Write grouped YAML without rounding scientific settings."""
    def mapping(settings):
        return {
            name: list(value) if isinstance(value, tuple) else value
            for name, value in asdict(settings).items() if value is not None
        }

    values = mapping(_settings_from_mapping(mapping(settings)))
    def nested(schema):
        return {key: nested(value) if isinstance(value, dict) else values[value]
                for key, value in schema.items() if isinstance(value, dict) or value in values}

    Path(path).write_text(yaml.safe_dump(nested(_SECTIONS), sort_keys=False))
