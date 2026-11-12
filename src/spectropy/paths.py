"""Locate the raw displacement container; accept legacy standalone datasets."""
from pathlib import Path
import os
import warnings


def displacement_root(directory='.'):
    root = Path(directory)
    if root.name == 'displacements':
        return root
    child = root / 'displacements'
    return child if child.is_dir() else root


def link_vasp_inputs(directories, source_directory='.'):
    """Link supplied optical inputs without replacing any existing inputs.

    Check every destination before creating links. Missing source inputs are
    reported, but do not prevent structure-only displacement generation.
    """
    source_directory = Path(source_directory).absolute()
    directories = list(dict.fromkeys(Path(d).absolute() for d in directories))
    sources = []
    missing = []
    for name in ('INCAR', 'KPOINTS', 'POTCAR'):
        source = source_directory / name
        if source.is_file():
            sources.append(source)
        elif source.exists():
            raise ValueError(f'Expected a file for VASP input: {source}')
        else:
            missing.append(name)
    planned = []
    for directory in directories:
        if directory.exists() and not directory.is_dir():
            raise FileExistsError(f'Calculation directory is not a directory: {directory}')
        for source in sources:
            target = directory / source.name
            if target.is_symlink() or target.exists():
                if target.is_symlink() and target.resolve() == source.resolve():
                    continue
                raise FileExistsError(f'Refusing to replace existing VASP input: {target}')
            planned.append((source, target))
    for source, target in planned:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.symlink_to(os.path.relpath(source, target.parent))
    if missing:
        warnings.warn(f'No supplied {", ".join(missing)} in {source_directory}; '
                      'links for these inputs were not created. Supply them before running VASP.',
                      RuntimeWarning, stacklevel=2)
    return [target for _, target in planned]
