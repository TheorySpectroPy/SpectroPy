"""Locate the raw displacement container; accept legacy standalone datasets."""
from pathlib import Path


def displacement_root(directory='.'):
    root = Path(directory)
    if root.name == 'displacements':
        return root
    child = root / 'displacements'
    return child if child.is_dir() else root
