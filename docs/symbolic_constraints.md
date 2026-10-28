# Symbolic atomic derivative constraints

This optional diagnostic solves the exact rank-three Cartesian covariance
equations. It does not modify derivatives, project tensors, change calculations,
or assume that a structurally allowed operation preserves magnetic order.

Install `spectropy[symbolic]`, then run:

```sh
spectropy-symbolic --poscar POSCAR --symmetric-dielectric --output constraints.md
```

Without an installed command, use `PYTHONPATH=SpectroPy/src python -m
spectropy_symbolic` with the same arguments. `--symmetric-dielectric` adds
symmetry under interchange of the two dielectric indices; this is an explicit
optional physical assumption, not a consequence of point-group symmetry.

Alternatively, `--operations group.json` accepts the full electronic-state
unitary group, including atom mappings:

```json
{"natoms": 2, "labels": ["A", "A"], "operations": [
  {"permutation": [0, 1], "rotation": [[1,0,0],[0,1,0],[0,0,1]]},
  {"permutation": [1, 0], "rotation": [[-1,0,0],[0,-1,0],[0,0,-1]]}
]}
```

Permutation indices are zero-based source-to-target mappings. Rotations act on
Cartesian column vectors; fractional rotations must first be converted using
`S = A W A^-1`. Exact symbolic strings such as `"sqrt(3)/2"` are accepted.
Floating entries are recovered using small radicals with a 1e-9 tolerance, then
exact orthogonality and group closure are required. If recovery fails, supply
exact matrices; the solver does not silently solve inconsistent approximate
symmetry equations. For large groups the exact closure check can be expensive.
Antiunitary magnetic operations are not supported.

The solver works with all 27 components jointly, ordered `(i, alpha, beta)`.
It reduces the homogeneous equations from each representative's site subgroup
by symbolic row reduction. The free columns become symbolic parameters;
operations exchanging atoms transfer the representative's form with the same
parameters. Every resulting atom form is independently checked against every
input operation, so no site-to-site constraint is omitted. The output includes
three 3x3 displacement slices per atom, component equations, independent
parameter counts, atom mappings, and exact rotation matrices.

Parameters may be complex, with the same linear constraints on real and
imaginary parts. No acoustic sum rule or mode-specific constraint is imposed.
Allowed components can still contain numerical error. This diagnostic derives
the admissible forms only; subtraction of a fitted symmetry-inconsistent
component is a separate step not implemented here.
