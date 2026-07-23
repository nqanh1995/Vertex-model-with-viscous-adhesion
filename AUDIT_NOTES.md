# Audit Notes

This document summarizes the code audit and cleanup performed on this
repository. **No physics/math or algorithmic behavior was changed** except
where explicitly called out below — every file with dense analytical
formulas (the Hessian/shear-modulus calculations) was diffed line-by-line
against the original to confirm only comments and whitespace changed.

All 28 `.m` files were reviewed. Every function now has a header comment
explaining its purpose, inputs, and outputs, plus inline comments on
non-obvious steps. Simple/self-explanatory lines were left uncommented
per the request. Dead/commented-out code blocks were removed throughout
(e.g. duplicate rotation logic in the T1 functions, an alternate
`convhull`-based area calculation, several duplicated derivative blocks).

## Bugs found

1. **`harmonic_oscillatory_shear.m` — undefined variable `T1_timer`
   (now resolved).** The function referenced `T1_timer` in two places
   (the emergency-stop branch, and the normal per-interval save branch)
   without ever defining it. Per clarification: this function is
   intentionally used to probe the tissue's **linear response regime**
   (small `strain_amp`), and T1 (neighbor-exchange) transitions are
   deliberately **not** performed during the main oscillatory-shear
   loop, since a discrete topological rearrangement would introduce
   nonlinear/plastic relaxation and invalidate the linear-response
   measurement. Given that, the `T1_timer` references were vestigial
   copy-paste leftovers from `avm_viscous.m` (which *does* track T1
   timers) rather than a real feature of this function, so they have
   been **removed**, and the function's header docstring now states the
   linear-response purpose and the T1-free design explicitly.

2. **`avm_viscous.m` — `toc` with no matching `tic` (fixed).** The
   original computed `elapsed_time = toc` at the end without ever
   calling `tic`, which reports a stale or incorrect elapsed time (or
   errors, depending on session state). Added `tic;` right before the
   main loop starts.

3. **`draw_state_LE.m` — missing-argument crash (fixed).** The function
   signature requires `c`, `lw`, `opac` as positional arguments with no
   defaults, but every call site (in both `avm_viscous.m` and
   `harmonic_oscillatory_shear.m`) passes only 4 arguments. This means
   the plotting calls in the emergency-stop diagnostic branches of both
   driver scripts would have crashed. Added default values
   (`c = 'k'`, `lw = 1`, `opac = 1`) so the function works both with and
   without these arguments.

4. **`harmonic_oscillatory_shear.m` — inconsistent T1 iteration
   counter (not fixed, noted here only).** During the initial relaxation
   loop, `n_t1` is incremented once per edge flipped *and* once per outer
   while-loop iteration, whereas the equivalent loop in `avm_viscous.m`
   only increments once per outer iteration. This isn't a crash, but it
   means the 500-iteration safety cap on the T1 clean-up loop behaves
   differently between the two scripts. Left as-is since it doesn't
   affect correctness, only how quickly the cap is reached.

5. **`vm_minimize.m` — missing dependency `get_VM_force.m` (now
   resolved).** This function previously called a non-sheared
   `get_VM_force`, which is not included anywhere in the repository. Per
   clarification: `get_VM_force_LE` evaluated at `strain = 0` is exactly
   equivalent to the flat (non-sheared) force calculation, so the call
   was replaced with `get_VM_force_LE(V, C, KP, P0, KA, A0, 0)`. This
   removes the missing dependency entirely — `vm_minimize.m` now only
   needs `lees_edwards_functions/` on the path (in addition to its own
   folder), which is noted in its header comment.

## Removed files

- **`lees_edwards_functions/unshear_initial_state_qvm.m`** — removed.
  This function was never called anywhere in the codebase (not from
  `run_sim.m`, `avm_viscous.m`, or `harmonic_oscillatory_shear.m`), so it
  wasn't needed to run the simulation. It also depended on a missing
  `get_qvm_stress.m` function that doesn't exist in this repository (a
  similarly named `get_viscous_qvm_stress.m` exists but has a different
  signature — it requires viscous force inputs `Fvx, Fvy` that this
  function never computed, so it wasn't a safe drop-in replacement).
  Since the function was both unused and unrunnable as written, it was
  deleted rather than patched. If you want this "find the strain that
  zeroes out shear stress" utility back in the future, it would need to
  be rewritten against `get_viscous_qvm_stress` (deciding what to pass
  for `Fvx, Fvy`, likely zeros if no viscous coupling forces are
  involved in that calculation) or against a newly written
  non-viscous stress function.

## Missing dependencies (not included in this repository)

These functions are called but not present anywhere in the codebase, so
they must exist on your local MATLAB path for these files to run:

- **`polygeom.m`** — a third-party MATLAB File Exchange utility used by
  `matlab_functions/cell_shape_info.m` for polygon geometry/principal
  moments.

## Minor cleanups (no behavior change)

- Replaced `(coeff_mat^-1) * A` with `coeff_mat \ A` in
  `get_vertex_velo_viscous.m` — mathematically identical, but avoids
  forming an explicit matrix inverse (faster and more numerically
  stable; standard MATLAB guidance).
- Normalized mixed tab/space indentation to spaces in the three files
  with dense analytical-derivative blocks (`get_dE2_mat_vm_LE.m`,
  `get_qvm_hessian.m`, `get_shear_modulus_LE.m`) without touching any
  formula content (verified by diff).
- Removed a handful of dead commented-out code blocks throughout (see
  individual files) that were either superseded alternatives or
  duplicated logic.
- Factored the repeated x/y unwrapping logic in
  `matlab_functions/get_unwrap_xy.m` into a small local helper function
  to avoid duplicating the same ~10 lines twice; output is identical.

## Files reviewed

**Top level:** `run_sim.m`, `avm_viscous.m`, `harmonic_oscillatory_shear.m`

**`matlab_functions/`:** `pbc_dist.m`, `vertices_to_cw_polygon.m`,
`drifting_correction.m`, `get_cell_shape.m`, `cell_shape_info.m`,
`get_unwrap_xy.m`, `get_theta_vertex.m`, `get_vertex_velo_viscous.m`,
`get_complex_modulus.m`, `get_adj_mat.m`, `draw_state.m`,
`T1_edgeswap.m`, `vm_minimize.m`, `make_voronoi_pbc_fast.m`

**`lees_edwards_functions/`:** `pbc_dist_LE.m`, `position_mod_LE.m`,
`get_VM_force_LE.m`, `get_cell_shape_LE.m`, `vm_minimize_LE.m`,
`T1_edgeswap_LE.m`, `draw_edges_LE.m`,
`draw_state_LE.m`, `get_dE2_mat_vm_LE.m`, `get_qvm_hessian.m`,
`calc_stress.m`, `get_shear_modulus_LE.m`, `get_viscous_qvm_stress.m`

(`unshear_initial_state_qvm.m` was reviewed but subsequently removed —
see "Removed files" above.)
