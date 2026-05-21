# Si2 VASP+LOBSTER Comparison Artifacts

This directory contains lightweight artifacts for the Si2 ABACUS LCAO-COHP vs
VASP+LOBSTER comparison reported in
`../../../docs/si2-abacus-vs-vasp-lobster-cohp.md`.

Included files:

- `summary.json`: quantitative comparison, LOBSTER spilling, and output paths.
- `si2_cohp_compare.csv`: aligned ABACUS and LOBSTER curves plus scaled LOBSTER
  curve.
- `si2_si_si_minus_cohp_overlay.png`: direct `-COHP/-pCOHP` overlay.
- `si2_si_si_minus_cohp_normalized_overlay.png`: best-fit normalized overlay.
- `abacus_si_si.dat` and `abacus_si_si.png`: regenerated ABACUS comparison curve.
- `jobs.json`: Slurm job ids for the VASP SCF and LOBSTER runs.

Heavy VASP and LOBSTER runtime outputs are intentionally not bundled.

