# Pt(111)-CO Top-Site VASP+LOBSTER Comparison Artifacts

This directory contains lightweight artifacts for the spin-polarized
Pt(111)-CO top-site ABACUS LCAO-COHP vs VASP+LOBSTER comparison reported in
`../../docs/pt111-co-abacus-vs-lobster-cohp.md`.

Included files:

- `analysis/summary.json`: quantitative Pt-C and C-O comparison, LOBSTER
  spilling, and output paths.
- `analysis/pt111_cohp_compare.csv`: aligned ABACUS and LOBSTER curves plus
  scaled LOBSTER curves.
- `analysis/*overlay.png`: direct and normalized overlays for Pt-C and C-O.
- `analysis/abacus_*`: regenerated ABACUS comparison curves used for the common
  energy-window comparison.
- `jobs.json`: Slurm job ids for the VASP SCF and LOBSTER runs.

Heavy VASP and LOBSTER runtime outputs are intentionally not bundled.

