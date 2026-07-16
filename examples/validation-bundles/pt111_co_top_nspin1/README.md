# Pt(111)-CO Top-Site nspin=1 COHP Example

This lightweight bundle retains the inputs and orbital mapping for the
historical nspin=1 Pt(111)-CO top-site run. Its pre-2026-07-17 processed curves
and summaries were withdrawn after correction of the COHP conventions.

COHP was computed with `src/cohp.py` from the final `OUT.ABACUS` matrices,
wavefunctions, k-point weights, and Fermi level. `mapping.json` records the
ABACUS global NAO ranges used to select top Pt and C orbitals for the total
Pt-C, Pt-d/C-p, and Pt-d/C-s curves.

The final SCF input uses the full required COHP output block:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

Included files:

- `INPUT`, `STRU`, `run_abacus.sbatch`: compact-basis relax input.
- `INPUT.final_scf`, `STRU.final_scf`: precision-basis final SCF input.
- `mapping.json`: atom and orbital-range mapping used for COHP selection.

Required APNS files are included in `../data/PP`, `../data/ORB`, and
`../data/apns-orbitals-precision-v1`; both `INPUT` files point to those
directories with relative paths.

The heavy `OUT.ABACUS` directory is intentionally not included.
