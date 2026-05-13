# Pt(111)-CO Top-Site nspin=1 COHP Example

This lightweight bundle stores the nspin=1 Pt(111)-CO top-site result.

The final SCF input enables `out_mat_hs 1 8` and `out_wfc_lcao 1`. COHP was then
computed with `src/cohp.py` from the final `OUT.ABACUS` matrices,
wavefunctions, k-point weights, and Fermi level. `mapping.json` records the
ABACUS global NAO ranges used to select top Pt and C orbitals for the total
Pt-C, Pt-d/C-p, and Pt-d/C-s curves.

Included files:

- `INPUT`, `STRU`, `run_abacus.sbatch`: compact-basis relax input.
- `INPUT.final_scf`, `STRU.final_scf`: precision-basis final SCF input.
- `mapping.json`: atom and orbital-range mapping used for COHP selection.
- `summary.md`, `summary.json`: COHP result summary.
- `top_Pt_C_total_sum.*`: total top Pt-C COHP result.
- `top_Pt_d_C_p_sum.*`: Pt-d / C-p component.
- `top_Pt_d_C_s_sum.*`: Pt-d / C-s component.

The heavy `OUT.ABACUS` directory is intentionally not included.
