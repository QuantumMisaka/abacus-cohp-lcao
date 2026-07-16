# LTS 3.10.x Si2 COHP Example

This lightweight bundle stores the input files for an ABACUS LTS 3.10.x LCAO
SCF validation run. The pre-2026-07-17 generated Si-Si COHP curve was withdrawn
after correction of the COHP unit, pair-counting, and ICOHP conventions.

The original SCF output directory is not included, but the result was generated
with `src/cohp.py` from an `OUT.ABACUS` directory containing `data-*-H/S`,
`WFC_NAO_K*.txt`, `kpoints`, and `running_scf.log`. The selected orbital groups
were Si atom 1 global NAO indices `0..12` and Si atom 2 indices `13..25`.

The SCF `INPUT` includes the required COHP output block:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

Included files:

- `INPUT`, `KPT`, `STRU`: calculation inputs used for the validation.

Regenerate the ABACUS curve and VASP+LOBSTER comparison before publishing a
quantitative result from this example.

Required pseudopotential and orbital files are included in
`../data/legacy-si`, and `INPUT` points to that directory with relative paths.

The heavy `OUT.ABACUS` directory is intentionally not included.
