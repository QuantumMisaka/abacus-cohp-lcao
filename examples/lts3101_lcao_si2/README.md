# LTS 3.10.x Si2 COHP Example

This lightweight bundle stores the input files and generated Si-Si COHP result
from an ABACUS LTS 3.10.x LCAO SCF validation run.

The original SCF output directory is not included, but the result was generated
with `src/cohp.py` from an `OUT.ABACUS` directory containing `data-*-H/S`,
`WFC_NAO_K*.txt`, `kpoints`, and `running_scf.log`. The selected orbital groups
were Si atom 1 global NAO indices `0..12` and Si atom 2 indices `13..25`.

Included files:

- `INPUT`, `KPT`, `STRU`: calculation inputs used for the validation.
- `si_si_COHP.dat`: energy-resolved Si-Si COHP data.
- `si_si_COHP.png`: plotted `-COHP` curve.

Required pseudopotential and orbital files are included in
`../data/legacy-si`, and `INPUT` points to that directory with relative paths.

The heavy `OUT.ABACUS` directory is intentionally not included.
