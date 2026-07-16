# Diamond ABACUS LCAO-COHP vs VASP+LOBSTER

This bundle retains reproduction inputs for the historical diamond comparison.
Its pre-2026-07-17 processed curves were withdrawn after correction of the COHP
normalization and integration conventions. The archived report is context only.

Included files:

- `abacus_relax/`: ABACUS LCAO cell-relax input using APNS efficiency orbitals.
- `abacus_scf/`: ABACUS LCAO precision-basis SCF input used for COHP
  post-processing.
- `POSCAR.initial`, `STRU.initial`, `mapping_initial.json`: initial structure
  and orbital mapping metadata.

The heavy ABACUS `OUT.ABACUS`, VASP `WAVECAR/CHGCAR/vasprun.xml`, and LOBSTER
raw output directories are intentionally not bundled.

To reproduce the ABACUS-side SCF output, run the inputs from their directories
after making sure `../../data/PP`, `../../data/ORB`, and
`../../data/apns-orbitals-precision-v1` are present.
