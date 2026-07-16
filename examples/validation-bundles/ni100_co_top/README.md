# Ni(100)-CO ABACUS LCAO-COHP vs VASP+LOBSTER

This bundle retains reproduction inputs for the historical spin-polarized
Ni(100)-CO comparison. Its pre-2026-07-17 processed curves were withdrawn after
correction of the COHP normalization and integration conventions. The archived
reports are context only:

- `../../docs/archive/ni100-co-abacus-vs-lobster-cohp.md`
- `../../docs/archive/ni100-co-magnetism-and-cohp-comparability-note.md`

Included files:

- `abacus_relax_efficiency/`: ABACUS LCAO relaxation input using APNS efficiency
  orbitals.
- `abacus_scf_efficiency/`: ABACUS LCAO SCF/COHP input using the same efficiency
  orbital set.
- `abacus_scf_precision/`: ABACUS LCAO SCF/COHP input using APNS precision
  orbitals; this is the primary ABACUS result used in the report.

The heavy ABACUS `OUT.ABACUS`, VASP `WAVECAR/CHGCAR/vasprun.xml`, and LOBSTER
raw output directories are intentionally not bundled.

The included ABACUS inputs use relative PP/ORB paths under `../../data/`.
