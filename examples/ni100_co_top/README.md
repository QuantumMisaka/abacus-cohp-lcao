# Ni(100)-CO ABACUS LCAO-COHP vs VASP+LOBSTER

This bundle contains the lightweight public artifacts for the spin-polarized
Ni(100)-CO comparison reported in:

- `../../docs/ni100-co-abacus-vs-lobster-cohp.md`
- `../../docs/ni100-co-magnetism-and-cohp-comparability-note.md`

Included files:

- `abacus_relax_efficiency/`: ABACUS LCAO relaxation input using APNS efficiency
  orbitals.
- `abacus_scf_efficiency/`: ABACUS LCAO SCF/COHP input using the same efficiency
  orbital set.
- `abacus_scf_precision/`: ABACUS LCAO SCF/COHP input using APNS precision
  orbitals; this is the primary ABACUS result used in the report.
- `analysis/`: ABACUS efficiency/precision COHP curves, LOBSTER comparison
  curves, overlay plots, and summary data.

The heavy ABACUS `OUT.ABACUS`, VASP `WAVECAR/CHGCAR/vasprun.xml`, and LOBSTER
raw output directories are intentionally not bundled. VASP+LOBSTER results are
included only through processed comparison curves and report tables.

The included ABACUS inputs use relative PP/ORB paths under `../../data/`.
