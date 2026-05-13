# Example Pseudopotentials and Orbitals

This directory contains the minimal pseudopotential and numerical-orbital files
needed by the bundled examples.

## APNS Source

The Pt(111)-CO example uses files from the ABACUS APNS PP/ORB library:

https://store.aissquare.com/datasets/dc875646-a526-41f1-a180-d54b218fc80a/ABACUS-APNS-PPORBs-v1.zip

Included APNS subsets:

- `PP/`: PBE pseudopotentials for Pt, C, and O.
- `ORB/`: compact APNS orbitals used by the relaxation input.
- `apns-orbitals-precision-v1/`: precision APNS orbitals used by final SCF and
  COHP post-processing.

## Legacy Si Validation Files

The Si2 validation example reproduces the original ABACUS LTS 3.10.x validation
run and therefore keeps the exact Si files used in that run:

- `legacy-si/Si.pz-vbc.UPF`
- `legacy-si/Si_lda_8.0au_50Ry_2s2p1d`

These files came from the ABACUS test pseudopotential/orbital set rather than
the Pt-CO APNS PBE subset.
