# ABACUS LCAO COHP

This repository contains a lightweight ABACUS LCAO-COHP post-processing workflow.
It extracts Hamiltonian, overlap, wavefunction, eigenvalue, Fermi-level, and k-point
weight information from ABACUS LCAO SCF outputs, then evaluates atom-pair or
orbital-group COHP curves.

This directory is the GitHub publication package for the local ABACUS-COHP
working tree. In the original workspace, it is also exposed as the shortcut
`abacus-cohp/abacus-cohp-lcao`, pointing to this standalone repository for
release and sharing.

The implementation is intended for ABACUS numerical atomic orbital analysis. Its
COHP values are ABACUS-NAO dependent and should not be assumed numerically
equivalent to LOBSTER pCOHP.

## Project Origin

This repository is based on the ABACUS LCAO-COHP development discussed in
deepmodeling/abacus-develop issue #3718:

https://github.com/deepmodeling/abacus-develop/issues/3718

The ABACUS-LCAO-COHP implementation and development work should be credited to
the GitHub user `kirk0830`. This repository does not claim original authorship
of that method implementation; it summarizes usage, packages the key scripts,
and provides reproducible ABACUS LCAO post-processing examples built on that
existing work.

## Contents

- `src/read_abacus_out.py`: ABACUS output readers for `data-*-H/S`, `WFC_NAO_K*.txt`,
  `kpoints`, `running_scf.log`, and helper utilities.
- `src/cohp.py`: COHP/COOP and pCOHP-like post-processing routines plus a CLI.
- `scripts/pt111_co_workflow.py`: reproducible Pt(111)-CO example workflow for SAI
  Slurm + ABACUS LTS v3.10.1.
- `scripts/diamond_abacus_pw_lobster_probe.py`: diagnostic workflow used to
  test whether ABACUS PW outputs can be consumed by LOBSTER or converted to the
  LOBSTER Generic interface.
- `docs/quickstart-abacus-scf-to-cohp.md`: fast user guide for the ABACUS
  SCF -> `src/cohp.py` post-processing workflow.
- `docs/`: method notes, validation notes, quick start guide, and Pt(111)-CO
  result report.
- `examples/`: lightweight result bundles for the two completed examples.
- `examples/data/`: pseudopotentials and numerical orbitals required by the
  bundled examples.

The `src/` directory is the core COHP calculation layer. In normal use,
`read_abacus_out.py` should not be called manually except for reader tests;
`cohp.py` is the user-facing post-processing command.

## Required ABACUS Output

This workflow starts from a completed ABACUS LCAO SCF calculation, not from a
plane-wave calculation. The SCF input must request the matrices and NAO
wavefunctions needed by the post-processor:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

All four lines above should be treated as the minimal COHP-output block for this
repository. The earlier documentation only emphasized `out_mat_hs 1 8` and
`out_wfc_lcao 1`; the completed Si2, Pt(111)-CO, diamond, and Ni(100)-CO tests
all used `out_app_flag 1` as well. ABACUS documents `out_app_flag` as the switch
controlling append-style output for LCAO `H(k)`, `S(k)`, and `wfc(k)` matrix
families together with `out_mat_hs` and `out_wfc_lcao`. Keeping it explicit makes
the generated `OUT.ABACUS` layout match the reader expectations used here.

Relevant ABACUS documentation:

- `out_mat_hs`: https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/hs_matrix.html
- `out_wfc_lcao` and `out_app_flag`: https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html

After SCF, the target `OUT.ABACUS` directory must contain:

- `data-*-H` and `data-*-S`: Hamiltonian and overlap matrices for each k point.
- `WFC_NAO_K*.txt` or `WFC_NAO_GAMMA*.txt`: LCAO eigenvectors in the ABACUS NAO basis.
- `kpoints`: k-point weights.
- `running_scf.log`: Fermi level and band information.

The user must provide orbital-index lists for the two atoms or orbital groups to
be analyzed. These indices are global NAO indices in the ABACUS basis order. For
example, if atom A owns orbitals `0..12` and atom B owns orbitals `13..25`, the
pair COHP is computed by passing those two lists to `--atom-i-orbs` and
`--atom-j-orbs`.

For a complete first-run workflow, start from
`docs/quickstart-abacus-scf-to-cohp.md`.

## Pseudopotentials and Orbitals

The Pt(111)-CO example uses the ABACUS APNS PP/ORB library. The upstream dataset
used for this repository is:

```text
https://store.aissquare.com/datasets/dc875646-a526-41f1-a180-d54b218fc80a/ABACUS-APNS-PPORBs-v1.zip
```

Only the files needed by the examples are included under `examples/data/`:

- `examples/data/PP`: Pt, C, and O pseudopotentials.
- `examples/data/ORB`: compact Pt, C, and O orbitals for relaxation.
- `examples/data/apns-orbitals-precision-v1`: precision Pt, C, and O orbitals
  for final SCF and COHP output.
- `examples/data/legacy-si`: the Si pseudopotential and orbital used by the Si2
  validation example.

Example `INPUT` files use relative paths such as `../data/PP` and
`../data/apns-orbitals-precision-v1`, so they can be run from their own example
directories without relying on a user-specific `~/PP_ORB` location.

## Minimal COHP Usage

Run on an ABACUS `OUT.ABACUS` directory containing:

- `data-*-H`
- `data-*-S`
- `WFC_NAO_K*.txt` or compatible text wavefunction files
- `kpoints`
- `running_scf.log`

The corresponding SCF `INPUT` must include:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

Example:

```bash
env MPLBACKEND=Agg python src/cohp.py \
  --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2,3 \
  --atom-j-orbs 13,14,15,16 \
  --method COHP \
  --de 0.05 \
  --smooth-nstddev 4 \
  --invert \
  --output-prefix my_pair_COHP
```

This command writes:

- `my_pair_COHP.dat`: two-column energy and COHP data.
- `my_pair_COHP.png`: plotted curve. With `--invert`, the figure follows the
  common `-COHP` plotting convention where occupied bonding contributions appear
  positive.

For spin-polarized ABACUS outputs:

```bash
python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2 \
  --atom-j-orbs 100,101,102 \
  --spin up --output-prefix pair_up

python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2 \
  --atom-j-orbs 100,101,102 \
  --spin down --output-prefix pair_down

python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2 \
  --atom-j-orbs 100,101,102 \
  --spin sum --output-prefix pair_sum
```

## Pt(111)-CO Workflow

The SAI-oriented workflow builds Pt(111) 2x2x4 + top-site CO with vacuum along
ABACUS B direction, writes ABACUS relax/final-SCF inputs, and summarizes COHP.

```bash
/home/pku-jianghong/liuzhaoqing/.conda/envs/ase/bin/python \
  scripts/pt111_co_workflow.py prepare
```

The generated Slurm scripts use:

- `module load abacus/LTSv3.10.1-sm70-auto`
- 4V100 partition
- 1 GPU
- no `mpirun`
- `OMP_NUM_THREADS=8`
- `ks_solver cusolver`

After relax finishes:

```bash
python scripts/pt111_co_workflow.py prepare-final --nspin 1
python scripts/pt111_co_workflow.py analyze --nspin 1
```

Use `--nspin 2` for the spin-polarized branch.

## Included Examples

### ABACUS LTS 3.10.x Si2

Located in `examples/lts3101_lcao_si2`.

This example confirms that ABACUS LTS 3.10.x LCAO SCF output can be processed
from `data-*-H/S` and `WFC_NAO_K*.txt` into a Si-Si COHP curve.

### Pt(111)-CO Top Site, nspin=1

Located in `examples/pt111_co_top_nspin1`.

Key result:

- Top Pt-C distance: about `1.8479 A`
- C-O distance: about `1.1500 A`
- top Pt-C `-ICOHP`: `0.058247`
- Pt-d/C-p contribution dominates the occupied Pt-C bonding contribution.

See `docs/pt111-co-cohp-results.md` for the full interpretation.

## ABACUS PW and LOBSTER Feasibility Note

The report `docs/diamond-abacus-pw-lobster-feasibility.md` records a diamond
test of ABACUS PW outputs against LOBSTER. The practical conclusion is:

- LOBSTER does not directly recognize ABACUS PW outputs.
- ABACUS `out_wfc_pw 2` can provide the plane-wave coefficients needed to build
  a LOBSTER Generic-style package, but this is only an interface path.
- LOBSTER's standard PW workflows require PAW pseudopotential data. QE itself
  supports NC, US, and PAW pseudopotentials, but LOBSTER's QE interface expects
  QE+PAW data rather than arbitrary QE PW output.
- The APNS `C.upf` used in the ABACUS PW probe is norm-conserving, so it lacks
  the PAW augmentation/projector/all-electron partial-wave data required for a
  scientifically strict LOBSTER COHP calculation.

This reinforces the intended scope of this repository: production use here is
ABACUS LCAO SCF -> `src/cohp.py` post-processing. ABACUS PW -> LOBSTER would be
a separate converter/exporter project with a consistent PAW data path.

## Development Checks

```bash
python -m py_compile src/cohp.py src/read_abacus_out.py scripts/pt111_co_workflow.py
python src/read_abacus_out.py
```
