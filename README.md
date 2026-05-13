# ABACUS LCAO COHP

This repository contains a lightweight ABACUS LCAO-COHP post-processing workflow.
It extracts Hamiltonian, overlap, wavefunction, eigenvalue, Fermi-level, and k-point
weight information from ABACUS LCAO SCF outputs, then evaluates atom-pair or
orbital-group COHP curves.

The implementation is intended for ABACUS numerical atomic orbital analysis. Its
COHP values are ABACUS-NAO dependent and should not be assumed numerically
equivalent to LOBSTER pCOHP.

## Contents

- `src/read_abacus_out.py`: ABACUS output readers for `data-*-H/S`, `WFC_NAO_K*.txt`,
  `kpoints`, `running_scf.log`, and helper utilities.
- `src/cohp.py`: COHP/COOP and pCOHP-like post-processing routines plus a CLI.
- `scripts/pt111_co_workflow.py`: reproducible Pt(111)-CO example workflow for SAI
  Slurm + ABACUS LTS v3.10.1.
- `docs/`: method notes, validation notes, and Pt(111)-CO result report.
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
```

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

## Development Checks

```bash
python -m py_compile src/cohp.py src/read_abacus_out.py scripts/pt111_co_workflow.py
python src/read_abacus_out.py
```
