# ABACUS LCAO COHP

Lightweight post-processing tools for ABACUS LCAO calculations. The main
command reads ABACUS Hamiltonian, overlap, wavefunction, eigenvalue, Fermi-level,
and k-point output files, then evaluates atom-pair or orbital-channel COHP/COOP
curves.

This repository targets **ABACUS numerical atomic orbital analysis**. Its COHP
values are ABACUS-NAO dependent and should not be treated as numerically
equivalent to LOBSTER pCOHP or ICOHP values.

## Contents

- `src/cohp.py`: user-facing COHP/COOP command-line post-processor.
- `src/read_abacus_out.py`: ABACUS output readers used by `cohp.py`.
- `scripts/scale_abacus_cohp_to_lobster.py`: optional empirical scaling helper
  for LOBSTER-like plot magnitudes.
- `docs/quickstart-abacus-scf-to-cohp.md`: full first-run guide.
- `docs/abacus-lcao-cohp-baseline.md`: method notes and interpretation limits.
- `docs/validation-summary.md`: compact summary of completed validation cases.
- `examples/lts3101_lcao_si2`: minimal Si2 example result.
- `examples/data`: pseudopotential and numerical-orbital subset used by bundled
  examples.
- `examples/validation-bundles`: processed validation artifacts; raw DFT output
  directories are intentionally excluded.

## Install

Use a Python environment with NumPy and Matplotlib:

```bash
python -m pip install -r requirements.txt
```

Optional example-generation workflows may require ASE and an ABACUS executable,
but the core post-processor starts from an already completed ABACUS LCAO SCF
output directory.

For experimental structure-generation or PW-to-LOBSTER diagnostic workflows,
install the optional dependencies separately:

```bash
python -m pip install -r requirements-experimental.txt
```

## Required ABACUS Output

The SCF `INPUT` must use an LCAO basis and request the matrix and wavefunction
files consumed by this post-processor:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

After SCF, `OUT.ABACUS` should contain at least:

```text
OUT.ABACUS/data-0-H
OUT.ABACUS/data-0-S
OUT.ABACUS/WFC_NAO_K1.txt
OUT.ABACUS/kpoints
OUT.ABACUS/running_scf.log
```

Multi-k calculations contain one `data-*-H/S` pair and one `WFC_NAO_K*.txt`
file per k point.

## Minimal COHP Usage

Inspect available atom shell channels:

```bash
python src/cohp.py --out-dir /path/to/OUT.ABACUS --list-orbitals
```

Run a pair COHP using 1-based ABACUS atom indices and shell labels:

```bash
env MPLBACKEND=Agg python src/cohp.py \
  --out-dir /path/to/OUT.ABACUS \
  --atom-i-index 1 \
  --atom-j-index 2 \
  --atom-i-orbs all \
  --atom-j-orbs all \
  --method COHP \
  --de 0.05 \
  --smooth-nstddev 4 \
  --invert \
  --output-prefix pair_COHP
```

For orbital-channel analysis, use labels such as `3d`, `2p`, or comma-separated
groups such as `3p,3d,4s`:

```bash
python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-index 95 \
  --atom-j-index 98 \
  --atom-i-orbs 3d \
  --atom-j-orbs 2p \
  --spin sum \
  --invert \
  --output-prefix Fe_O_d_p
```

If atom indices are omitted, `--atom-i-orbs` and `--atom-j-orbs` keep the legacy
zero-based global NAO-index mode, for example `0,1,2`.

## Outputs

The command writes:

- `pair_COHP.dat`: raw two-column absolute-energy and COHP data.
- `pair_COHP_EminusEf.dat`: two-column `E - E_Fermi` and COHP data, written by
  default for easier comparison with common COHP plotting conventions.
- `pair_COHP.meta.json`: Fermi energy, output file paths, method, spin, and
  smoothing settings, plus occupied `ICOHP/-ICOHP`.
- `pair_COHP.png`: plotted curve. With `--invert`, the figure follows the
  common `-COHP` convention where positive occupied area is usually interpreted
  as bonding contribution. The plot labels `-ICOHP` by default.

Add `--no-shift-to-efermi` to suppress `*_EminusEf.dat` and plot on the absolute
energy axis.
Add `--no-icohp-label` to hide the plot annotation.

The reported `ICOHP/-ICOHP` is integrated before Gaussian smoothing is applied
for plotting; it is not inferred from the smoothed filled area. Existing
two-column curves can be integrated directly:

```bash
python scripts/integrate_icohp.py pair_COHP.dat --efermi 7.111283804
python scripts/integrate_icohp.py pair_COHP_EminusEf.dat
```

## Performance

`src/cohp.py` uses the streaming COHP/COOP path by default. It reads only the
selected H/S sub-block and selected WFC rows, which keeps memory use low for
large ABACUS outputs.

Use `--workers N` to parallelize k-point parsing. Use `--legacy-full-read` only
for debugging or old full-matrix comparisons.

## Optional LOBSTER-Like Scale

After generating a raw two-column `.dat` curve, apply an empirical scale preset
if you want a LOBSTER-like plot magnitude:

```bash
python scripts/scale_abacus_cohp_to_lobster.py pair_COHP.dat \
  --preset Si-Si \
  --efermi 0.0 \
  --output-prefix pair_lobster_like
```

Use `--list-presets` to inspect available channels. This helper is for
readability and validation plots only; it does not make ABACUS NAO-COHP and
LOBSTER pCOHP the same observable.

## Scientific Boundary

Use ABACUS LCAO-COHP for trends inside a fixed ABACUS setup: same ABACUS
version, pseudopotentials, NAO basis, and SCF settings. The robust interpretation
is the sign, occupied/unoccupied energy distribution, relative bond trends, and
orbital-channel contributions. Do not compare absolute ICOHP values directly
against LOBSTER without an explicit benchmark and clear caveat.

## Documentation

- Start with `docs/quickstart-abacus-scf-to-cohp.md`.
- Read `docs/abacus-lcao-cohp-baseline.md` for method positioning and limits.
- Read `docs/validation-summary.md` for completed Si2, Pt(111)-CO, diamond,
  Ni(100)-CO, and Fe/O streaming benchmark summaries.
- Historical case reports and exploratory notes are kept in `docs/archive/`.
