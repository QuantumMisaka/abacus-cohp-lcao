# ABACUS SCF To COHP Quick Start

This guide shows the shortest reliable path from an ABACUS LCAO SCF calculation
to a COHP curve with this repository.

## 1. Prepare An ABACUS LCAO SCF

The COHP script does not run DFT. It post-processes an already completed ABACUS
LCAO SCF directory.

Add the following output block to the SCF `INPUT`:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

The four lines above are the minimal output settings for this workflow. In
particular, do not omit `out_app_flag 1`: the completed validation examples in
this repository use it, and it keeps ABACUS matrix/wavefunction outputs in the
layout expected by `src/read_abacus_out.py`.

ABACUS references for these output switches:

- `out_mat_hs`: https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/hs_matrix.html
- `out_wfc_lcao` and `out_app_flag`: https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html

After SCF finishes, check that `OUT.ABACUS` contains these files:

```text
OUT.ABACUS/data-0-H
OUT.ABACUS/data-0-S
OUT.ABACUS/WFC_NAO_K1.txt
OUT.ABACUS/kpoints
OUT.ABACUS/running_scf.log
```

For multiple k points, ABACUS writes one `data-*-H/S` pair and one
`WFC_NAO_K*.txt` file per k point.

## 2. Select Atoms And Orbital Channels

The recommended interface uses 1-based ABACUS atom indices and shell labels:

```bash
python src/cohp.py --out-dir /path/to/OUT.ABACUS --list-orbitals
```

The listing is inferred from `STRU`, `INPUT`, and `orbital_dir`. A metal-oxygen
`d-p` channel can then be requested as:

```bash
--atom-i-index 95 --atom-j-index 98 --atom-i-orbs 3d --atom-j-orbs 2p
```

Labels such as `3d`, `2p`, and `4s` are chemistry-friendly aliases. ABACUS NAOs
do not carry a strict principal quantum number in the matrix ordering, so the
script maps them to all NAOs of the matching angular-momentum channel on that
atom. Multiple channels can be comma-separated, for example
`--atom-i-orbs 3p,3d,4s`. If `--atom-*-orbs` is omitted with atom-index mode,
that atom defaults to `all`.

The legacy global-NAO mode remains available. If `--atom-i-index` and
`--atom-j-index` are omitted, `--atom-i-orbs` and `--atom-j-orbs` are interpreted
as zero-based global ABACUS NAO indices.

## 3. Run COHP Post-Processing

Install the Python requirements if needed:

```bash
python -m pip install -r requirements.txt
```

Run a total pair COHP:

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
  --output-prefix si_si_COHP
```

The command writes:

```text
si_si_COHP.dat
si_si_COHP_EminusEf.dat
si_si_COHP.meta.json
si_si_COHP.png
```

`--invert` plots `-COHP`, the common convention where positive occupied area is
usually interpreted as bonding contribution. Without `--invert`, the raw COHP
sign is written.
`si_si_COHP.dat` keeps the raw absolute-energy COHP data. `si_si_COHP_EminusEf.dat`
uses `E - E_Fermi`, matching the usual VASP+LOBSTER COHP energy reference.
`si_si_COHP.meta.json` records the Fermi energy, output files, method, spin, and
smoothing settings. The command-line log also prints `E_Fermi` and the generated
raw, shifted, metadata, and plot paths.

To keep only the absolute-energy output and plot without writing
`*_EminusEf.dat`, add:

```bash
--no-shift-to-efermi
```

## 4. Spin-Polarized Output

For spin-polarized ABACUS output, run separate spin channels or the summed curve:

```bash
python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-index 95 \
  --atom-j-index 98 \
  --atom-i-orbs 3d \
  --atom-j-orbs 2p \
  --spin up \
  --invert \
  --output-prefix pair_up

python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-index 95 \
  --atom-j-index 98 \
  --atom-i-orbs 3d \
  --atom-j-orbs 2p \
  --spin down \
  --invert \
  --output-prefix pair_down

python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-index 95 \
  --atom-j-index 98 \
  --atom-i-orbs 3d \
  --atom-j-orbs 2p \
  --spin sum \
  --invert \
  --output-prefix pair_sum
```

## 5. Interpret The Result

Use this script for ABACUS-internal, NAO-basis COHP analysis. The robust
interpretation is:

- sign and occupied/unoccupied energy distribution of `-COHP`;
- comparison between bonds computed with the same ABACUS version, pseudopotential
  set, NAO basis, and SCF settings;
- orbital-channel trends inside the same calculation.

Do not assume that ABACUS LCAO-COHP values are numerically identical to LOBSTER
pCOHP values. Tests in this project show consistent qualitative bonding trends
but method-dependent absolute ICOHP scales.

## 6. LOBSTER-like Empirical Scale

After generating a raw two-column `.dat` curve with `src/cohp.py`, users who
want a LOBSTER-like reading magnitude can apply the benchmark scale helper:

```bash
python scripts/scale_abacus_cohp_to_lobster.py si_si_COHP.dat \
  --preset Si-Si \
  --efermi 7.111283804 \
  --output-prefix si_si_lobster_like
```

For a file whose second column is already `-COHP` on an `E-E_F` axis:

```bash
python scripts/scale_abacus_cohp_to_lobster.py pair_minus_cohp.dat \
  --preset Pt-C \
  --input-convention minus-cohp \
  --efermi 0.0
```

If the input is a `*_EminusEf.dat` file, use `--efermi 0.0` to avoid subtracting
the Fermi energy twice.

Use `--list-presets` to inspect available channels. The script writes:

- `*_lobster_like.dat`: `energy_ev source_value minus_cohp lobster_like_minus_cohp`
- `*_lobster_like.json`: selected preset, scale factor, integrated `-ICOHP`,
  and the interpretation warning

The scale is empirical and intended for readability and plot comparison only.
It does not make ABACUS NAO-COHP and LOBSTER pCOHP strict numerical equivalents.

## 7. High-Performance Post-Processing

`src/cohp.py` defaults to the streaming COHP/COOP parser. It reads only the
selected H/S sub-block and the selected WFC rows, which keeps memory use low for
large ABACUS outputs.

Use this form for a large post-processing job:

```bash
env MPLBACKEND=Agg python src/cohp.py \
  --out-dir /path/to/OUT.ABACUS \
  --atom-i-index 95 \
  --atom-j-index 98 \
  --atom-i-orbs 3d \
  --atom-j-orbs 2p \
  --method COHP \
  --workers 4 \
  --invert \
  --output-prefix Fe131_O366
```

Use `--workers N` to parallelize k-point parsing. Keep `--legacy-full-read` for
debugging or old-style comparisons only.

## 8. Working Examples

- `examples/lts3101_lcao_si2`: minimal Si-Si validation result.
- `examples/data`: pseudopotential and numerical-orbital files needed by the
  bundled examples.
- `examples/validation-bundles`: processed validation artifacts for Si2,
  Pt(111)-CO, diamond, Ni(100)-CO, and a large Fe/O streaming benchmark.

The bundled examples include final ABACUS input files, lightweight COHP outputs,
and selected comparison summaries. Heavy `OUT.ABACUS`, VASP, and LOBSTER runtime
directories are intentionally excluded.
