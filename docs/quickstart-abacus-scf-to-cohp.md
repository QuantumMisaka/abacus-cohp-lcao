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

## 2. Identify Global NAO Indices

`src/cohp.py` needs two orbital lists:

- `--atom-i-orbs`: global ABACUS NAO indices for atom or orbital group I.
- `--atom-j-orbs`: global ABACUS NAO indices for atom or orbital group J.

The indices are zero-based Python indices in the full ABACUS NAO order used by
the output matrices. For a two-atom Si example with 13 NAOs per atom, the first
Si atom is `0..12` and the second Si atom is `13..25`.

For orbital-resolved analysis, pass only the subset of orbital indices belonging
to the channel of interest, for example metal `d` orbitals against adsorbate `p`
orbitals. The Pt(111)-CO example stores such ranges in
`examples/pt111_co_top_nspin1/mapping.json`.

## 3. Run COHP Post-Processing

Install the Python requirements if needed:

```bash
python -m pip install -r requirements.txt
```

Run a total pair COHP:

```bash
env MPLBACKEND=Agg python src/cohp.py \
  --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2,3,4,5,6,7,8,9,10,11,12 \
  --atom-j-orbs 13,14,15,16,17,18,19,20,21,22,23,24,25 \
  --method COHP \
  --de 0.05 \
  --smooth-nstddev 4 \
  --invert \
  --output-prefix si_si_COHP
```

The command writes:

```text
si_si_COHP.dat
si_si_COHP.png
```

`--invert` plots `-COHP`, the common convention where positive occupied area is
usually interpreted as bonding contribution. Without `--invert`, the raw COHP
sign is written.

## 4. Spin-Polarized Output

For spin-polarized ABACUS output, run separate spin channels or the summed curve:

```bash
python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2 \
  --atom-j-orbs 100,101,102 \
  --spin up \
  --invert \
  --output-prefix pair_up

python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2 \
  --atom-j-orbs 100,101,102 \
  --spin down \
  --invert \
  --output-prefix pair_down

python src/cohp.py --out-dir /path/to/OUT.ABACUS \
  --atom-i-orbs 0,1,2 \
  --atom-j-orbs 100,101,102 \
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

## 6. Working Examples

- `examples/lts3101_lcao_si2`: minimal Si-Si validation result.
- `examples/pt111_co_top_nspin1`: Pt(111)-CO top-site result with total Pt-C and
  Pt-d/C-p channel examples.
- `examples/lts3101_lcao_si2/lobster_compare`: Si2 ABACUS LCAO-COHP vs
  fixed-geometry VASP+LOBSTER comparison artifacts.
- `examples/pt111_co_top_nspin2_lobster_compare`: spin-polarized Pt(111)-CO
  top-site ABACUS LCAO-COHP vs fixed-geometry VASP+LOBSTER comparison artifacts.

The bundled examples include final ABACUS input files, lightweight COHP outputs,
and selected VASP+LOBSTER comparison summaries. Heavy `OUT.ABACUS`,
VASP, and LOBSTER runtime directories are intentionally excluded.
