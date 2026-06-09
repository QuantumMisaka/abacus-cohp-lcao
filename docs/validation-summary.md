# ABACUS LCAO-COHP Validation Summary

This note summarizes the public validation evidence for the ABACUS LCAO-COHP
post-processing workflow in this repository. Detailed run notes are archived in
`docs/archive/`; processed lightweight artifacts are under
`examples/validation-bundles/`.

## Scope

The validated production path is:

```text
ABACUS LCAO SCF output -> src/cohp.py -> COHP/COOP curves
```

The required ABACUS output block is:

```text
basis_type lcao
out_mat_hs 1 8
out_wfc_lcao 1
out_app_flag 1
```

## Completed Cases

| Case | Public artifacts | Main purpose |
|---|---|---|
| Si2 | `examples/lts3101_lcao_si2`, `examples/validation-bundles/si2_lobster_compare` | Minimal LCAO SCF to Si-Si COHP validation |
| Pt(111)-CO | `examples/validation-bundles/pt111_co_top_nspin1`, `examples/validation-bundles/pt111_co_top_nspin2_lobster_compare` | Atom-pair and orbital-channel COHP on a surface adsorbate |
| Diamond | `examples/validation-bundles/diamond_cohp_compare` | C-C bonding comparison against VASP+LOBSTER processed curves |
| Ni(100)-CO | `examples/validation-bundles/ni100_co_top` | Spin-polarized metal-adsorbate and internal CO bond analysis |
| Fe/O large case | `examples/validation-bundles/fe131_o366_lts3101_scf_122` | Streaming parser performance and memory validation |

## Practical Conclusions

- ABACUS LCAO-COHP gives chemically useful bonding and antibonding trends within
  a fixed ABACUS NAO setup.
- Si-Si, C-C, Pt-C, Ni-C, and C-O channels show qualitatively consistent
  bonding assignments against processed VASP+LOBSTER references.
- Absolute `-ICOHP` magnitudes are method- and representation-dependent.
  LOBSTER values are often tens to more than one hundred times larger than the
  current ABACUS NAO-COHP values in these benchmarks.
- The reliable cross-code comparison is sign, relative trend, and normalized
  curve shape; direct equality of absolute `-ICOHP` is not a supported claim.

## Streaming Benchmark

The default streaming parser was checked against cached reference curves for
small and medium validation cases. Reported fast-vs-reference curve differences
were at numerical roundoff level.

The large Fe/O case used ABACUS LTS 3.10.1 precision orbitals, `nspin 2`,
`nbands 3700`, and a `1 2 2` k-point mesh. The raw ABACUS output is not bundled,
but the public summary records:

- channel: `Fe131_3d_O366_2p_sum`
- elapsed time: `481.000 s`
- Python `ru_maxrss`: `86684 KB`
- Fermi energy: `2.1538541769 eV`

The corresponding processed summary is in
`examples/validation-bundles/fe131_o366_lts3101_scf_122/cohp_streaming_benchmark_summary/latest/`.

## Archive

Detailed historical reports are kept in `docs/archive/`. They preserve run
context and exploratory notes, but the public user entrypoints are the README,
quickstart, method baseline, and this validation summary.
