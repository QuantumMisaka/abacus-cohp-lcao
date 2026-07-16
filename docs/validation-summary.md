# ABACUS LCAO-COHP Validation Summary

This note records the current validation status for the ABACUS LCAO-COHP
post-processing workflow. On 2026-07-17 the implementation corrected the
Hamiltonian unit conversion, unordered Hermitian pair factor, and ICOHP
integration. Processed curves made with the earlier convention have therefore
been withdrawn. `examples/validation-bundles/` currently retains reproduction
inputs only; the reports in `docs/archive/` are historical and must not be used
as quantitative validation of the corrected implementation.

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

## Cases requiring regeneration

| Case | Public artifacts | Main purpose |
|---|---|---|
| Si2 | `examples/lts3101_lcao_si2` | Minimal LCAO SCF to Si-Si COHP validation |
| Pt(111)-CO | `examples/validation-bundles/pt111_co_top_nspin1` | Atom-pair and orbital-channel COHP on a surface adsorbate |
| Diamond | `examples/validation-bundles/diamond_cohp_compare` | C-C bonding comparison against VASP+LOBSTER |
| Ni(100)-CO | `examples/validation-bundles/ni100_co_top` | Spin-polarized metal-adsorbate and internal CO bond analysis |
| Fe/O large case | `examples/validation-bundles/fe131_o366_lts3101_scf_122` | Streaming parser performance and memory validation |

## Supported conclusions

- ABACUS LCAO-COHP gives chemically useful bonding and antibonding trends within
  a fixed ABACUS NAO setup.
- Absolute `-ICOHP` magnitudes are method- and representation-dependent.
- The reliable cross-code comparison is sign, relative trend, and normalized
  curve shape; direct equality of absolute `-ICOHP` is not a supported claim.
- Quantitative claims for Si-Si, C-C, Pt-C, Ni-C, and C-O remain pending until
  these cases are regenerated with the corrected implementation.

## Streaming Benchmark

The streaming and legacy readers agree in unit tests. The historical large-case
benchmark summary was withdrawn with the old processed curves and must be rerun
before publishing performance or numerical-equivalence claims.

## Archive

Detailed historical reports are kept in `docs/archive/`. They preserve run
context and exploratory notes, but the public user entrypoints are the README,
quickstart, method baseline, and this validation summary.
