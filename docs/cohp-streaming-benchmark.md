# COHP Streaming Benchmark

> Validation status (2026-07-17): historical processed summaries were
> withdrawn after correction of the COHP normalization and ICOHP conventions.
> Commands below describe how to regenerate them; previous numerical results
> must not be used as evidence for the corrected implementation.

This note records the cached ABACUS COHP post-processing benchmark for the
streaming COHP/COOP path in `src/cohp.py`.

## Cached Cases

The public `examples/` directory contains lightweight inputs, mappings, and
reference curves. The full benchmark was run in the development repository
against reusable ABACUS `OUT.ABACUS` caches:

- `runs/lts3101_lcao_si2/OUT.ABACUS`
- `runs/diamond_cohp_compare/abacus/scf/OUT.ABACUS`
- `runs/pt111_co_top/nspin1/final_scf/OUT.ABACUS`
- `runs/pt111_co_top/nspin2/final_scf/OUT.ABACUS`
- `runs/ni100_co_top/abacus/scf_efficiency/OUT.ABACUS`
- `runs/ni100_co_top/abacus/scf_precision/OUT.ABACUS`

The benchmark command was:

```bash
python scripts/benchmark_cohp_streaming.py --skip-large --workers 1 --output-label latest --summary-dir runs/cohp_streaming_benchmark
```

The generated development summary is:

- `runs/cohp_streaming_benchmark/latest/benchmark_summary.json`
- `runs/cohp_streaming_benchmark/latest/benchmark_summary.md`

All cached `runs/` cases completed successfully. The largest observed
fast-vs-reference curve difference was `8.882e-16` for diamond and `1.665e-16`
for the Ni(100)-CO efficiency case; all reported `minus_icohp_diff` values were
at numerical roundoff level.

## Large Fe/O Case

The large full-SCF case is documented in this public repository at:

- `examples/validation-bundles/fe131_o366_lts3101_scf_122`

It uses ABACUS LTS 3.10.1, APNS precision orbitals, `ks_solver cusolver`,
`nbands 3700`, Gaussian smearing with `sigma 0.004`, `nspin 2`, and a `1 2 2`
K-point mesh. The raw ABACUS output cache is about 20 GB and must not be
committed to git.

The large streaming benchmark command was:

```bash
python scripts/benchmark_cohp_streaming.py --case fe131_o366_full --workers 1 --output-label latest --summary-dir user_cases/fe131_o366_lts3101_scf_122/cohp_streaming_benchmark_summary
```

Result:

- channel: `Fe131_3d_O366_2p_sum`
- status: `ok`
- elapsed: `481.000 s`
- Python `ru_maxrss`: `86684 KB`
- Fermi energy: `2.1538541769 eV`
- publish a new summary only after rerunning with the corrected implementation

The previous large probe baseline for the legacy full-matrix path took about
12 minutes and required tens of GB of memory. The streaming full-SCF result is
therefore suitable as a performance regression case. The project should commit
only lightweight inputs, documentation, checksums, and benchmark summaries for
this case; the raw H/S/WFC text files should remain in external cache storage.
