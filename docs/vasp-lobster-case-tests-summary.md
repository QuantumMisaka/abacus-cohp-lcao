# VASP+LOBSTER Comparison Case Test Summary

This note summarizes the completed ABACUS LCAO-COHP vs VASP+LOBSTER comparison
cases. All curves are aligned to each method's own Fermi level before comparison.
Absolute integrated values are reported, but the robust cross-code comparison is
the sign, energy location, and normalized shape of the bonding/antibonding
features because ABACUS NAO COHP and LOBSTER pCOHP use different local
representations and normalization conventions.

## Completed Calculations

| Case | VASP job | LOBSTER job | LOBSTER basis | LOBSTER charge spilling | Status |
|---|---:|---:|---|---:|---|
| Si2 Si-Si | 440499 | 440568 | Bunge | 1.46 % | completed |
| Pt(111)-CO top, Pt-C and C-O | 440498 | 440620 | pbeVaspFit2015 | 2.38 % | completed |

## Quantitative Comparison

| Case/channel | ABACUS full-window -ICOHP | ABACUS common-window -ICOHP | LOBSTER -ICOHP | LOBSTER / full-window ABACUS | Pearson | Normalized RMSE | Best-fit LOBSTER scale |
|---|---:|---:|---:|---:|---:|---:|---:|
| Si2 Si-Si | 0.061985 | 0.040358 | 4.622330 | 74.572216 | 0.585982 | 0.866088 | 0.037246 |
| Pt(111)-CO Pt-C | 0.116498 | 0.079853 | 4.812040 | 41.305893 | 0.598226 | 0.831035 | 0.024720 |
| Pt(111)-CO C-O | 0.480546 | 0.189705 | 18.102810 | 37.671307 | 0.766804 | 0.630763 | 0.019123 |

## Qualitative Result

- All tested channels show net positive occupied `-COHP/-pCOHP` area in both
  ABACUS and LOBSTER, so the qualitative bonding assignment is consistent.
- LOBSTER gives a much larger absolute -ICOHP scale than the current ABACUS
  NAO-COHP postprocessor. This is treated as a projection/normalization
  difference rather than a direct bond-strength disagreement.
- Normalized overlays and best-fit scale factors are included for shape-level
  comparison after removing the dominant absolute-scale mismatch.
- Si2 Bunge spilling is below the 5% threshold, so no pbeVaspFit2015 rerun was
  required for that case.

## Reports And Artifacts

- `docs/si2-abacus-vs-vasp-lobster-cohp.md`
- `docs/pt111-co-abacus-vs-lobster-cohp.md`
- `runs/lts3101_lcao_si2/analysis/summary.json`
- `runs/pt111_co_top/nspin2/analysis/summary.json`

