# Pt(111)-CO ABACUS LCAO-COHP vs VASP+LOBSTER COHP Report

## Workflow

- ABACUS: existing nspin=2 precision-basis LCAO final SCF and COHP output.
- VASP: fixed-geometry spin-polarized PBE PAW SCF generated from the same ABACUS final SCF STRU.
- LOBSTER: `pbeVaspFit2015` basis on the VASP SCF outputs.

## Quantitative Comparison

| Channel | Existing/full ABACUS -ICOHP | Common-window ABACUS -ICOHP | LOBSTER -ICOHP | LOBSTER / existing-full ABACUS | LOBSTER / common-window ABACUS | Pearson | Normalized RMSE | Best-fit scale |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Pt-C_total | 0.116498 | 0.079853 | 4.812040 | 41.305893 | 60.260967 | 0.598226 | 0.831035 | 0.024720 |
| C-O_total | 0.480546 | 0.189705 | 18.102810 | 37.671307 | 95.425869 | 0.766804 | 0.630763 | 0.019123 |

## Qualitative Comparison

- Curves are aligned to their own Fermi levels and plotted as `-COHP/-pCOHP`.
- Positive occupied area is interpreted as net bonding within each method's own representation.
- The fixed-geometry protocol isolates COHP representation/projection differences from structural relaxation differences.

## LOBSTER Quality

- Charge spilling: 2.3800 %

## Output Files

- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/pt111_co_top/nspin2/analysis/pt111_Pt-C_total_minus_cohp_overlay.png
- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/pt111_co_top/nspin2/analysis/pt111_C-O_total_minus_cohp_overlay.png
- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/pt111_co_top/nspin2/analysis/pt111_cohp_compare.csv
- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/pt111_co_top/nspin2/analysis/summary.json
