# Si2 ABACUS LCAO-COHP vs VASP+LOBSTER COHP Report

## Workflow

- ABACUS: existing LTS 3.10.x LCAO Si2 validation output.
- VASP: fixed-geometry PBE PAW SCF generated from the same ABACUS STRU geometry.
- LOBSTER: basis `Bunge` selected from completed results.

## Quantitative Comparison

| Quantity | Value |
|---|---:|
| Existing ABACUS Si-Si full-window -ICOHP | 0.061985 |
| ABACUS Si-Si common-window -ICOHP | 0.040358 |
| LOBSTER Si-Si -ICOHP | 4.622330 |
| LOBSTER / existing ABACUS -ICOHP | 74.572216 |
| LOBSTER / common-window ABACUS -ICOHP | 114.533745 |
| Pearson after E-grid alignment | 0.585982 |
| Best-fit normalized RMSE | 0.866088 |
| Best-fit LOBSTER scale to ABACUS | 0.037246 |

## Qualitative Comparison

- Curves are aligned to their own Fermi levels and plotted as `-COHP/-pCOHP`.
- Positive occupied area is treated as net bonding within each method's definition.
- Absolute -ICOHP magnitudes are reported but not interpreted as the same observable because ABACUS NAO COHP and LOBSTER projected pCOHP use different local representations.

## LOBSTER Quality

- Charge spilling: 1.4600 %
- pbeVaspFit2015 rerun recommended by 5% threshold: False

## Output Files

- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/lts3101_lcao_si2/analysis/si2_si_si_minus_cohp_overlay.png
- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/lts3101_lcao_si2/analysis/si2_si_si_minus_cohp_normalized_overlay.png
- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/lts3101_lcao_si2/analysis/si2_cohp_compare.csv
- /home/pku-jianghong/liuzhaoqing/work/sidereus/abacus-cohp/runs/lts3101_lcao_si2/analysis/summary.json
