# Ni(100)-CO Magnetism And COHP Comparability Note

## 1. Magnetic Moments

| Source | Total moment (Bohr/cell) | Absolute moment (Bohr/cell) | Notes |
|---|---:|---:|---|
| ABACUS efficiency SCF | 10.4264 | 11.0162 | From `OUT.ABACUS/running_scf.log` and `abacus.json`. |
| ABACUS precision SCF | 10.4553 | 11.7916 | From `OUT.ABACUS/running_scf.log` and `abacus.json`. |
| VASP relax final | 10.4881 | not printed | From final `number of electron ... magnetization` in `OUTCAR`. |
| VASP SCF final | 10.4887 | not printed | From final `number of electron ... magnetization` in `OUTCAR`. |

ABACUS and VASP therefore converge to essentially the same total spin state, around 10.4-10.5 Bohr per 2x2 four-layer slab plus CO cell.

The current VASP SCF used `LORBIT = 0`, so `OUTCAR` does not contain a native VASP `magnetization (x)` per-ion table. LOBSTER did write projected Mulliken/Loewdin gross populations from the VASP wavefunction; these provide a useful local-moment estimate, but should be labeled as LOBSTER-projected moments rather than raw VASP sphere moments.

| Atom | Element | Mulliken moment | Loewdin moment |
|---:|---|---:|---:|
| 1 | Ni | 0.720 | 0.720 |
| 2 | Ni | 0.720 | 0.720 |
| 3 | Ni | 0.720 | 0.720 |
| 4 | Ni | 0.730 | 0.720 |
| 5 | Ni | 0.710 | 0.700 |
| 6 | Ni | 0.650 | 0.650 |
| 7 | Ni | 0.650 | 0.660 |
| 8 | Ni | 0.700 | 0.690 |
| 9 | Ni | 0.660 | 0.660 |
| 10 | Ni | 0.660 | 0.660 |
| 11 | Ni | 0.660 | 0.660 |
| 12 | Ni | 0.660 | 0.660 |
| 13 | Ni | 0.720 | 0.730 |
| 14 | Ni | 0.710 | 0.710 |
| 15 | Ni | 0.700 | 0.700 |
| 16 | Ni(top, bonded to CO) | 0.180 | 0.200 |
| 17 | C | -0.040 | -0.040 |
| 18 | O | -0.030 | -0.030 |

The projected Ni moment sum is 10.55 Bohr by Mulliken population and 10.49 Bohr by Loewdin population, consistent with the VASP total moment. The top Ni atom bonded to CO is strongly quenched relative to other Ni atoms.

The current ABACUS outputs do not contain an atom-resolved magnetic moment table. They reliably provide total and absolute cell magnetization, but not per-atom moments in the present run artifacts.

## 2. ABACUS Efficiency vs Precision COHP

| Channel | Efficiency sum -ICOHP | Precision sum -ICOHP | Precision / Efficiency |
|---|---:|---:|---:|
| Ni-C total | 0.042415 | 0.057212 | 1.349 |
| Ni-d / C-p | 0.054402 | 0.047431 | 0.872 |
| C-O total | 0.502162 | 0.541192 | 1.078 |
| C-p / O-p | 0.367043 | 0.342366 | 0.933 |

The two ABACUS orbital sets agree well on the chemical hierarchy and signs: C-O is much stronger than Ni-C, and the occupied integrated `-ICOHP` is net bonding for both bonds. Quantitatively:

- C-O is stable: total differs by about 8%, and C-p/O-p differs by about 7%.
- Ni-d/C-p is stable as an orbital-resolved descriptor: about 13% difference and curve correlation about 0.97 after interpolation to the same energy grid.
- Ni-C total is less stable: precision gives about 35% larger integrated value, and total-curve shape correlation is low, indicating sensitivity to non-d/p channels, basis completeness, and the weak adsorption bond.

Conclusion: efficiency and precision ABACUS SCF wavefunctions give the same qualitative COHP chemistry, but the weak Ni-C total `-ICOHP` should not be treated as tightly converged quantitatively. For reporting, use precision as the primary ABACUS result and use efficiency as a robustness check.

## 3. ABACUS vs VASP+LOBSTER Semi-Quantitative Comparability

| Bond | ABACUS precision sum -ICOHP | LOBSTER sum -ICOHP | LOBSTER / ABACUS |
|---|---:|---:|---:|
| Ni-C | 0.057212 | 3.275980 | 57.26 |
| C-O | 0.541192 | 16.541540 | 30.57 |

The scale factor is not transferable between bonds: Ni-C needs a factor near 57, while C-O needs a factor near 31. This alone argues against a single semi-quantitative normalization between ABACUS and LOBSTER integrated `-ICOHP`.

Curve-shape comparison after interpolation and normalization also gives mixed results:

| Curve pair | Pearson correlation | Best-fit normalized RMSE |
|---|---:|---:|
| ABACUS precision Ni-C total vs LOBSTER Ni-C | 0.281 | 0.961 |
| ABACUS precision C-O total vs LOBSTER C-O | 0.879 | 0.473 |

C-O has meaningful normalized shape similarity, so qualitative energy-window comparison is reasonable for that strong internal bond. Ni-C does not: the weak adsorption bond is too projection- and basis-sensitive, and the normalized curve shape is not robustly comparable.

Conclusion: ABACUS and VASP+LOBSTER COHP are chemically comparable in sign, bonding hierarchy, and broad qualitative trends, but they are not generally semi-quantitatively comparable on a shared normalized numerical scale for this system. Any numerical calibration would need bond-specific and basis-specific benchmarking; a single global normalization is not justified.
