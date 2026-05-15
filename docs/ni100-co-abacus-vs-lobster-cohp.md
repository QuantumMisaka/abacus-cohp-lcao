# Ni(100)-CO ABACUS LCAO-COHP vs VASP+LOBSTER COHP Report

## Scope

This report compares COHP analysis for a spin-polarized Ni(100)-CO top-site adsorption model using:

- ABACUS LCAO-COHP with APNS efficiency orbitals for relaxation and both efficiency/precision orbitals for COHP SCF.
- VASP 6.4.2 SCF followed by LOBSTER 5.1.1 with the `pbeVaspFit2015` basis.

All plotted and tabulated COHP values below use the `-COHP` or `-pCOHP` convention: positive occupied integrated `-ICOHP` is interpreted as net bonding.

## Workflow Status

| Step | Job | Status | Output directory |
|---|---:|---|---|
| ABACUS relax, efficiency basis | 408431 | COMPLETED | `examples/ni100_co_top/abacus_relax_efficiency` |
| ABACUS SCF/COHP, efficiency basis | 408561 | COMPLETED | `examples/ni100_co_top/abacus_scf_efficiency` |
| ABACUS SCF/COHP, precision basis | 408562 | COMPLETED | `examples/ni100_co_top/abacus_scf_precision` |
| VASP relax, 16 V100 GPUs | 408715 | COMPLETED | raw output not bundled |
| VASP SCF, 16 V100 GPUs | 408799 | COMPLETED | raw output not bundled |
| LOBSTER COHP | 408938 | COMPLETED | processed curves bundled |

## Settings Audit

- Structure: Ni(100) 2x2 four-layer slab, CO adsorbed at the top site, bottom two Ni layers fixed.
- Spin: both ABACUS and VASP are spin-polarized; Ni initial magnetic moment is 1.
- ABACUS relax: efficiency APNS orbitals, force threshold 0.05 eV/Ang.
- ABACUS SCF/COHP: `symmetry -1`, `mixing_beta 0.4`, `mixing_ndim 20`, `scf_thr 1e-7`; efficiency and precision APNS orbital sets were both evaluated.
- VASP relax/SCF: `ENCUT=520`, `ISYM=-1`, `NELM=300`, `EDIFF=1E-6`, `ISPIN=2`, `MAGMOM=16*1.0 1*0.0 1*0.0`, `ISMEAR=0`, `SIGMA=0.05`.
- VASP relax-specific settings: `IBRION=2`, `POTIM=0.2`, `EDIFFG=-0.05`.
- VASP Slurm: 4 nodes on `4V100`, 16 GPUs, `mpirun -np ${SLURM_NTASKS} --map-by ${MAP_OPT} vasp_std`, with SAI rank mapping.
- LOBSTER: `COHPstartEnergy -12`, `COHPendEnergy 8`, `basisSet pbeVaspFit2015`; pairs are Ni16-C17 and C17-O18.

## Geometry And Magnetism

| Method | Ni-C distance (A) | C-O distance (A) | Total magnetic moment (Bohr) |
|---|---:|---:|---:|
| ABACUS relax + efficiency SCF | 1.7448 | 1.1564 | 10.4264 |
| ABACUS relax + precision SCF | 1.7448 | 1.1564 | 10.4553 |
| VASP relax + SCF / LOBSTER pair length | 1.7353 | 1.1631 | 10.4887 |

The ABACUS and VASP optimized local geometries agree closely: Ni-C differs by about 0.01 A and C-O differs by about 0.007 A. Both workflows preserve a clear spin-polarized Ni slab state with total moment near 10.5 Bohr.

## Integrated COHP Comparison

| Method | Pair/channel | up -ICOHP | down -ICOHP | sum -ICOHP |
|---|---|---:|---:|---:|
| ABACUS efficiency | Ni-C total | 0.020877 | 0.021740 | 0.042415 |
| ABACUS efficiency | Ni-d / C-p | 0.026132 | 0.028362 | 0.054402 |
| ABACUS efficiency | C-O total | 0.251662 | 0.250651 | 0.502162 |
| ABACUS efficiency | C-p / O-p | 0.185192 | 0.181870 | 0.367043 |
| ABACUS precision | Ni-C total | 0.028276 | 0.028952 | 0.057212 |
| ABACUS precision | Ni-d / C-p | 0.022689 | 0.024719 | 0.047431 |
| ABACUS precision | C-O total | 0.271206 | 0.269937 | 0.541192 |
| ABACUS precision | C-p / O-p | 0.172745 | 0.169606 | 0.342366 |
| LOBSTER pbeVaspFit2015 | Ni16-C17 | 1.612400 | 1.663580 | 3.275980 |
| LOBSTER pbeVaspFit2015 | C17-O18 | 8.392680 | 8.148860 | 16.541540 |

The absolute integrated values are not on the same numerical scale. LOBSTER reconstructs a local basis from PAW wavefunctions and reports much larger `-ICOHP` values; ABACUS LCAO-COHP uses the native numerical atomic-orbital Hamiltonian/overlap projection. The robust comparison is therefore qualitative: sign, relative bond strength, spin symmetry, and energy-window trends.

## Curve-Level Comparison

Relevant plotted data are in `examples/ni100_co_top/analysis/`:

- ABACUS precision Ni-C: `abacus_precision_Ni_C_total_sum.png`
- ABACUS precision C-O: `abacus_precision_C_O_total_sum.png`
- LOBSTER Ni-C: `lobster_pair1_minus_pcohp.png`
- LOBSTER C-O: `lobster_pair2_minus_pcohp.png`
- Normalized overlays: `overlay_abacus_precision_lobster_Ni_C.png`, `overlay_abacus_precision_lobster_C_O.png`

For Ni-C, both methods give a net bonding occupied integral under the `-COHP` convention. LOBSTER shows broad occupied bonding contributions in the deeper occupied window and clear antibonding weight close to and above the Fermi level. ABACUS gives the same net bonding sign, but the magnitude is much smaller and the energy distribution is smoother. The precision-basis ABACUS Ni-C total `-ICOHP` is larger than the efficiency-basis value, while the explicit Ni-d/C-p channel is similar but slightly smaller in the precision basis. This indicates that the Ni-C adsorption bond is weak and sensitive to basis/projector details, with Ni-d/C-p interactions still being the chemically relevant channel.

For C-O, both methods identify a much stronger internal CO bond than the Ni-C adsorption bond. LOBSTER gives a very large C-O `-ICOHP`, with strong occupied bonding intensity mainly in the lower occupied manifold and antibonding intensity near the upper occupied/unoccupied region. ABACUS also gives a stable positive C-O total `-ICOHP` in both efficiency and precision bases. The C-p/O-p channel accounts for most, but not all, of the ABACUS C-O bonding signal, consistent with the expected dominant p-p character of the CO bond.

The spin-up and spin-down integrated values are close for both Ni-C and C-O in ABACUS and LOBSTER. This is consistent with the total magnetism being mostly hosted by the Ni slab rather than producing a strongly spin-split CO internal bond.

## LOBSTER Reliability Notes

LOBSTER completed successfully and recovered 169.9953 of 170 electrons. The reported charge spillings are 2.62% for spin channel 1 and 3.67% for spin channel 2, which are acceptable for a qualitative COHP comparison. It also warned that 30 of 150 k-points could not be orthonormalized to `1.0E-5`; this is not automatically fatal, but it means the LOBSTER curves should be treated as a strong reference for chemical trends rather than an exact absolute benchmark.

LOBSTER also used recommended basis functions because no explicit orbital list was provided: C/O `2s 2p` and Ni `4s 3d`. PAW bands from 105 upward were ignored because there were more PAW bands than local basis functions.

## Conclusions

1. Both ABACUS LCAO-COHP and VASP+LOBSTER identify the same qualitative bonding hierarchy: C-O is much stronger than the Ni-C adsorption bond.
2. Both approaches give positive occupied `-ICOHP` for Ni-C and C-O, so the sign of the integrated bonding conclusion is consistent.
3. Both approaches show weak spin asymmetry in the COHP integrals despite a sizable total magnetic moment near 10.5 Bohr, so the slab is spin-polarized while the analyzed Ni-C-O bonds are not strongly spin-imbalanced.
4. Absolute `-ICOHP` values differ by large factors, so ABACUS LCAO-COHP should not be used as a drop-in quantitative substitute for LOBSTER without calibration.
5. ABACUS LCAO-COHP is chemically meaningful for qualitative bond analysis in this case: it captures the correct net bonding sign, the weak Ni-C versus strong C-O contrast, and the relevant orbital-channel trend. Its main value is fast, native post-processing of ABACUS LCAO results; its current limitation is quantitative comparability of integrated COHP magnitudes and detailed peak positions against LOBSTER.
