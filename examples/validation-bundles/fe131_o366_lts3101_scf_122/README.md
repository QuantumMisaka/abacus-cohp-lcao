# Fe131-O366 COHP ABACUS SCF Case

Purpose: generate ABACUS LTS 3.10.1 LCAO H/S/WFC output for profiling `cohp.py` on Fe131 3d and O366 2p.

Key settings:
- KPT mesh: `1 2 2`
- Solver: `ks_solver cusolver`
- Parallelism: 4 V100 GPUs, 4 MPI ranks
- Basis: APNS precision orbitals
- Calculation: `scf`, no structure relaxation
- Bands: `nbands 3700`
- Smearing: Gaussian, `smearing_sigma 0.004`
- Spin: `nspin 2`
- Fe initial magnetic moments: existing Fe signs preserved, magnitudes set to `4.0`

Expected ABACUS outputs:
- `OUT.ABACUS/data-*-H`
- `OUT.ABACUS/data-*-S`
- `OUT.ABACUS/WFC_NAO_*`
- `OUT.ABACUS/running_scf.log`
