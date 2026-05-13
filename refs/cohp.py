import read_abacus_out as rao
import numpy as np
from pathlib import Path

"""COHP"""
def cal_COHPmatskIJ_e_ij(Hk, Sk, Ck, atomI_orbs, atomJ_orbs, mode="COHP"):
    """Calculate the energy-resolved COHP or COOP matrix for all bands, every selected orbital pair of the atom I and J, and a k-point.
    
    Args: 
    Hk: Hamiltonian matrix for a k-point, H(k),
    Sk: Overlap matrix for a k-point, S(k),
    Ck: Wavefunction set of one k-point, arranged as C(k)[i, n], where i is the local index and n is the band index.
    atomI_orbs: indices of the orbitals of the atom I,
    atomJ_orbs: indices of the orbitals of the atom J.
    mode: "COHP" or "COOP".

    Returns:
    np.ndarray: energy-resolved COHP or COOP matrix for a k-point.
    """
    nrows, ncols = Hk.shape
    nlocal, nbands = Ck.shape
    assert nrows == ncols
    assert nrows == nlocal
    # allocate memory for value
    value = [np.zeros((len(atomI_orbs), len(atomJ_orbs)), dtype=np.float64) for i in range(nbands)]

    for ib in range(nbands):
        for i_inval, i_inorb in enumerate(atomI_orbs):
            for j_inval, j_inorb in enumerate(atomJ_orbs):
                power = Sk[i_inorb, j_inorb] if mode == "COOP" else Hk[i_inorb, j_inorb]
                value[ib][i_inval, j_inval] = (Ck[i_inorb, ib].conj()*power*Ck[j_inorb, ib]).real
    return value

def cal_COHPvalskIJ_e(Hk, Sk, Ek, Ck, atomI_orbs, atomJ_orbs, mode="COHP"):
    """Compute the sum_{ij} over c*{Ii,n}(k)H{IiJj}(k)c{Jj,n}(k), 
    where I, J are atom indexes and i, j are orbitals indexes.
    n is the band index and this quantity corresponds to an eigenenergy.

    Args:
        Hk (np.ndarray): Hamiltonian matrix for a k-point.
        Ck (np.ndarray): Wavefunction set of one k-point.
        atomI_orbs (list): indices of the orbitals of the atom I.
        atomJ_orbs (list): indices of the orbitals of the atom J.

    Returns:
        tuple: eigenenergies and the sum_{ij} on c*{Ii,n}(k)H{IiJj}(k)c{Jj,n}(k) for each band.
    """
    
    nlocal, nband = Ck.shape
    matskIJ_ij_e = cal_COHPmatskIJ_e_ij(Hk, Sk, Ck, atomI_orbs, atomJ_orbs, mode=mode)
    # dimension assertation
    assert len(Ek) == nband
    assert len(matskIJ_ij_e) == nband
    for ib in range(nband):
        assert matskIJ_ij_e[ib].shape == (len(atomI_orbs), len(atomJ_orbs))
    # compute the sum_{ij} over c*{Ii,n}(k)H{IiJj}(k)c{Jj,n}(k), i.e., for atom-pair IJ
    valskIJ_e = [np.sum(matskIJ_ij_e[i]) for i in range(nband)]
    return Ek, valskIJ_e

def cal_COHPvalsIJ_e(Hks, Sks, Eks, Cks, wk = None, 
                     atomI_orbs = None, atomJ_orbs = None,
                     mode: str = "COHP"):
    nks = len(Hks)
    nbands = len(Eks[0])

    wk = [1/nks for i in range(nks)] if wk is None else wk

    Evals = []
    COHPvalsIJ_e = []
    for ik in range(nks):
        Evalsk, COHPvalskIJ_e = cal_COHPvalskIJ_e(Hk=Hks[ik], Sk=Sks[ik], Ek=Eks[ik], Ck=Cks[ik], 
                                                  atomI_orbs=atomI_orbs, atomJ_orbs=atomJ_orbs, mode=mode)
        # dimension assertation
        assert len(Evalsk) == nbands
        assert len(COHPvalskIJ_e) == nbands
        # add k-point weight
        COHPvalskIJ_e = [wk[ik]*COHPvalskIJ_e[i] for i in range(nbands)]
        # add to the global list
        Evals += Evalsk.tolist()
        COHPvalsIJ_e += COHPvalskIJ_e
    # dimension assertation
    assert len(Evals) == nks*nbands
    assert len(COHPvalsIJ_e) == nks*nbands
    # sort COHPvalsIJ_e according to the Evals
    Evals, COHPvalsIJ_e = zip(*sorted(zip(Evals, COHPvalsIJ_e)))
    # energy unit conversion
    Evals = np.array([rao.unit_conversion(e, "Ry", "eV") for e in Evals])
    return dos_integral(Evals, COHPvalsIJ_e)

"""pCOHP: why always positive?"""
def cal_pCOHPmatskIJ_e_ij(Ek, Sk, Ak, Ck, atomI_orbs, atomJ_orbs, mode="COHP"):
    """Extract the energy-resolved COHP or COOP matrix for a k-point.
    
    Args: 
    Hk: Hamiltonian matrix for a k-point, H(k),
    Ak: Overlap matrix for a k-point, A(k),
    Ck: Wavefunction set of one k-point, arranged as C(k)[i, n], where i is the local index and n is the band index.
    atomI_orbs: indices of the orbitals of the atom I,
    atomJ_orbs: indices of the orbitals of the atom J.
    mode: "COHP" or "COOP".

    Returns:
    np.ndarray: energy-resolved COHP or COOP matrix for a k-point.
    """
    Ek = np.diag(Ek)

    ndim1, ndim2 = Ek.shape
    nrows, ncols = Sk.shape
    nlocal, nbands = Ck.shape
    nlocal, nao = Ak.shape
    # dimension assertion
    assert ndim1 == ndim2
    assert ndim1 == nbands
    assert nrows == ncols
    assert nrows == nlocal
    assert nao == len(atomI_orbs) + len(atomJ_orbs)
    # allocate memory for value
    value = [np.zeros((len(atomI_orbs), len(atomJ_orbs)), dtype=np.float64) for i in range(nbands)]

    Tk = Ak.conj().T@Ck
    Eprimek = Tk@Ek@Tk.conj().T if mode == "COHP" else np.zeros((nao, nao))
    Sprimek = Tk@Tk.conj().T if mode == "COOP" else np.zeros((nao, nao)) # this is always to be positive

    for ib in range(nbands):
        for imu in range(len(atomI_orbs)):
            for inu in range(len(atomJ_orbs)):
                power = Sprimek[inu, imu] if mode == "COOP" else Eprimek[inu, imu]
                value[ib][imu, inu] = ((Tk.conj().T)[ib, inu]*power*Tk[imu, ib]).real
    return value

def cal_pCOHPvalskIJ_e(Hk, Sk, Ek, Ck, Ak, atomI_orbs, atomJ_orbs, mode="COHP"):
    """Compute the sum_{ij} over c*{Ii,n}(k)H{IiJj}(k)c{Jj,n}(k), 
    where I, J are atom indexes and i, j are orbitals indexes.
    n is the band index and this quantity corresponds to an eigenenergy.

    Args:
        Hk (np.ndarray): Hamiltonian matrix for a k-point.
        Ck (np.ndarray): Wavefunction set of one k-point.
        atomI_orbs (list): indices of the orbitals of the atom I.
        atomJ_orbs (list): indices of the orbitals of the atom J.

    Returns:
        tuple: eigenenergies and the sum_{ij} on c*{Ii,n}(k)H{IiJj}(k)c{Jj,n}(k) for each band.
    """
    
    nlocal, nband = Ck.shape
    assert len(Ek) == nband

    matskIJ_ij_e = cal_pCOHPmatskIJ_e_ij(Ek, Sk, Ak, Ck, atomI_orbs, atomJ_orbs, mode=mode)
    assert len(matskIJ_ij_e) == nband

    valskIJ_e = [np.sum(matskIJ_ij_e[i]) for i in range(nband)]
    return Ek, valskIJ_e

def cal_pCOHPvalsIJ_e(Hks, Sks, Eks, Cks, Aks, wk = None,
                      atomI_orbs = None, atomJ_orbs = None,
                      mode: str = "COHP"):
        nks = len(Hks)
        wk = [1/nks for i in range(nks)] if wk is None else wk
    
        Es = []
        pCOHPvalsIJ_e = []
        for ik in range(nks):
            Evalsk, pCOHPvalskIJ_e = cal_pCOHPvalskIJ_e(Hks[ik], Sks[ik], Eks[ik], Cks[ik], Aks[ik], atomI_orbs, atomJ_orbs, mode=mode)
            # for each band there will be value pair of (energy, pCOHP)
            assert len(Evalsk) == len(pCOHPvalskIJ_e)
            # add k-point weight
            pCOHPvalskIJ_e = [wk[ik]*pCOHPvalskIJ_e[i] for i in range(len(Evalsk))]
            # add to the global list
            Es += Evalsk.tolist()
            pCOHPvalsIJ_e += pCOHPvalskIJ_e
        # sort COHPvalsIJ_e according to the energy
        Es, pCOHPvalsIJ_e = zip(*sorted(zip(Es, pCOHPvalsIJ_e)))
        Es = np.array([rao.unit_conversion(e, "Ry", "eV") for e in Es])
        # discard the imaginary part
        pCOHPvalsIJ_e = np.array(pCOHPvalsIJ_e, dtype=np.float64)
        return dos_integral(Es, pCOHPvalsIJ_e)

"""utils"""
def dos_integral(x, y):
    x, y = np.array(x), np.array(y)
    indices = np.argsort(x)
    x, y = x[indices], y[indices]
    # get the unique energies
    unique_x = np.unique(x)
    # accumulate the occ with the same energy
    unique_y = np.array([np.sum(y[x == e]) for e in unique_x])

    return unique_x, unique_y

"""There is a bug in the following function"""
import matplotlib.pyplot as plt
def draw_COHP(x, y,
              testmethod = "COHP", 
              emin = None, emax = None,
              width = None,
              shift_toefermi: bool = True, efermi = None,
              invert_COHP: bool = True):
    # set x and y
    e_shifted = x - efermi if shift_toefermi and efermi is not None else x
    COHPvalsIJ_e = -y if invert_COHP else y
    # set ylim
    emin = min(e_shifted) if emin is None else emin
    emax = max(e_shifted) if emax is None else emax
    # set xlim
    width = max(abs(min(COHPvalsIJ_e)), abs(max(COHPvalsIJ_e))) if width is None else width

    plt.figure(figsize=(6, 18))
    plt.plot(COHPvalsIJ_e, e_shifted)
    #plt.ylim(emin, emax)
    plt.xlim(-width, width)

    plt.axvline(0, color='black', lw=0.5)
    plt.axhline(0, color='black', lw=0.5, linestyle='--', label='E_Fermi')
    plt.text(width*1.05, 0.5, r'$\epsilon_F$', fontsize=15)

    plt.fill_betweenx(e_shifted, COHPvalsIJ_e, 0, where=(e_shifted <= 0),
                      interpolate=True, alpha=0.3)
    plt.ylabel("Energy (eV)")

    xlabel = testmethod if not invert_COHP else "-"+testmethod
    plt.xlabel(xlabel, fontsize=15)
    # if width < 1e-3, use scientific notation
    if width <= 1e-3:
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    plt.savefig(f"{testmethod}.png")

"""about test"""
def test_selection(testcase, minimal_basis = False):
    if testcase == 1: # 2s2p1d
        path = "diamond/OUT.ABACUS"
        if minimal_basis:
            atomI_orbs = [0, 2, 3, 4]
            atomJ_orbs = [13, 15, 16, 17]
        else:
            atomI_orbs = list(range(0, 13))
            atomJ_orbs = list(range(13, 26))
    elif testcase == 2: # 3s3p3d2f, 3s3p2d
        path = "GaAs/OUT.ABACUS"
        if minimal_basis:
            atomI_orbs = [0, 3, 4, 5, 12, 13, 14, 15, 16]
            atomJ_orbs = [41, 44, 45, 46]
        else:
            atomI_orbs = list(range(0, 41))
            atomJ_orbs = list(range(41, 63))
    elif testcase == 3: # 6s3p2d, 3s3p2d
        path = "CsCl/OUT.ABACUS"
        if minimal_basis:
            atomI_orbs = [0, 6, 7, 8]
            atomJ_orbs = [25, 31, 32, 33]
        else:
            atomI_orbs = list(range(0, 25))
            atomJ_orbs = list(range(25, 47))
    else:
        raise ValueError("Invalid testcase")
    return path, atomI_orbs, atomJ_orbs

def test_initialize(testcase, minimal_basis):
    path, atomI_orbs, atomJ_orbs = test_selection(testcase=testcase, minimal_basis=minimal_basis)
    kptwts = rao.read_kpoints(path + "/kpoints", as_dict=False)[0][:, -1]
    nks = len(kptwts)

    Hks = [rao.read_mat_hs(path + f"/data-{ik}-H") for ik in range(nks)]
    Sks = [rao.read_mat_hs(path + f"/data-{ik}-S") for ik in range(nks)]
    temp = [rao.read_lowf(path + f"/LOWF_K_{ik+1}.txt") for ik in range(nks)]
    
    Cks, kvecs, Eks, occs = tuple(map(list, zip(*temp)))
    assert len(Cks) == nks
    assert len(kvecs) == nks
    assert len(Eks) == nks
    assert len(occs) == nks

    efermi = rao.read_etraj_fromlog(path + "/running_scf.log", term="fermi")[-1]
    return Hks, Sks, Cks, Eks, kptwts, efermi, atomI_orbs, atomJ_orbs

def _first_existing(candidates):
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None

def initialize_from_outdir(out_dir, atomI_orbs, atomJ_orbs, spin="sum"):
    """Initialize COHP post-processing from an ABACUS OUT.* directory."""
    out_dir = Path(out_dir)
    h_files = sorted(out_dir.glob("data-*-H"), key=lambda p: int(p.stem.split("-")[1]))
    s_files = sorted(out_dir.glob("data-*-S"), key=lambda p: int(p.stem.split("-")[1]))
    if len(h_files) == 0 or len(h_files) != len(s_files):
        raise FileNotFoundError(f"Cannot find matching data-*-H/data-*-S files in {out_dir}")

    nmatrix = len(h_files)
    Hks = [rao.read_mat_hs(str(path)) for path in h_files]
    Sks = [rao.read_mat_hs(str(path)) for path in s_files]

    wfc_files = []
    for ik in range(nmatrix):
        wfc = _first_existing([
            out_dir / f"WFC_NAO_K{ik + 1}.txt",
            out_dir / f"WFC_NAO_K{ik + 1}_ION1.txt",
            out_dir / f"LOWF_K_{ik + 1}.txt",
            out_dir / f"LOWF_K_{ik + 1}.dat",
            out_dir / f"WFC_NAO_GAMMA{ik + 1}.txt",
            out_dir / f"WFC_NAO_GAMMA{ik + 1}_ION1.txt",
            out_dir / f"LOWF_GAMMA_S{ik + 1}.dat",
        ])
        if wfc is None:
            raise FileNotFoundError(f"Cannot find wavefunction text file for k index {ik + 1} in {out_dir}")
        wfc_files.append(wfc)

    temp = [rao.read_lowf(str(path)) for path in wfc_files]
    Cks, kvecs, Eks, occs = tuple(map(list, zip(*temp)))

    kpoints_file = out_dir / "kpoints"
    if kpoints_file.exists():
        base_wts = rao.read_kpoints(str(kpoints_file), as_dict=False)[0][:, -1]
    else:
        base_wts = np.full(nmatrix, 1.0 / nmatrix)

    if nmatrix == len(base_wts):
        nspin = 1
        kptwts = base_wts
    elif nmatrix == 2 * len(base_wts):
        nspin = 2
        kptwts = np.concatenate([base_wts, base_wts])
    else:
        raise ValueError(
            f"Cannot map {nmatrix} H/S matrices onto {len(base_wts)} k-point weights in {out_dir}"
        )

    spin = spin.lower()
    if spin not in {"sum", "up", "down"}:
        raise ValueError("spin must be one of: sum, up, down")
    if nspin == 2 and spin in {"up", "down"}:
        nk = len(base_wts)
        selected = slice(0, nk) if spin == "up" else slice(nk, 2 * nk)
        Hks = Hks[selected]
        Sks = Sks[selected]
        Cks = Cks[selected]
        Eks = Eks[selected]
        kptwts = kptwts[selected]
    elif nspin == 1 and spin in {"up", "down"}:
        raise ValueError(f"Requested spin={spin}, but {out_dir} contains nspin=1 output")

    log_file = out_dir / "running_scf.log"
    efermi_values = rao.read_etraj_fromlog(str(log_file), term="fermi") if log_file.exists() else []
    efermi = efermi_values[-1] if len(efermi_values) else 0.0
    return Hks, Sks, Cks, Eks, kptwts, efermi, atomI_orbs, atomJ_orbs

def main(testcase, testmethod,
         de = 0.1, smooth = True, smooth_nstddev = 3,
         shift_toefermi = True, invert_COHP = False,
         emin = -10, emax = 10, width = 2,
         minimal_basis = False):
    # Initialize the test
    Hks, Sks, Cks, Eks, wks, efermi, atomI_orbs, atomJ_orbs = test_initialize(testcase=testcase, 
                                                                              minimal_basis=minimal_basis)
    # Compute the COHP/COOP/pCOHP/pCOOP for the IJ atom-pair
    if testmethod.startswith("pCO"):
        Aks = [Sk[:, atomI_orbs + atomJ_orbs] for Sk in Sks]
        e, COHPvalsIJ_e = cal_pCOHPvalsIJ_e(Hks=Hks, Sks=Sks, Eks=Eks, Cks=Cks, Aks=Aks,
                                            wk=wks,
                                            atomI_orbs=atomI_orbs, 
                                            atomJ_orbs=atomJ_orbs, 
                                            mode=testmethod[1:])
    elif testmethod.startswith("CO"):
        e, COHPvalsIJ_e = cal_COHPvalsIJ_e(Hks=Hks, Sks=Sks, Eks=Eks, Cks=Cks,
                                           wk=wks,
                                           atomI_orbs=atomI_orbs, 
                                           atomJ_orbs=atomJ_orbs, 
                                           mode=testmethod)
    else:
        raise ValueError("Invalid testmethod")
    # zero-padding COHP, necessary for recovering the original expression including \delta(\epsilon - \epsilon')
    e, COHPvalsIJ_e = rao.zero_padding(xmin=np.min(e)*1.1, 
                                       xmax=np.max(e)*1.1, 
                                       dx=de, 
                                       x=e, y=COHPvalsIJ_e)
    # Smoothing COHP, optional
    COHPvalsIJ_e = rao.Gauss_smoothing(x=e, 
                                       y=COHPvalsIJ_e, 
                                       sigma=smooth_nstddev*de, 
                                       normalize=False) if smooth else COHPvalsIJ_e
    # Draw the COHP
    draw_COHP(e, COHPvalsIJ_e, 
              testmethod=testmethod,
              emin=emin, emax=emax,
              width=width,
              shift_toefermi=shift_toefermi, efermi=efermi, 
              invert_COHP=invert_COHP)
    return e, COHPvalsIJ_e

def run_outdir(out_dir, atomI_orbs, atomJ_orbs, testmethod="COHP",
               de=0.1, smooth=True, smooth_nstddev=3,
               shift_toefermi=True, invert_COHP=False,
               emin=-10, emax=10, width=None, output_prefix=None, spin="sum"):
    Hks, Sks, Cks, Eks, wks, efermi, atomI_orbs, atomJ_orbs = initialize_from_outdir(
        out_dir=out_dir,
        atomI_orbs=atomI_orbs,
        atomJ_orbs=atomJ_orbs,
        spin=spin,
    )
    if testmethod.startswith("pCO"):
        Aks = [Sk[:, atomI_orbs + atomJ_orbs] for Sk in Sks]
        e, vals = cal_pCOHPvalsIJ_e(
            Hks=Hks, Sks=Sks, Eks=Eks, Cks=Cks, Aks=Aks, wk=wks,
            atomI_orbs=atomI_orbs, atomJ_orbs=atomJ_orbs, mode=testmethod[1:],
        )
    elif testmethod.startswith("CO"):
        e, vals = cal_COHPvalsIJ_e(
            Hks=Hks, Sks=Sks, Eks=Eks, Cks=Cks, wk=wks,
            atomI_orbs=atomI_orbs, atomJ_orbs=atomJ_orbs, mode=testmethod,
        )
    else:
        raise ValueError("Invalid testmethod")

    e, vals = rao.zero_padding(xmin=np.min(e) * 1.1, xmax=np.max(e) * 1.1, dx=de, x=e, y=vals)
    vals = rao.Gauss_smoothing(x=e, y=vals, sigma=smooth_nstddev * de, normalize=False) if smooth else vals
    if output_prefix:
        np.savetxt(f"{output_prefix}.dat", np.column_stack([e, vals]), header="Energy(eV) COHP")
        testmethod_for_plot = output_prefix
    else:
        testmethod_for_plot = testmethod
    draw_COHP(
        e, vals, testmethod=testmethod_for_plot, emin=emin, emax=emax, width=width,
        shift_toefermi=shift_toefermi, efermi=efermi, invert_COHP=invert_COHP,
    )
    return e, vals

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="ABACUS LCAO COHP post-processing")
    parser.add_argument("--out-dir", help="ABACUS OUT.* directory")
    parser.add_argument("--atom-i-orbs", help="Comma-separated orbital indices for atom/group I")
    parser.add_argument("--atom-j-orbs", help="Comma-separated orbital indices for atom/group J")
    parser.add_argument("--method", default="COHP", choices=["COHP", "COOP", "pCOHP", "pCOOP"])
    parser.add_argument("--de", type=float, default=0.1)
    parser.add_argument("--no-smooth", action="store_true")
    parser.add_argument("--smooth-nstddev", type=float, default=3)
    parser.add_argument("--emin", type=float, default=-10)
    parser.add_argument("--emax", type=float, default=10)
    parser.add_argument("--width", type=float, default=None)
    parser.add_argument("--invert", action="store_true")
    parser.add_argument("--output-prefix")
    parser.add_argument("--spin", default="sum", choices=["sum", "up", "down"])
    args = parser.parse_args()
    if args.out_dir:
        if not args.atom_i_orbs or not args.atom_j_orbs:
            raise ValueError("--atom-i-orbs and --atom-j-orbs are required with --out-dir")
        run_outdir(
            out_dir=args.out_dir,
            atomI_orbs=[int(x) for x in args.atom_i_orbs.split(",") if x],
            atomJ_orbs=[int(x) for x in args.atom_j_orbs.split(",") if x],
            testmethod=args.method,
            de=args.de,
            smooth=not args.no_smooth,
            smooth_nstddev=args.smooth_nstddev,
            emin=args.emin,
            emax=args.emax,
            width=args.width,
            invert_COHP=args.invert,
            output_prefix=args.output_prefix,
            spin=args.spin,
        )
        raise SystemExit(0)
    
    testcase = 1
    testmethod = "pCOHP"

    de = 0.05 # eV

    smooth = True
    smooth_nstddev = 5

    shift_toefermi = True
    invert_COHP = True

    emin = -10
    emax = 10
    width = 0.025

    minimal_basis = True

    import unittest
    class TestCOHP(unittest.TestCase):
        def test_dos_integral(self):
            x = [1, 2, 3, 4, 5, 6]
            y = [1, 2, 3, 4, 5, 6]
            x, y = dos_integral(x, y)
            self.assertEqual(x.tolist(), [1, 2, 3, 4, 5, 6])
            self.assertEqual(y.tolist(), [1, 2, 3, 4, 5, 6])

            x = [1, 2, 3, 3, 4, 4, 4, 5, 6, 6]
            y = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            x, y = dos_integral(x, y)
            self.assertEqual(x.tolist(), [1, 2, 3, 4, 5, 6])
            self.assertEqual(y.tolist(), [1, 1, 2, 3, 1, 2])

            x = [4, 4, 4, 2, 3, 1, 5, 8, 7, 6]
            y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            x, y = dos_integral(x, y)
            self.assertEqual(x.tolist(), [1, 2, 3, 4, 5, 6, 7, 8])
            self.assertEqual(y.tolist(), [6, 4, 5, 6, 7, 10, 9, 8])

            x = [-5, -10, 1, 6, -3, 0, 3, 5, 7, 8]
            y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            x, y = dos_integral(x, y)
            self.assertEqual(x.tolist(), [-10, -5, -3, 0, 1, 3, 5, 6, 7, 8])
            self.assertEqual(y.tolist(), [2, 1, 5, 6, 3, 7, 8, 4, 9, 10])
        def test_cal_COHPmatskIJ_e_ij(self):
            Hks = [rao.read_mat_hs("diamond/OUT.ABACUS/data-0-H")]
            Sks = [rao.read_mat_hs("diamond/OUT.ABACUS/data-0-S")]
            temp = [rao.read_lowf("diamond/OUT.ABACUS/LOWF_K_1.txt")]
            Cks, kvecs, Eks, occs = tuple(map(list, zip(*temp)))
            atomI_orbs = [0, 2, 3, 4]
            atomJ_orbs = [13, 15, 16, 17]
            value = cal_COHPmatskIJ_e_ij(Hk=Hks[0], Sk=Sks[0], Ck=Cks[0], atomI_orbs=atomI_orbs, atomJ_orbs=atomJ_orbs)
            
            nlocal, nbands = Cks[0].shape
            for ib in range(nbands):
                self.assertEqual(value[ib].shape, (len(atomI_orbs), len(atomJ_orbs)))
            value_totest_byhand = value[0]
            # first, is iorb, jorb, therefore H 0, 13
            print("Cks[0][0, 0] = ", Cks[0][0, 0])
            print("Hks[0][0, 13] = ", Hks[0][0, 13])
            print("Cks[0][13, 0] = ", Cks[0][13, 0])
            self.assertEqual(value_totest_byhand[0, 0], Cks[0][0, 0]*Hks[0][0, 13]*Cks[0][13, 0])
            self.assertEqual(value_totest_byhand[0, 1], Cks[0][0, 0]*Hks[0][0, 15]*Cks[0][15, 0])
            self.assertEqual(value_totest_byhand[0, 2], Cks[0][0, 0]*Hks[0][0, 16]*Cks[0][16, 0])
            self.assertEqual(value_totest_byhand[0, 3], Cks[0][0, 0]*Hks[0][0, 17]*Cks[0][17, 0])

    if "test" not in testmethod:
        main(testcase, testmethod,
            de = de, smooth = smooth, smooth_nstddev = smooth_nstddev,
            shift_toefermi = shift_toefermi, invert_COHP = invert_COHP,
            emin = emin, emax = emax, width = width,
            minimal_basis = minimal_basis)
    else:
        unittest.main()
