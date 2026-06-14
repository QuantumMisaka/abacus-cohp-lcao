from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys

import read_abacus_out as rao
import numpy as np

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

def cal_COHPvalskIJ_e_selected(block, Ek, Ck_selected, n_i, mode="COHP"):
    """Compute COHP/COOP for one k point from selected matrix/WFC blocks."""
    c_i = Ck_selected[:n_i, :]
    c_j = Ck_selected[n_i:, :]
    vals = np.einsum("ib,ij,jb->b", c_i.conj(), block, c_j, optimize=True).real
    return Ek, vals

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

def _sorted_matrix_files(out_dir, suffix):
    return sorted(out_dir.glob(f"data-*-{suffix}"), key=lambda p: int(p.stem.split("-")[1]))

def _wfc_file_for_index(out_dir, ik):
    return _first_existing([
        out_dir / f"WFC_NAO_K{ik + 1}.txt",
        out_dir / f"WFC_NAO_K{ik + 1}_ION1.txt",
        out_dir / f"LOWF_K_{ik + 1}.txt",
        out_dir / f"LOWF_K_{ik + 1}.dat",
        out_dir / f"WFC_NAO_GAMMA{ik + 1}.txt",
        out_dir / f"WFC_NAO_GAMMA{ik + 1}_ION1.txt",
        out_dir / f"LOWF_GAMMA_S{ik + 1}.dat",
    ])

def _matrix_file_selection(out_dir, spin="sum"):
    out_dir = Path(out_dir)
    h_files = _sorted_matrix_files(out_dir, "H")
    s_files = _sorted_matrix_files(out_dir, "S")
    if len(h_files) == 0 or len(h_files) != len(s_files):
        raise FileNotFoundError(f"Cannot find matching data-*-H/data-*-S files in {out_dir}")

    nmatrix = len(h_files)
    kpoints_file = out_dir / "kpoints"
    if kpoints_file.exists():
        base_wts = rao.read_kpoints(str(kpoints_file), as_dict=False)[0][:, -1]
    else:
        base_wts = np.full(nmatrix, 1.0 / nmatrix)

    if nmatrix == len(base_wts):
        nspin = 1
        kptwts = np.asarray(base_wts, dtype=float)
    elif nmatrix == 2 * len(base_wts):
        nspin = 2
        kptwts = np.concatenate([base_wts, base_wts]).astype(float)
    else:
        raise ValueError(
            f"Cannot map {nmatrix} H/S matrices onto {len(base_wts)} k-point weights in {out_dir}"
        )

    spin = spin.lower()
    if spin not in {"sum", "up", "down"}:
        raise ValueError("spin must be one of: sum, up, down")
    if nspin == 2 and spin == "up":
        indices = range(0, len(base_wts))
    elif nspin == 2 and spin == "down":
        indices = range(len(base_wts), nmatrix)
    elif nspin == 1 and spin in {"up", "down"}:
        raise ValueError(f"Requested spin={spin}, but {out_dir} contains nspin=1 output")
    else:
        indices = range(nmatrix)

    selected = []
    for ik in indices:
        wfc = _wfc_file_for_index(out_dir, ik)
        if wfc is None:
            raise FileNotFoundError(f"Cannot find wavefunction text file for k index {ik + 1} in {out_dir}")
        selected.append((ik, h_files[ik], s_files[ik], wfc, float(kptwts[ik])))
    return selected

def _cohp_streaming_one(item, atomI_orbs, atomJ_orbs, mode="COHP"):
    _ik, h_file, s_file, wfc_file, weight = item
    selected_orbs = list(atomI_orbs) + list(atomJ_orbs)
    mat_file = s_file if mode == "COOP" else h_file
    block = rao.read_mat_hs_submatrix(mat_file, atomI_orbs, atomJ_orbs)
    c_selected, _kvec, ek, _occ = rao.read_lowf_selected(wfc_file, selected_orbs)
    ek, vals = cal_COHPvalskIJ_e_selected(block, ek, c_selected, len(atomI_orbs), mode=mode)
    return ek, vals * weight

def cal_COHPvalsIJ_e_streaming(out_dir, atomI_orbs, atomJ_orbs, wk=None, mode="COHP", spin="sum", workers=1):
    """Stream COHP/COOP from selected H/S and WFC blocks without full dense matrices."""
    if mode not in {"COHP", "COOP"}:
        raise ValueError("Streaming path supports COHP and COOP only")
    items = _matrix_file_selection(out_dir, spin=spin)
    if wk is not None:
        if len(wk) != len(items):
            raise ValueError(f"Expected {len(items)} k weights, got {len(wk)}")
        items = [(ik, h, s, wfc, float(wk_i)) for (ik, h, s, wfc, _), wk_i in zip(items, wk)]

    if workers is None or workers <= 1 or len(items) <= 1:
        pieces = [_cohp_streaming_one(item, atomI_orbs, atomJ_orbs, mode=mode) for item in items]
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(_cohp_streaming_one, item, list(atomI_orbs), list(atomJ_orbs), mode)
                for item in items
            ]
            pieces = [future.result() for future in futures]

    Evals = []
    COHPvalsIJ_e = []
    for ek, vals in pieces:
        Evals += ek.tolist()
        COHPvalsIJ_e += vals.tolist()
    Evals, COHPvalsIJ_e = zip(*sorted(zip(Evals, COHPvalsIJ_e)))
    Evals = np.array([rao.unit_conversion(e, "Ry", "eV") for e in Evals])
    COHPvalsIJ_e = np.array(COHPvalsIJ_e, dtype=np.float64)
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

def _trapezoid(y: np.ndarray, x: np.ndarray) -> float:
    integrator = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(integrator(y, x))


def integrate_icohp(energy, values, efermi=0.0, input_convention="cohp"):
    """Integrate occupied COHP and report both ICOHP and -ICOHP conventions."""
    if input_convention not in {"cohp", "minus-cohp"}:
        raise ValueError("input_convention must be 'cohp' or 'minus-cohp'")

    energy = np.asarray(energy, dtype=float)
    values = np.asarray(values, dtype=float)
    if energy.shape != values.shape:
        raise ValueError("energy and values must have the same shape")

    indices = np.argsort(energy)
    energy = energy[indices]
    values = values[indices]
    cohp_values = -values if input_convention == "minus-cohp" else values
    occupied = energy <= float(efermi)
    occupied_energy = energy[occupied]
    occupied_cohp = cohp_values[occupied]

    if len(occupied_energy) < 2:
        icohp = float("nan")
    else:
        icohp = _trapezoid(occupied_cohp, occupied_energy)

    return {
        "icohp": float(icohp),
        "minus_icohp": float(-icohp),
        "efermi_ev": float(efermi),
        "input_convention": input_convention,
        "integration": "trapezoid",
        "energy_window_ev": {
            "min": float(occupied_energy[0]) if len(occupied_energy) else float("nan"),
            "max": float(occupied_energy[-1]) if len(occupied_energy) else float("nan"),
            "occupied_condition": "energy <= efermi",
        },
        "points_total": int(len(energy)),
        "points_occupied": int(len(occupied_energy)),
        "uses_plot_smoothing": False,
    }


def _icohp_annotation(icohp):
    return f"-ICOHP = {icohp['minus_icohp']:.6f} eV"


"""There is a bug in the following function"""
import matplotlib.pyplot as plt
def _visible_energy_window(x, y, emin=None, emax=None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    emin = float(np.min(x)) if emin is None else float(emin)
    emax = float(np.max(x)) if emax is None else float(emax)
    mask = (x >= emin) & (x <= emax)
    if not np.any(mask):
        raise ValueError(f"No COHP data points in energy window [{emin}, {emax}] eV")
    return x[mask], y[mask], emin, emax


def draw_COHP(x, y,
              testmethod = "COHP", 
              emin = None, emax = None,
              width = None,
              shift_toefermi: bool = True, efermi = None,
              invert_COHP: bool = True,
              annotation: str | None = None):
    # set x and y
    e_shifted = x - efermi if shift_toefermi and efermi is not None else x
    COHPvalsIJ_e = -y if invert_COHP else y
    e_plot, cohp_plot, emin, emax = _visible_energy_window(e_shifted, COHPvalsIJ_e, emin=emin, emax=emax)
    # set xlim
    width = max(abs(min(cohp_plot)), abs(max(cohp_plot))) if width is None else width
    if width == 0:
        width = 1.0

    plt.figure(figsize=(6, 18))
    plt.plot(cohp_plot, e_plot)
    plt.ylim(emin, emax)
    plt.xlim(-width, width)

    plt.axvline(0, color='black', lw=0.5)
    if emin <= 0 <= emax:
        plt.axhline(0, color='black', lw=0.5, linestyle='--', label='E_Fermi')
        plt.text(width*1.05, 0.5, r'$\epsilon_F$', fontsize=15)

    plt.fill_betweenx(e_plot, cohp_plot, 0, where=(e_plot <= 0),
                      interpolate=True, alpha=0.3)
    plt.ylabel("Energy (eV)")

    xlabel = testmethod if not invert_COHP else "-"+testmethod
    plt.xlabel(xlabel, fontsize=15)
    if annotation:
        plt.text(
            0.03,
            0.97,
            annotation,
            transform=plt.gca().transAxes,
            va="top",
            ha="left",
            fontsize=12,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.75, "edgecolor": "0.8"},
        )
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

L_MULTIPLICITY = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}
SHELL_ORDER = ["s", "p", "d", "f", "g"]


@dataclass
class OrbitalAtom:
    atom_index: int
    symbol: str
    orbital_start: int
    orbital_stop: int
    shells: dict


@dataclass
class OrbitalMap:
    atoms: list
    total_orbitals: int
    stru_path: Path
    input_path: Path | None
    orbital_dir: Path | None


@dataclass
class OrbitalSelection:
    atom_index: int
    symbol: str
    selector: str
    normalized_selectors: list
    indices: list


def _strip_comment(line):
    return line.split("#", 1)[0].strip()


def _section_lines(lines, section_name):
    section_headers = {
        "ATOMIC_SPECIES",
        "NUMERICAL_ORBITAL",
        "LATTICE_CONSTANT",
        "LATTICE_VECTORS",
        "ATOMIC_POSITIONS",
    }
    start = None
    for idx, line in enumerate(lines):
        if _strip_comment(line).upper() == section_name:
            start = idx + 1
            break
    if start is None:
        return []

    out = []
    for line in lines[start:]:
        cleaned = _strip_comment(line)
        if not cleaned:
            if out:
                break
            continue
        if cleaned.upper() in section_headers:
            break
        out.append(cleaned)
    return out


def _atomic_position_lines(lines):
    section_headers = {
        "ATOMIC_SPECIES",
        "NUMERICAL_ORBITAL",
        "LATTICE_CONSTANT",
        "LATTICE_VECTORS",
        "ATOMIC_POSITIONS",
    }
    start = None
    for idx, line in enumerate(lines):
        if _strip_comment(line).upper() == "ATOMIC_POSITIONS":
            start = idx + 1
            break
    if start is None:
        return []

    out = []
    for line in lines[start:]:
        cleaned = _strip_comment(line)
        if not cleaned:
            continue
        if cleaned.upper() in section_headers:
            break
        out.append(cleaned)
    return out


def parse_input_orbital_dir(input_path):
    """Read orbital_dir from an ABACUS INPUT file."""
    if input_path is None:
        return None
    input_path = Path(input_path)
    for line in input_path.read_text().splitlines():
        fields = _strip_comment(line).split()
        if len(fields) >= 2 and fields[0].lower() == "orbital_dir":
            path = Path(fields[1]).expanduser()
            return path if path.is_absolute() else (input_path.parent / path).resolve()
    return None


def parse_stru_metadata(stru_path):
    """Parse species order, orbital filenames, and atom symbols from an ABACUS STRU."""
    stru_path = Path(stru_path)
    lines = stru_path.read_text().splitlines()

    species = []
    for line in _section_lines(lines, "ATOMIC_SPECIES"):
        fields = line.split()
        if fields:
            species.append(fields[0])
    if not species:
        raise ValueError(f"Cannot parse ATOMIC_SPECIES from {stru_path}")

    orbital_files = _section_lines(lines, "NUMERICAL_ORBITAL")
    if len(orbital_files) != len(species):
        raise ValueError(
            f"NUMERICAL_ORBITAL in {stru_path} has {len(orbital_files)} entries, "
            f"but ATOMIC_SPECIES has {len(species)} species"
        )

    atom_position_lines = _atomic_position_lines(lines)
    if len(atom_position_lines) < 2:
        raise ValueError(f"Cannot parse ATOMIC_POSITIONS from {stru_path}")

    atoms = []
    cursor = 1
    while cursor < len(atom_position_lines):
        symbol = atom_position_lines[cursor].split()[0]
        if cursor + 2 >= len(atom_position_lines):
            raise ValueError(f"Incomplete ATOMIC_POSITIONS block for {symbol} in {stru_path}")
        natom = int(float(atom_position_lines[cursor + 2].split()[0]))
        atoms.extend([symbol] * natom)
        cursor += 3 + natom

    orbital_by_symbol = dict(zip(species, orbital_files))
    return species, orbital_by_symbol, atoms


def _shell_counts_from_orbital_file(orbital_path):
    text = Path(orbital_path).read_text(errors="ignore")
    shells = {}
    for shell in SHELL_ORDER:
        pattern = rf"Number\s+of\s+{shell.upper()}orbital\s*-->\s*(\d+)"
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            shells[shell] = int(match.group(1))
    return shells


def _shell_counts_from_filename(orbital_path):
    name = Path(orbital_path).name.lower()
    matches = re.findall(r"(\d+)([spdfg])", name)
    return {shell: int(count) for count, shell in matches}


def parse_orbital_shells(orbital_path):
    """Return zeta counts per angular-momentum shell from an ABACUS .orb file."""
    orbital_path = Path(orbital_path)
    if orbital_path.exists():
        shells = _shell_counts_from_orbital_file(orbital_path)
        if shells:
            return shells
    shells = _shell_counts_from_filename(orbital_path)
    if shells:
        return shells
    raise ValueError(f"Cannot determine orbital shells from {orbital_path}")


def _resolve_orbital_path(entry, stru_path, orbital_dir):
    path = Path(entry).expanduser()
    candidates = []
    if path.is_absolute():
        candidates.append(path)
    else:
        candidates.append((Path(stru_path).parent / path).resolve())
        if orbital_dir is not None:
            candidates.append((Path(orbital_dir) / path.name).resolve())
            candidates.append((Path(orbital_dir) / path).resolve())
    existing = _first_existing(candidates)
    return existing if existing is not None else candidates[0]


def _shell_ranges(start, shells):
    ranges = {}
    cursor = start
    for shell in SHELL_ORDER:
        if shell not in shells:
            continue
        width = shells[shell] * L_MULTIPLICITY[shell]
        ranges[shell] = list(range(cursor, cursor + width))
        cursor += width
    ranges["all"] = list(range(start, cursor))
    return ranges


def build_orbital_map(stru_path, input_path=None, orbital_dir=None):
    """Build a 1-based atom to global-NAO shell map from STRU/INPUT/orbital files."""
    stru_path = Path(stru_path)
    input_path = Path(input_path) if input_path is not None else None
    orbital_dir = Path(orbital_dir).expanduser() if orbital_dir is not None else parse_input_orbital_dir(input_path)
    if orbital_dir is not None and not orbital_dir.is_absolute():
        orbital_dir = (stru_path.parent / orbital_dir).resolve()

    _, orbital_by_symbol, atom_symbols = parse_stru_metadata(stru_path)
    shells_by_symbol = {}
    for symbol, orbital_entry in orbital_by_symbol.items():
        orbital_path = _resolve_orbital_path(orbital_entry, stru_path, orbital_dir)
        shells_by_symbol[symbol] = parse_orbital_shells(orbital_path)

    atoms = []
    cursor = 0
    for idx, symbol in enumerate(atom_symbols, start=1):
        if symbol not in shells_by_symbol:
            raise ValueError(f"No NUMERICAL_ORBITAL entry for atom {idx} symbol {symbol}")
        ranges = _shell_ranges(cursor, shells_by_symbol[symbol])
        atoms.append(
            OrbitalAtom(
                atom_index=idx,
                symbol=symbol,
                orbital_start=cursor,
                orbital_stop=cursor + len(ranges["all"]),
                shells=ranges,
            )
        )
        cursor += len(ranges["all"])
    return OrbitalMap(
        atoms=atoms,
        total_orbitals=cursor,
        stru_path=stru_path,
        input_path=input_path,
        orbital_dir=orbital_dir,
    )


def parse_global_orbital_indices(text):
    """Parse a comma-separated 0-based global NAO index list."""
    indices = [int(token.strip()) for token in text.split(",") if token.strip()]
    if not indices:
        raise ValueError("Orbital index list is empty")
    if any(index < 0 for index in indices):
        raise ValueError("Global NAO indices must be non-negative")
    return indices


def _selector_tokens(selector):
    return [token.strip().lower() for token in selector.split(",") if token.strip()]


def _selector_shell(token):
    if token == "all":
        return "all"
    match = re.fullmatch(r"(?:\d+)?([spdfg])", token)
    return match.group(1) if match else None


def resolve_atom_orbitals(orbital_map, atom_index, selector="all"):
    """Resolve a 1-based atom index and shell selector to global NAO indices."""
    if atom_index < 1 or atom_index > len(orbital_map.atoms):
        raise ValueError(f"atom index {atom_index} is outside 1..{len(orbital_map.atoms)}")

    tokens = _selector_tokens(selector or "all")
    numeric = [re.fullmatch(r"\d+", token) is not None for token in tokens]
    if any(numeric):
        if not all(numeric):
            raise ValueError("Do not mix global NAO indices with shell labels in --atom-*-orbs")
        raise ValueError("Use global NAO indices without --atom-*-index")

    atom = orbital_map.atoms[atom_index - 1]
    normalized = []
    selected = set()
    for token in tokens:
        shell = _selector_shell(token)
        if shell is None:
            raise ValueError(f"Invalid orbital selector '{token}'")
        if shell not in atom.shells:
            available = ", ".join([s for s in SHELL_ORDER + ["all"] if s in atom.shells])
            raise ValueError(
                f"Atom {atom_index} ({atom.symbol}) has no {shell} shell; "
                f"available shells: {available}"
            )
        if shell not in normalized:
            normalized.append(shell)
        selected.update(atom.shells[shell])

    return OrbitalSelection(
        atom_index=atom_index,
        symbol=atom.symbol,
        selector=selector,
        normalized_selectors=normalized,
        indices=sorted(selected),
    )


def _infer_companion_file(out_dir, explicit_path, filename):
    if explicit_path:
        return Path(explicit_path)
    out_dir = Path(out_dir)
    candidates = [out_dir / filename, out_dir.parent / filename]
    return _first_existing(candidates)


def resolve_cli_orbitals(out_dir, atom_i_orbs, atom_j_orbs,
                         atom_i_index=None, atom_j_index=None,
                         stru_path=None, input_path=None, orbital_dir=None):
    """Resolve CLI orbital arguments, preserving legacy global-index mode."""
    if atom_i_index is None and atom_j_index is None:
        if not atom_i_orbs or not atom_j_orbs:
            raise ValueError("--atom-i-orbs and --atom-j-orbs are required with --out-dir")
        return parse_global_orbital_indices(atom_i_orbs), parse_global_orbital_indices(atom_j_orbs), None
    if atom_i_index is None or atom_j_index is None:
        raise ValueError("--atom-i-index and --atom-j-index must be used together")

    stru_path = _infer_companion_file(out_dir, stru_path, "STRU")
    input_path = _infer_companion_file(out_dir, input_path, "INPUT")
    if stru_path is None:
        raise FileNotFoundError("Cannot find STRU; pass --stru explicitly")

    orbital_map = build_orbital_map(stru_path=stru_path, input_path=input_path, orbital_dir=orbital_dir)
    selection_i = resolve_atom_orbitals(orbital_map, atom_i_index, atom_i_orbs or "all")
    selection_j = resolve_atom_orbitals(orbital_map, atom_j_index, atom_j_orbs or "all")
    return selection_i.indices, selection_j.indices, (selection_i, selection_j, orbital_map)


def print_orbital_resolution(info):
    if info is None:
        return
    selection_i, selection_j, orbital_map = info
    for label, selection in [("I", selection_i), ("J", selection_j)]:
        preview = ",".join(str(index) for index in selection.indices[:12])
        if len(selection.indices) > 12:
            preview += ",..."
        shells = ",".join(selection.normalized_selectors)
        print(
            f"atom {label}: {selection.symbol} #{selection.atom_index}, "
            f"selector {selection.selector} -> {shells}, "
            f"{len(selection.indices)} NAOs, global indices {preview}",
            file=sys.stderr,
        )
    print(f"orbital map: {len(orbital_map.atoms)} atoms, {orbital_map.total_orbitals} NAOs", file=sys.stderr)


def print_orbital_map(orbital_map):
    print("atom_index symbol orbital_start orbital_stop shells")
    for atom in orbital_map.atoms:
        shells = ",".join([shell for shell in SHELL_ORDER + ["all"] if shell in atom.shells])
        print(f"{atom.atom_index} {atom.symbol} {atom.orbital_start} {atom.orbital_stop} {shells}")

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


def _cohp_output_paths(output_prefix):
    prefix = Path(output_prefix)
    raw_path = prefix.with_suffix(".dat")
    shifted_path = prefix.with_name(f"{prefix.stem}_EminusEf").with_suffix(".dat")
    metadata_path = prefix.with_suffix(".meta.json")
    return raw_path, shifted_path, metadata_path


def _write_cohp_outputs(output_prefix, energy, values, efermi, *,
                        shift_toefermi=True, testmethod="COHP", spin="sum",
                        de=0.1, smooth=True, smooth_nstddev=3, icohp=None):
    raw_path, shifted_path, metadata_path = _cohp_output_paths(output_prefix)
    np.savetxt(
        raw_path,
        np.column_stack([energy, values]),
        header="Energy(eV) COHP shifted_by_efermi=false",
    )

    files = {
        "raw": {
            "path": str(raw_path),
            "energy_column": "Energy(eV)",
            "shifted_by_efermi": False,
        }
    }
    if shift_toefermi:
        shifted_energy = np.asarray(energy, dtype=float) - efermi
        np.savetxt(
            shifted_path,
            np.column_stack([shifted_energy, values]),
            header="Energy-E_Fermi(eV) COHP shifted_by_efermi=true",
        )
        files["shifted"] = {
            "path": str(shifted_path),
            "energy_column": "Energy-E_Fermi(eV)",
            "shifted_by_efermi": True,
        }
    files["metadata"] = {
        "path": str(metadata_path),
    }

    metadata = {
        "efermi_ev": float(efermi),
        "shift_toefermi": bool(shift_toefermi),
        "method": testmethod,
        "spin": spin,
        "de_ev": float(de),
        "smooth": bool(smooth),
        "smooth_nstddev": float(smooth_nstddev),
        "files": files,
    }
    if icohp is not None:
        metadata["icohp"] = icohp
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def _print_cohp_output_summary(metadata, plot_path):
    print(f"E_Fermi = {metadata['efermi_ev']:.6f} eV", file=sys.stderr)
    print(f"raw COHP: {metadata['files']['raw']['path']}", file=sys.stderr)
    if "shifted" in metadata["files"]:
        print(f"E-E_Fermi COHP: {metadata['files']['shifted']['path']}", file=sys.stderr)
    if "icohp" in metadata:
        print(f"ICOHP = {metadata['icohp']['icohp']:.6f} eV", file=sys.stderr)
        print(f"-ICOHP = {metadata['icohp']['minus_icohp']:.6f} eV", file=sys.stderr)
    print(f"metadata: {metadata['files']['metadata']['path']}", file=sys.stderr)
    print(f"plot: {plot_path}", file=sys.stderr)

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
               emin=-10, emax=10, width=None, output_prefix=None, spin="sum",
               workers=1, legacy_full_read=False, icohp_label=True):
    out_dir = Path(out_dir)
    log_file = out_dir / "running_scf.log"
    efermi_values = rao.read_etraj_fromlog(str(log_file), term="fermi") if log_file.exists() else []
    efermi = efermi_values[-1] if len(efermi_values) else 0.0

    if testmethod.startswith("pCO") or legacy_full_read:
        Hks, Sks, Cks, Eks, wks, efermi, atomI_orbs, atomJ_orbs = initialize_from_outdir(
            out_dir=out_dir,
            atomI_orbs=atomI_orbs,
            atomJ_orbs=atomJ_orbs,
            spin=spin,
        )
        nlocal = Hks[0].shape[0]
        requested = atomI_orbs + atomJ_orbs
        if requested and max(requested) >= nlocal:
            raise ValueError(f"Requested global NAO index {max(requested)} but ABACUS output has {nlocal} orbitals")

    if testmethod.startswith("pCO"):
        Aks = [Sk[:, atomI_orbs + atomJ_orbs] for Sk in Sks]
        e, vals = cal_pCOHPvalsIJ_e(
            Hks=Hks, Sks=Sks, Eks=Eks, Cks=Cks, Aks=Aks, wk=wks,
            atomI_orbs=atomI_orbs, atomJ_orbs=atomJ_orbs, mode=testmethod[1:],
        )
    elif testmethod.startswith("CO") and legacy_full_read:
        e, vals = cal_COHPvalsIJ_e(
            Hks=Hks, Sks=Sks, Eks=Eks, Cks=Cks, wk=wks,
            atomI_orbs=atomI_orbs, atomJ_orbs=atomJ_orbs, mode=testmethod,
        )
    elif testmethod.startswith("CO"):
        e, vals = cal_COHPvalsIJ_e_streaming(
            out_dir=out_dir,
            atomI_orbs=atomI_orbs,
            atomJ_orbs=atomJ_orbs,
            mode=testmethod,
            spin=spin,
            workers=workers,
        )
    else:
        raise ValueError("Invalid testmethod")

    e, vals_unsmoothed = rao.zero_padding(xmin=np.min(e) * 1.1, xmax=np.max(e) * 1.1, dx=de, x=e, y=vals)
    icohp = integrate_icohp(e, vals_unsmoothed, efermi=efermi, input_convention="cohp")
    vals = (
        rao.Gauss_smoothing(x=e, y=vals_unsmoothed, sigma=smooth_nstddev * de, normalize=False)
        if smooth
        else vals_unsmoothed
    )
    output_prefix = output_prefix or testmethod
    metadata = _write_cohp_outputs(
        output_prefix,
        e,
        vals,
        efermi,
        shift_toefermi=shift_toefermi,
        testmethod=testmethod,
        spin=spin,
        de=de,
        smooth=smooth,
        smooth_nstddev=smooth_nstddev,
        icohp=icohp,
    )
    testmethod_for_plot = output_prefix
    draw_COHP(
        e, vals, testmethod=testmethod_for_plot, emin=emin, emax=emax, width=width,
        shift_toefermi=shift_toefermi, efermi=efermi, invert_COHP=invert_COHP,
        annotation=_icohp_annotation(icohp) if icohp_label else None,
    )
    _print_cohp_output_summary(metadata, Path(f"{testmethod_for_plot}.png"))
    return e, vals

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="ABACUS LCAO COHP post-processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python src/cohp.py --out-dir OUT.ABACUS --atom-i-orbs 0,1,2 --atom-j-orbs 100,101\n"
            "  python src/cohp.py --out-dir OUT.ABACUS --atom-i-index 95 --atom-j-index 98 "
            "--atom-i-orbs 3d --atom-j-orbs 2p\n"
            "  python src/cohp.py --out-dir OUT.ABACUS --atom-i-index 95 --atom-j-index 98 "
            "--atom-i-orbs 3d --atom-j-orbs 2p --workers 4\n"
            "  python src/cohp.py --out-dir OUT.ABACUS --atom-i-index 95 --atom-j-index 98 "
            "--atom-i-orbs 3d --atom-j-orbs 2p --legacy-full-read\n"
            "  python src/cohp.py --out-dir OUT.ABACUS --list-orbitals\n"
            "\nPerformance tips:\n"
            "  - COHP/COOP uses the streaming parser by default and only reads the selected H/S block and WFC rows.\n"
            "  - Use --workers N to parallelize k-point parsing when the output contains many files.\n"
            "  - Use --legacy-full-read only when you need the older full-matrix path for debugging or comparison.\n"
        ),
    )
    parser.add_argument("--out-dir", help="ABACUS OUT.* directory")
    parser.add_argument(
        "--atom-i-orbs",
        help="Comma-separated global NAO indices, or shell labels with --atom-i-index, e.g. 3d or 3p,3d,4s",
    )
    parser.add_argument(
        "--atom-j-orbs",
        help="Comma-separated global NAO indices, or shell labels with --atom-j-index, e.g. 2p or all",
    )
    parser.add_argument("--atom-i-index", type=int, help="1-based atom index for atom/group I")
    parser.add_argument("--atom-j-index", type=int, help="1-based atom index for atom/group J")
    parser.add_argument("--stru", help="ABACUS STRU path for atom-index shell selection")
    parser.add_argument("--input", dest="input_path", help="ABACUS INPUT path for orbital_dir discovery")
    parser.add_argument("--orbital-dir", help="Directory containing ABACUS numerical orbital files")
    parser.add_argument("--list-orbitals", action="store_true", help="List atom shell channels and exit")
    parser.add_argument("--method", default="COHP", choices=["COHP", "COOP", "pCOHP", "pCOOP"])
    parser.add_argument("--de", type=float, default=0.1)
    parser.add_argument("--no-smooth", action="store_true")
    parser.add_argument("--smooth-nstddev", type=float, default=3)
    parser.add_argument("--emin", type=float, default=-10)
    parser.add_argument("--emax", type=float, default=10)
    parser.add_argument("--width", type=float, default=None)
    parser.add_argument("--invert", action="store_true")
    parser.add_argument(
        "--no-shift-to-efermi",
        action="store_false",
        dest="shift_toefermi",
        help="Do not write E-E_Fermi data or shift the automatic plot energy axis",
    )
    parser.add_argument("--output-prefix")
    parser.add_argument("--spin", default="sum", choices=["sum", "up", "down"])
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Parallel worker count for the streaming COHP/COOP parser; 1 keeps it serial",
    )
    parser.add_argument(
        "--legacy-full-read",
        action="store_true",
        help="Force the legacy full-matrix reader instead of the default streaming COHP/COOP path",
    )
    parser.add_argument(
        "--no-icohp-label",
        action="store_false",
        dest="icohp_label",
        help="Do not annotate the automatic COHP plot with the occupied -ICOHP value",
    )
    parser.set_defaults(shift_toefermi=True)
    parser.set_defaults(icohp_label=True)
    args = parser.parse_args()
    if args.out_dir:
        if args.list_orbitals:
            stru_path = _infer_companion_file(args.out_dir, args.stru, "STRU")
            input_path = _infer_companion_file(args.out_dir, args.input_path, "INPUT")
            if stru_path is None:
                raise FileNotFoundError("Cannot find STRU; pass --stru explicitly")
            orbital_map = build_orbital_map(
                stru_path=stru_path,
                input_path=input_path,
                orbital_dir=args.orbital_dir,
            )
            print_orbital_map(orbital_map)
            raise SystemExit(0)
        atomI_orbs, atomJ_orbs, resolution = resolve_cli_orbitals(
            out_dir=args.out_dir,
            atom_i_orbs=args.atom_i_orbs,
            atom_j_orbs=args.atom_j_orbs,
            atom_i_index=args.atom_i_index,
            atom_j_index=args.atom_j_index,
            stru_path=args.stru,
            input_path=args.input_path,
            orbital_dir=args.orbital_dir,
        )
        print_orbital_resolution(resolution)
        run_outdir(
            out_dir=args.out_dir,
            atomI_orbs=atomI_orbs,
            atomJ_orbs=atomJ_orbs,
            testmethod=args.method,
            de=args.de,
            smooth=not args.no_smooth,
            smooth_nstddev=args.smooth_nstddev,
            emin=args.emin,
            emax=args.emax,
            width=args.width,
            invert_COHP=args.invert,
            shift_toefermi=args.shift_toefermi,
            output_prefix=args.output_prefix,
            spin=args.spin,
            workers=args.workers,
            legacy_full_read=args.legacy_full_read,
            icohp_label=args.icohp_label,
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
