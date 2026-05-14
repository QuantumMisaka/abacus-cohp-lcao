#!/usr/bin/env python3
"""Probe whether ABACUS PW output can be consumed by LOBSTER COHP.

This workflow intentionally separates three questions:
1. Can ABACUS PW generate the plane-wave data needed by LOBSTER Generic?
2. Can LOBSTER directly detect ABACUS outputs? (expected: no)
3. Can a diagnostic Generic package be assembled from ABACUS WAVEFUNC*.dat?
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT_DEFAULT = REPO_ROOT / "runs" / "diamond_abacus_pw_lobster"
REPORT = REPO_ROOT / "docs" / "diamond-abacus-pw-lobster-feasibility.md"

PP_ORB = Path.home() / "PP_ORB"
APNS_C_UPF = PP_ORB / "PP" / "C.upf"
LOBSTER_BIN = REPO_ROOT / "tools" / "lobster-5.1.1" / "lobster-5.1.1"
LOBSTER_QE_C_PAW = REPO_ROOT / "tools" / "lobster-5.1.1" / "QE" / "diamond" / "C.pbe-n-kjpaw_psl.0.1.UPF"

BOHR_TO_ANG = 0.529177210903
ANG_TO_BOHR = 1.0 / BOHR_TO_ANG
RY_TO_EV = 13.605693009
INITIAL_A = 3.567
K_MESH = (7, 7, 7)
NBANDS = 16
ECUT_RY = 100
ABACUS_SUFFIX = "ABACUS"


@dataclass
class WfcPw:
    ik: int
    nk: int
    k_cart: np.ndarray
    weight: float
    ng: int
    nbands: int
    ecut_ry: float
    lat0: float
    tpiba: float
    recip: np.ndarray
    g_indices: np.ndarray
    coeffs: np.ndarray


def run(cmd: list[str], cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=check, text=True, capture_output=True)


def primitive_cell(a: float = INITIAL_A) -> np.ndarray:
    return np.array(
        [[0.0, a / 2.0, a / 2.0], [a / 2.0, 0.0, a / 2.0], [a / 2.0, a / 2.0, 0.0]],
        dtype=float,
    )


def layout(root: Path) -> dict[str, Path]:
    return {
        "root": root,
        "abacus_pw": root / "abacus_pw_scf",
        "direct": root / "lobster" / "direct_abacus_outputs",
        "generic": root / "lobster" / "generic_abacus_pw",
        "analysis": root / "analysis",
        "jobs": root / "jobs.json",
        "state": root / "workflow_state.json",
    }


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def read_json(path: Path, default: Any) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text())


def write_abacus_stru(path: Path) -> None:
    cell = primitive_cell()
    frac = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]], dtype=float)
    lines = [
        "ATOMIC_SPECIES",
        "C 12.011 C.upf",
        "",
        "LATTICE_CONSTANT",
        f"{ANG_TO_BOHR:.15f}",
        "",
        "LATTICE_VECTORS",
        *[" ".join(f"{x:.16f}" for x in row) for row in cell],
        "",
        "ATOMIC_POSITIONS",
        "Direct",
        "",
        "C",
        "0.0",
        "2",
        *[" ".join(f"{x:.16f}" for x in row) + " 1 1 1" for row in frac],
    ]
    path.write_text("\n".join(lines) + "\n")


def write_kpt(path: Path) -> None:
    path.write_text("K_POINTS\n0\nGamma\n" + f"{K_MESH[0]} {K_MESH[1]} {K_MESH[2]} 0 0 0\n")


def write_abacus_input(path: Path, solver: str) -> None:
    lines = [
        "INPUT_PARAMETERS",
        "calculation scf",
        "basis_type pw",
        f"pseudo_dir {APNS_C_UPF.parent}",
        f"suffix {ABACUS_SUFFIX}",
        f"ecutwfc {ECUT_RY}",
        f"nbands {NBANDS}",
        "nspin 1",
        "symmetry -1",
        "device gpu",
        f"ks_solver {solver}",
        "scf_thr 1e-8",
        "scf_nmax 120",
        "mixing_beta 0.4",
        "smearing_method fixed",
        "out_wfc_pw 2",
        "out_chg 1",
        "out_band 1",
        "out_bandgap 1",
    ]
    if solver == "dav_subspace":
        lines.append("pw_diag_ndim 2")
    path.write_text("\n".join(lines) + "\n")


def write_abacus_sbatch(path: Path, solver: str) -> None:
    lines = [
        "#!/bin/bash",
        "#SBATCH --job-name=dia-ab-pw",
        "#SBATCH --partition=4V100",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --gpus-per-node=1",
        "#SBATCH --qos=rush-1o2gpu",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
        "",
        "set -euo pipefail",
        "module purge",
        "module load abacus/LTSv3.10.1-sm70-auto",
        "export OMP_NUM_THREADS=8",
        "export OMP_PROC_BIND=spread",
        "export OMP_PLACES=cores",
        "nvidia-smi dmon -s pucvmte -o T > nvdmon_job-${SLURM_JOB_ID}.log &",
        "abacus > abacus.stdout 2> abacus.stderr",
        "exit",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_lobsterin(path: Path, direct: bool) -> None:
    lines = [
        "COHPstartEnergy -22",
        "COHPendEnergy 18",
        "basisSet Bunge",
        "includeOrbitals sp",
        "cohpbetween atom 1 and atom 2",
    ]
    if direct:
        lines.insert(0, "# Expected to fail: LOBSTER has no native ABACUS detector.")
    path.write_text("\n".join(lines) + "\n")


def write_lobster_sbatch(path: Path, job_name: str) -> None:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --partition=CPU-MISC",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=16",
        "#SBATCH --qos=huge-cpu",
        "#SBATCH --output=slurm-%j.out",
        "#SBATCH --error=slurm-%j.err",
        "",
        "set -euo pipefail",
        "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}",
        f"'{LOBSTER_BIN}' > lobster.out 2> lobster.err",
        "exit",
    ]
    path.write_text("\n".join(lines) + "\n")


def prepare(root: Path, solver: str) -> None:
    if not APNS_C_UPF.exists():
        raise FileNotFoundError(APNS_C_UPF)
    if not LOBSTER_BIN.exists():
        raise FileNotFoundError(LOBSTER_BIN)
    paths = layout(root)
    for key in ("abacus_pw", "direct", "generic", "analysis"):
        paths[key].mkdir(parents=True, exist_ok=True)
    write_abacus_stru(paths["abacus_pw"] / "STRU")
    write_kpt(paths["abacus_pw"] / "KPT")
    write_abacus_input(paths["abacus_pw"] / "INPUT", solver)
    write_abacus_sbatch(paths["abacus_pw"] / "run_abacus.sbatch", solver)
    write_lobsterin(paths["direct"] / "lobsterin", direct=True)
    write_lobster_sbatch(paths["direct"] / "run_lobster.sbatch", "dia-lob-dir")
    write_lobsterin(paths["generic"] / "lobsterin", direct=False)
    write_lobster_sbatch(paths["generic"] / "run_lobster.sbatch", "dia-lob-gen")
    write_json(paths["jobs"], {"jobs": {}})
    write_json(
        paths["state"],
        {
            "root": str(root),
            "abacus_pw": str(paths["abacus_pw"]),
            "direct_lobster": str(paths["direct"]),
            "generic_lobster": str(paths["generic"]),
            "abacus_solver": solver,
            "kmesh": list(K_MESH),
            "nbands": NBANDS,
            "ecut_ry": ECUT_RY,
            "apns_c_upf": str(APNS_C_UPF),
            "diagnostic_paw_upf": str(LOBSTER_QE_C_PAW),
        },
    )
    print(f"prepared {root}")


def submit(root: Path, step: str) -> int:
    paths = layout(root)
    run_dirs = {"abacus-pw": paths["abacus_pw"], "lobster-direct": paths["direct"], "lobster-generic": paths["generic"]}
    scripts = {"abacus-pw": "run_abacus.sbatch", "lobster-direct": "run_lobster.sbatch", "lobster-generic": "run_lobster.sbatch"}
    if step not in run_dirs:
        raise ValueError(f"unknown step: {step}")
    result = run(["sbatch", scripts[step]], cwd=run_dirs[step])
    match = re.search(r"Submitted batch job\s+(\d+)", result.stdout)
    if not match:
        raise RuntimeError(result.stdout + result.stderr)
    job_id = int(match.group(1))
    jobs = read_json(paths["jobs"], {"jobs": {}})
    jobs["jobs"][step] = {"job_id": job_id, "run_dir": str(run_dirs[step])}
    write_json(paths["jobs"], jobs)
    print(f"{step}: submitted job {job_id}")
    return job_id


def slurm_state(job_id: int) -> str:
    sq = run(["squeue", "-j", str(job_id), "-h", "-o", "%T"], check=False)
    if sq.stdout.strip():
        return sq.stdout.strip()
    sa = run(["sacct", "-X", "-j", str(job_id), "--format=JobIDRaw,State,ExitCode", "-n", "-P"], check=False)
    for line in sa.stdout.splitlines():
        cols = line.split("|")
        if len(cols) >= 3 and cols[0] == str(job_id):
            return f"{cols[1]} ({cols[2]})"
    return "UNKNOWN"


def status(root: Path) -> None:
    paths = layout(root)
    jobs = read_json(paths["jobs"], {"jobs": {}})["jobs"]
    markers = {
        "abacus-pw": [paths["abacus_pw"] / "OUT.ABACUS" / "running_scf.log", paths["abacus_pw"] / "OUT.ABACUS" / "WAVEFUNC1.dat"],
        "lobster-direct": [paths["direct"] / "lobsterout"],
        "lobster-generic": [paths["generic"] / "lobsterout"],
    }
    for step, files in markers.items():
        ready = all(path.exists() for path in files)
        job = jobs.get(step)
        job_text = "not submitted" if not job else f"job {job['job_id']}: {slurm_state(int(job['job_id']))}"
        print(f"[{step}] {job_text}; outputs_ready={ready}")


def parse_stru(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    lat0 = float(lines[lines.index("LATTICE_CONSTANT") + 1])
    idx = lines.index("LATTICE_VECTORS")
    cell = np.array([[float(x) for x in lines[idx + 1 + i].split()[:3]] for i in range(3)]) * lat0 * BOHR_TO_ANG
    idx = lines.index("ATOMIC_POSITIONS")
    mode = lines[idx + 1].lower()
    natom = int(lines[idx + 4])
    coords = np.array([[float(x) for x in lines[idx + 5 + i].split()[:3]] for i in range(natom)])
    frac = coords if mode.startswith("direct") else coords @ np.linalg.inv(cell)
    return cell, frac


def read_record_header(fh: Any) -> int:
    raw = fh.read(4)
    if len(raw) != 4:
        raise EOFError
    return struct.unpack("<i", raw)[0]


def read_wfc_pw(path: Path) -> WfcPw:
    with path.open("rb") as fh:
        n = read_record_header(fh)
        if n != 72:
            raise ValueError(f"{path}: unexpected header record length {n}")
        payload = fh.read(72)
        end = read_record_header(fh)
        if end != 72:
            raise ValueError(f"{path}: bad header trailer {end}")
        ik, nk = struct.unpack("<ii", payload[:8])
        kx, ky, kz, weight = struct.unpack("<dddd", payload[8:40])
        ng, nbands = struct.unpack("<ii", payload[40:48])
        ecut_ry, lat0, tpiba = struct.unpack("<ddd", payload[48:72])

        n = read_record_header(fh)
        recip = np.frombuffer(fh.read(n), dtype="<f8").reshape(3, 3).copy()
        end = read_record_header(fh)
        if n != 72 or end != 72:
            raise ValueError(f"{path}: bad reciprocal lattice record")

        n = read_record_header(fh)
        if n != ng * 12:
            raise ValueError(f"{path}: bad G-vector record length {n}; expected {ng * 12}")
        g_indices = np.frombuffer(fh.read(n), dtype="<i4").reshape(ng, 3).copy()
        end = read_record_header(fh)
        if end != n:
            raise ValueError(f"{path}: bad G-vector trailer")

        coeffs = np.zeros((nbands, ng), dtype=np.complex128)
        for ib in range(nbands):
            n = read_record_header(fh)
            if n != ng * 16:
                raise ValueError(f"{path}: bad band {ib + 1} record length {n}; expected {ng * 16}")
            raw = np.frombuffer(fh.read(n), dtype="<f8").reshape(ng, 2)
            coeffs[ib] = raw[:, 0] + 1j * raw[:, 1]
            end = read_record_header(fh)
            if end != n:
                raise ValueError(f"{path}: bad band {ib + 1} trailer")
    return WfcPw(ik, nk, np.array([kx, ky, kz]), weight, ng, nbands, ecut_ry, lat0, tpiba, recip, g_indices, coeffs)


def parse_istate(path: Path) -> dict[int, dict[str, list[float]]]:
    data: dict[int, dict[str, list[float]]] = {}
    current = None
    for line in path.read_text().splitlines():
        match = re.search(r"Kpoint\s*=\s*(\d+)", line)
        if match:
            current = int(match.group(1))
            coord_matches = re.findall(r"\(([^)]*)\)", line)
            coords = []
            if coord_matches:
                try:
                    coords = [float(x) for x in coord_matches[-1].split()[:3]]
                except ValueError:
                    coords = []
            data[current] = {"eigen_values": [], "occupations": [], "coordinates": coords}
            continue
        cols = line.split()
        if current and len(cols) >= 3 and cols[0].isdigit():
            data[current]["eigen_values"].append(float(cols[1]))
            data[current]["occupations"].append(float(cols[2]))
    return data


def parse_fermi(log_path: Path) -> float:
    text = log_path.read_text(errors="ignore")
    matches = re.findall(r"EFERMI\s*=\s*([-+0-9.Ee]+)|Fermi energy is\s*([-+0-9.Ee]+)", text)
    vals = [float(a or b) for a, b in matches if a or b]
    return vals[-1] if vals else 0.0


def parse_cube_grid(path: Path) -> list[int]:
    lines = path.read_text(errors="ignore").splitlines()
    if len(lines) < 6:
        return [0, 0, 0]
    return [abs(int(lines[i].split()[0])) for i in range(3, 6)]


def parse_paw_upf(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="ignore")
    header_match = re.search(r"<PP_HEADER\s+([^>]*)>", text, flags=re.S)
    if not header_match:
        raise ValueError(f"{path} lacks PP_HEADER")

    def attrs(attr_text: str) -> dict[str, str]:
        return {key: val for key, val in re.findall(r'([A-Za-z0-9_]+)="([^"]*)"', attr_text)}

    def tag_vector(tag: str, body: str = text) -> list[float]:
        match = re.search(rf"<{re.escape(tag)}\b[^>]*>(.*?)</{re.escape(tag)}>", body, flags=re.S)
        if not match:
            raise ValueError(f"{path} lacks {tag}")
        return [float(x) for x in match.group(1).split()]

    def numbered(prefix: str, body: str = text) -> list[tuple[int, dict[str, str], list[float]]]:
        out = []
        for match in re.finditer(rf"<{prefix}\.(\d+)\b([^>]*)>(.*?)</{prefix}\.\1>", body, flags=re.S):
            out.append((int(match.group(1)), attrs(match.group(2)), [float(x) for x in match.group(3).split()]))
        return sorted(out, key=lambda item: item[0])

    header = attrs(header_match.group(1))
    beta_elems = numbered("PP_BETA")
    ae_elems = numbered("PP_AEWFC")
    ps_elems = numbered("PP_PSWFC")
    if not beta_elems or not ae_elems or not ps_elems:
        raise ValueError(f"{path} does not contain complete PAW beta/AE/PS arrays")
    q_match = re.search(r"<PP_AUGMENTATION\b.*?<PP_Q\b[^>]*>(.*?)</PP_Q>", text, flags=re.S)
    if q_match is None:
        raise ValueError(f"{path} lacks PP_Q augmentation charge matrix")
    nproj = len(beta_elems)
    labels = [meta.get("label", "").strip().lower() for _, meta, _ in beta_elems]
    valence = []
    for label in labels:
        if label and label not in valence:
            valence.append(label)
    return {
        "element": header["element"].strip(),
        "title": header.get("generated", "PAW UPF"),
        "n_electrons": float(header["z_valence"]),
        "n_projections": nproj,
        "angular_momentum": [int(meta["angular_momentum"]) for _, meta, _ in beta_elems],
        "valence": valence,
        "augmentation_charges": np.array([float(x) for x in q_match.group(1).split()], dtype=float)
        .reshape(nproj, nproj)
        .tolist(),
        "grid": tag_vector("PP_R"),
        "drgrid": tag_vector("PP_RAB"),
        "projectors": [vec for _, _, vec in beta_elems],
        "all_electron_wfc": [vec for _, _, vec in ae_elems],
        "ps_partialwave": [vec for _, _, vec in ps_elems],
    }


def g_to_miller(g_indices: np.ndarray, grid: list[int]) -> np.ndarray:
    grid_arr = np.array(grid, dtype=np.int32)
    if np.any(grid_arr <= 0):
        grid_arr = g_indices.max(axis=0).astype(np.int32) + 1
    out = g_indices.astype(np.int32).copy()
    half = grid_arr // 2
    for i in range(3):
        mask = out[:, i] > half[i]
        out[mask, i] -= grid_arr[i]
    return out


def build_generic(root: Path) -> None:
    paths = layout(root)
    out = paths["abacus_pw"] / "OUT.ABACUS"
    wfc_files = sorted(out.glob("WAVEFUNC*.dat"), key=lambda p: int(re.search(r"(\d+)", p.name).group(1)))
    if not wfc_files:
        raise FileNotFoundError(f"no WAVEFUNC*.dat in {out}")
    istate = parse_istate(out / "istate.info")
    grid = parse_cube_grid(out / "SPIN1_CHG.cube")
    fermi = parse_fermi(out / "running_scf.log")
    cell, frac = parse_stru(paths["abacus_pw"] / "STRU")
    paw = parse_paw_upf(LOBSTER_QE_C_PAW)

    paths["generic"].mkdir(parents=True, exist_ok=True)
    kdir = paths["generic"] / "LOBSTER_Kpoints"
    kdir.mkdir(exist_ok=True)
    k_points = []
    n_bands = None
    ecut_ev = None
    for wfc_path in wfc_files:
        wfc = read_wfc_pw(wfc_path)
        n_bands = wfc.nbands
        ecut_ev = wfc.ecut_ry * RY_TO_EV
        with h5py.File(kdir / f"kPoint{wfc.ik}.hdf5", "w") as h5:
            h5.create_dataset("Miller", data=g_to_miller(wfc.g_indices, grid).astype(np.int32))
            h5.create_dataset("PWCoeffs", data=wfc.coeffs.astype(np.complex128))
        state = istate.get(wfc.ik, {"eigen_values": [0.0] * wfc.nbands, "occupations": [0.0] * wfc.nbands})
        occupations = state["occupations"][: wfc.nbands]
        if wfc.weight > 0:
            # ABACUS istate.info reports occupations multiplied by the k-point weight.
            occupations = [min(1.0, max(0.0, x / wfc.weight)) for x in occupations]
        k_points.append(
            {
                "coordinates": state.get("coordinates") or wfc.k_cart.tolist(),
                "weight": wfc.weight,
                "occupations": occupations,
                "eigen_values": state["eigen_values"][: wfc.nbands],
                "n_plane_waves": wfc.ng,
            }
        )

    data = {
        "file": {"used_program": "ABACUS", "no_sym": 1, "no_inv": 1},
        "cell": {
            "lattice_vectors_real": cell.tolist(),
            "cell_volume": float(abs(np.linalg.det(cell))),
            "atomic_structure": [{"element": "C", "coordinates": x.tolist()} for x in frac],
        },
        "wave_function": {
            "energy_minimum": -22.0,
            "energy_maximum": 18.0,
            "energy_step": 1000,
            "size_of_g_vector_grid": [float(x) for x in grid],
            "cutoff_energy": float(ecut_ev),
            "fermi_energy": float(fermi),
            "n_bands": int(n_bands or NBANDS),
            "n_spins": 1,
            "spin_channels": [{"k_points": k_points}],
        },
        "paw": [paw],
    }
    write_json(paths["generic"] / "LobsterInput.json", data)
    write_json(
        paths["analysis"] / "generic_manifest.json",
        {
            "status": "generated",
            "n_kpoints": len(k_points),
            "fft_grid": grid,
            "warning": "PAW data were taken from the LOBSTER QE diamond example, not from the APNS NC pseudopotential used by ABACUS. This is a format diagnostic only.",
        },
    )
    print(f"generated Generic input for {len(k_points)} k-points")


def prepare_direct(root: Path) -> None:
    paths = layout(root)
    out = paths["abacus_pw"] / "OUT.ABACUS"
    if not out.exists():
        raise FileNotFoundError(out)
    paths["direct"].mkdir(parents=True, exist_ok=True)
    for name in ["INPUT", "STRU", "KPT"]:
        shutil.copy2(paths["abacus_pw"] / name, paths["direct"] / name)
    for pattern in ["running_scf.log", "istate.info", "kpoints", "WAVEFUNC*.dat", "SPIN1_CHG.cube"]:
        for src in out.glob(pattern):
            if src.is_file():
                shutil.copy2(src, paths["direct"] / src.name)
    print(f"prepared direct LOBSTER probe in {paths['direct']}")


def summarize_outputs(root: Path) -> dict[str, Any]:
    paths = layout(root)
    out = paths["abacus_pw"] / "OUT.ABACUS"
    wfc_files = sorted(out.glob("WAVEFUNC*.dat"))
    direct_out = paths["direct"] / "lobsterout"
    generic_out = paths["generic"] / "lobsterout"
    direct_text = direct_out.read_text(errors="ignore") if direct_out.exists() else ""
    generic_text = generic_out.read_text(errors="ignore") if generic_out.exists() else ""
    direct_err = (paths["direct"] / "lobster.err").read_text(errors="ignore") if (paths["direct"] / "lobster.err").exists() else ""
    generic_err = (paths["generic"] / "lobster.err").read_text(errors="ignore") if (paths["generic"] / "lobster.err").exists() else ""
    return {
        "abacus_converged": "charge density convergence is achieved" in (out / "running_scf.log").read_text(errors="ignore") if (out / "running_scf.log").exists() else False,
        "wfc_file_count": len(wfc_files),
        "has_istate": (out / "istate.info").exists(),
        "has_charge_cube": (out / "SPIN1_CHG.cube").exists(),
        "direct_lobsterout": direct_out.exists(),
        "generic_lobsterout": generic_out.exists(),
        "direct_error_excerpt": (direct_err or direct_text)[-2000:],
        "generic_error_excerpt": (generic_err or generic_text)[-2000:],
        "generic_reference_lobsterout": (root / "lobster" / "generic_reference" / "lobsterout").exists(),
        "generic_reference_excerpt": ((root / "lobster" / "generic_reference" / "lobsterout").read_text(errors="ignore")[-2000:] if (root / "lobster" / "generic_reference" / "lobsterout").exists() else ""),
    }


def write_report(root: Path) -> None:
    paths = layout(root)
    summary = summarize_outputs(root)
    manifest = read_json(paths["analysis"] / "generic_manifest.json", {})
    report = f"""# ABACUS PW 输出用于 LOBSTER COHP 的可行性测试

## 测试对象

- 体系：金刚石 primitive cell，2 个 C 原子。
- ABACUS：LTS v3.10.1，PW 基组，`ecutwfc {ECUT_RY}` Ry，`nbands {NBANDS}`，`KPT {K_MESH[0]} {K_MESH[1]} {K_MESH[2]}`，`symmetry -1`。
- ABACUS 关键输出参数：`out_wfc_pw 2`、`out_chg 1`、`out_band 1`。
- ABACUS 赝势：APNS `C.upf`，该文件为 norm-conserving UPF，`is_paw=\"F\"`。
- LOBSTER：`{LOBSTER_BIN}`。

## 依据文件

- ABACUS 波函数输出文档：`abacus-develop/docs/advanced/elec_properties/wfc.md`。
- ABACUS INPUT 参数文档：`abacus-develop/docs/advanced/input_files/input-main.md`。
- ABACUS PW 波函数二进制写出源码：`abacus-develop/source/module_io/write_wfc_pw.cpp`。
- LOBSTER Generic 格式说明：`tools/lobster-5.1.1/Generic/README.md`。
- LOBSTER 用户手册：`tools/lobster-5.1.1/Lobster_Users_Guide_5.1.1.pdf`。
- LOBSTER 方法论文页面：<https://www.oqi.ox.ac.uk/publication/1114472/europe-pubmed-central>，其摘要说明 LOBSTER 的投影方法基于 PAW DFT 计算。
- Quantum ESPRESSO 官方赝势页：<https://www.quantum-espresso.org/pseudopotentials/>，确认 QE 本身支持 NC、US 和 PAW 多类赝势。
- Quantum ESPRESSO UPF 说明：<https://pseudopotentials.quantum-espresso.org/home/unified-pseudopotential-format>，确认 UPF 容器可存储 NC、US、PAW 等不同类型赝势。

## 已完成输出检查

- ABACUS SCF 收敛：`{summary['abacus_converged']}`。
- `WAVEFUNC*.dat` 数量：`{summary['wfc_file_count']}`。
- `istate.info`：`{summary['has_istate']}`。
- `SPIN1_CHG.cube`：`{summary['has_charge_cube']}`。
- Generic manifest：`{json.dumps(manifest, ensure_ascii=False)}`。

## LOBSTER 对 PW 赝势类型的限制

需要区分两个问题：

1. QE 作为 PW-DFT 程序支持多类赝势。QE 官方赝势页明确列出 PAW、USPP 和 NC PP；UPF 官方说明也表明 UPF 是一个可容纳 NC、US、PAW 等类型的统一格式。因此，“LOBSTER 支持 QE 输出”并不等价于“LOBSTER 支持任意 QE 赝势类型输出”。
2. LOBSTER 对其可读取的 PW 结果有更强限制。LOBSTER 用户手册在 VASP 准备部分明确要求使用 PAW potential，而不是 ultrasoft pseudopotential；ABINIT 部分说明 VASP 的要求同样适用，并进一步限定为 PAW-XML 数据集；QE 部分说明 VASP 和 ABINIT 部分的要求同样适用于 QE，并在 FAQ 中要求读取本次计算使用的 PAW 数据 UPF 文件。

因此，对 LOBSTER 的标准 VASP/ABINIT/QE 接口而言，PW 计算结果需要来自 PAW 赝势/PAW 数据集。QE 虽然能用 NC、US、PAW 计算，但 LOBSTER 的 QE 接口要求的是 QE+PAW 数据，而不是 QE+NC 或 QE+US 的一般 PW 输出。

这一点也和 LOBSTER 的方法论文一致：其局域轨道投影方法建立在从 PAW DFT 计算重构化学信息的框架上。对本项目而言，APNS `C.upf` 的 `PP_HEADER` 为 `pseudo_type="NC"` 且 `is_paw="F"`，不含 LOBSTER Generic/PAW 重构所需的 augmentation charges、projectors、all-electron partial waves 等 PAW 数据。因此，当前 ABACUS PW+APNS NC 输出即使具备 PW coefficients，也缺少 LOBSTER 所需的 PAW 赝势侧信息。

## 直接读取 ABACUS PW 输出

LOBSTER 原生示例和二进制字符串显示其可自动识别 VASP、Quantum ESPRESSO、ABINIT 和 Generic 输入，但没有 ABACUS 原生检测路径。本测试将 ABACUS 的 `INPUT/STRU/KPT/WAVEFUNC*.dat/istate.info` 放入 LOBSTER 目录后运行。

结果文件存在：`{summary['direct_lobsterout']}`。

错误摘要：

```text
{summary['direct_error_excerpt'].strip()}
```

结论：ABACUS PW 输出不能像 VASP 的 `WAVECAR/vasprun.xml` 那样被 LOBSTER 直接消费。

## Generic 转换测试

转换脚本从 ABACUS `WAVEFUNC*.dat` 解析 k 点、PW 数量、G 网格索引、band 数、ecut 和 complex128 波函数系数，并写入 LOBSTER Generic 的 `LOBSTER_Kpoints/kPoint*.hdf5`：

- `Miller`: `(nPW, 3) int32`
- `PWCoeffs`: `(nBands, nPW) complex128`

能级和占据来自 `istate.info`，FFT grid 来自 `SPIN1_CHG.cube`，cell/atomic structure 来自 `STRU`。

重要限制：APNS `C.upf` 是 NC 赝势，不包含 LOBSTER Generic 文档要求的 PAW augmentation、projectors、all-electron wavefunctions。当前 Generic 诊断包为了测试格式通路，使用了 LOBSTER 自带 QE diamond 示例中的 C PAW UPF 填充 `paw` 字段。这与 ABACUS SCF 的赝势不一致，因此即使 LOBSTER 成功输出 COHP，也只能证明接口格式可能打通，不能作为严格科学结果。

Generic LOBSTER 结果文件存在：`{summary['generic_lobsterout']}`。

错误摘要：

```text
{summary['generic_error_excerpt'].strip()}
```

额外校验：将 LOBSTER 发行包自带的 Generic 示例原样复制到本测试目录并运行，本机 `lobster-5.1.1` 同样停在程序识别阶段，未打开 `LobsterInput.json`。这说明当前二进制/运行方式下 Generic detector 没有被触发；该现象独立于 ABACUS 转换数据本身。

发行包 Generic 示例摘要：

```text
{summary['generic_reference_excerpt'].strip()}
```

## 结论

1. **直接答案：不能直接做。** LOBSTER 不原生识别 ABACUS PW 输出。
2. **工程路径：理论上可以通过 LOBSTER Generic 做。** ABACUS `out_wfc_pw 2` 给出了构造 `PWCoeffs`/`Miller` 的核心数据，`out_band 1` 给出 `istate.info`，`out_chg 1` 可辅助获得 FFT grid。
3. **科学阻塞：PAW 数据一致性。** 严格 LOBSTER COHP 需要与 PW 波函数同源的 PAW projectors、augmentation charges、AE/PS partial waves。APNS 当前用于 ABACUS PW 的 C.upf 是 NC，不满足这一要求。
4. **推荐实践：** 当前生产级 ABACUS COHP 仍应使用本项目的 LCAO-COHP 后处理；ABACUS PW 到 LOBSTER 需要新增稳定的 Generic 导出器，并配套 ABACUS 可用且与 LOBSTER 数据结构一致的 PAW 势库。
"""
    REPORT.write_text(report)
    print(f"wrote {REPORT}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT_DEFAULT)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_prepare = sub.add_parser("prepare")
    p_prepare.add_argument("--solver", choices=["bpcg", "dav_subspace"], default="bpcg")
    p_submit = sub.add_parser("submit")
    p_submit.add_argument("step", choices=["abacus-pw", "lobster-direct", "lobster-generic"])
    sub.add_parser("status")
    sub.add_parser("prepare-direct")
    sub.add_parser("build-generic")
    sub.add_parser("report")
    args = parser.parse_args()
    if args.cmd == "prepare":
        prepare(args.root, args.solver)
    elif args.cmd == "submit":
        submit(args.root, args.step)
    elif args.cmd == "status":
        status(args.root)
    elif args.cmd == "prepare-direct":
        prepare_direct(args.root)
    elif args.cmd == "build-generic":
        build_generic(args.root)
    elif args.cmd == "report":
        write_report(args.root)


if __name__ == "__main__":
    main()
