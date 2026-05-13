#!/usr/bin/env python3
"""End-to-end helpers for Pt(111)-CO ABACUS LCAO-COHP runs."""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "refs"))

import cohp  # noqa: E402


PP_ORB = Path.home() / "PP_ORB"
PP_DIR = PP_ORB / "PP"
COMPACT_ORB_DIR = PP_ORB / "ORB"
PRECISION_ORB_DIR = PP_ORB / "apns-orbitals-precision-v1"

SPECIES = {
    "Pt": {
        "mass": 195.084,
        "pseudo": "Pt_ONCV_PBE-1.2.upf",
        "compact_orb": "Pt_gga_8au_100Ry_4s2p2d1f.orb",
        "precision_orb": "Pt_gga_10au_100Ry_6s3p3d2f.orb",
        "compact_shells": {"s": 4, "p": 2, "d": 2, "f": 1},
        "precision_shells": {"s": 6, "p": 3, "d": 3, "f": 2},
    },
    "C": {
        "mass": 12.011,
        "pseudo": "C.upf",
        "compact_orb": "C_gga_8au_100Ry_2s2p1d.orb",
        "precision_orb": "C_gga_10au_100Ry_3s3p2d.orb",
        "compact_shells": {"s": 2, "p": 2, "d": 1},
        "precision_shells": {"s": 3, "p": 3, "d": 2},
    },
    "O": {
        "mass": 15.999,
        "pseudo": "O.upf",
        "compact_orb": "O_gga_6au_100Ry_2s2p1d.orb",
        "precision_orb": "O_gga_10au_100Ry_3s3p2d1f.orb",
        "compact_shells": {"s": 2, "p": 2, "d": 1},
        "precision_shells": {"s": 3, "p": 3, "d": 2, "f": 1},
    },
}

L_MULT = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}


@dataclass
class AbacusAtom:
    symbol: str
    position: np.ndarray
    movable: tuple[int, int, int]
    source_index: int


def shell_count(shells: dict[str, int]) -> int:
    return sum(nzeta * L_MULT[label] for label, nzeta in shells.items())


def shell_ranges(start: int, shells: dict[str, int]) -> dict[str, list[int]]:
    ranges = {}
    cursor = start
    for label in ["s", "p", "d", "f", "g"]:
        if label not in shells:
            continue
        width = shells[label] * L_MULT[label]
        ranges[label] = list(range(cursor, cursor + width))
        cursor += width
    ranges["all"] = list(range(start, cursor))
    return ranges


def build_atoms() -> tuple[np.ndarray, list[AbacusAtom], dict]:
    from ase.build import add_adsorbate, fcc111

    slab = fcc111("Pt", size=(2, 2, 4), a=3.92, vacuum=15.0, orthogonal=True)
    pt_positions = slab.get_positions()
    top_z = np.max(pt_positions[:, 2])
    top_candidates = np.where(np.isclose(pt_positions[:, 2], top_z, atol=1e-3))[0]
    center_xy = np.mean(slab.cell.array[:2, :2], axis=0)
    top_pt = min(
        top_candidates,
        key=lambda idx: np.linalg.norm(pt_positions[idx, :2] - center_xy),
    )
    add_adsorbate(slab, "C", height=1.85, position=pt_positions[top_pt, :2])
    add_adsorbate(slab, "O", height=1.85 + 1.15, position=pt_positions[top_pt, :2])

    old_cell = slab.cell.array
    old_pos = slab.get_positions()
    new_cell = np.array([old_cell[0], old_cell[2], old_cell[1]], dtype=float)[:, [0, 2, 1]]
    new_pos = old_pos[:, [0, 2, 1]]

    pt_z = old_pos[:16, 2]
    layers = sorted({round(z, 6) for z in pt_z})
    fixed_z = set(layers[:2])

    atoms = []
    for idx, (symbol, pos) in enumerate(zip(slab.get_chemical_symbols(), new_pos)):
        if symbol == "Pt" and round(old_pos[idx, 2], 6) in fixed_z:
            movable = (0, 0, 0)
        else:
            movable = (1, 1, 1)
        atoms.append(AbacusAtom(symbol=symbol, position=pos, movable=movable, source_index=idx))

    atoms = [a for a in atoms if a.symbol == "Pt"] + [a for a in atoms if a.symbol == "C"] + [a for a in atoms if a.symbol == "O"]
    source_to_abacus = {atom.source_index: i for i, atom in enumerate(atoms)}
    mapping = {
        "top_pt_atom_index": source_to_abacus[int(top_pt)],
        "carbon_atom_index": source_to_abacus[16],
        "oxygen_atom_index": source_to_abacus[17],
        "vacuum_axis": "B",
        "fixed_pt_layers": 2,
        "initial_distances_angstrom": {"Pt-C": 1.85, "C-O": 1.15},
    }
    return new_cell, atoms, mapping


def orbital_mapping(atoms: list[AbacusAtom], basis: str) -> dict:
    key = f"{basis}_shells"
    cursor = 0
    atoms_out = []
    for i, atom in enumerate(atoms):
        shells = SPECIES[atom.symbol][key]
        ranges = shell_ranges(cursor, shells)
        atom_entry = {
            "atom_index": i,
            "symbol": atom.symbol,
            "source_index": atom.source_index,
            "orbital_start": cursor,
            "orbital_stop": cursor + len(ranges["all"]),
            "shells": ranges,
        }
        atoms_out.append(atom_entry)
        cursor = atom_entry["orbital_stop"]
    return {
        "basis": basis,
        "total_orbitals": cursor,
        "atoms": atoms_out,
        "orbitals_per_element": {
            element: shell_count(SPECIES[element][key]) for element in SPECIES
        },
    }


def write_stru(path: Path, cell: np.ndarray, atoms: list[AbacusAtom], basis: str) -> None:
    orb_key = "compact_orb" if basis == "compact" else "precision_orb"
    grouped = {symbol: [atom for atom in atoms if atom.symbol == symbol] for symbol in ["Pt", "C", "O"]}
    lines = [
        "ATOMIC_SPECIES",
        *[f"{symbol} {SPECIES[symbol]['mass']:.6f} {SPECIES[symbol]['pseudo']}" for symbol in ["Pt", "C", "O"]],
        "",
        "NUMERICAL_ORBITAL",
        *[SPECIES[symbol][orb_key] for symbol in ["Pt", "C", "O"]],
        "",
        "LATTICE_CONSTANT",
        "1.889726125457828",
        "",
        "LATTICE_VECTORS",
        *[" ".join(f"{x:.12f}" for x in row) for row in cell],
        "",
        "ATOMIC_POSITIONS",
        "Cartesian_angstrom",
    ]
    for symbol in ["Pt", "C", "O"]:
        lines += [symbol, "0.0", str(len(grouped[symbol]))]
        for atom in grouped[symbol]:
            x, y, z = atom.position
            mx, my, mz = atom.movable
            lines.append(f"{x:.12f} {y:.12f} {z:.12f} {mx} {my} {mz}")
    path.write_text("\n".join(lines) + "\n")


def write_input(path: Path, nspin: int, calculation: str, basis_dir: Path, out_hs: bool) -> None:
    lines = [
        "INPUT_PARAMETERS",
        f"calculation {calculation}",
        "basis_type lcao",
        f"pseudo_dir {PP_DIR}",
        f"orbital_dir {basis_dir}",
        "ecutwfc 100",
        "device gpu",
        "ks_solver cusolver",
        f"nspin {nspin}",
        "smearing_method gaussian",
        "smearing_sigma 0.004",
        "mixing_beta 0.4",
        "kspacing 0.14",
        "scf_thr 1e-7",
        "scf_nmax 120",
    ]
    if calculation == "relax":
        lines += [
            "force_thr_ev 0.03",
            "relax_method bfgs_trad",
            "relax_nmax 80",
            "out_stru 1",
        ]
    if out_hs:
        lines += ["out_mat_hs 1 8", "out_wfc_lcao 1", "out_app_flag 1"]
    path.write_text("\n".join(lines) + "\n")


def write_sbatch(path: Path, job_name: str) -> None:
    path.write_text(
        "\n".join(
            [
                "#!/bin/bash",
                f"#SBATCH --job-name={job_name}",
                "#SBATCH --partition=4V100",
                "#SBATCH --nodes=1",
                "#SBATCH --ntasks=1",
                "#SBATCH --gpus-per-node=1",
                "#SBATCH --qos=rush-1o2gpu",
                "#SBATCH --output=slurm-%j.out",
                "#SBATCH --error=slurm-%j.err",
                "",
                "set -euo pipefail",
                "module load abacus/LTSv3.10.1-sm70-auto",
                "export OMP_NUM_THREADS=8",
                "export OMP_PROC_BIND=spread",
                "export OMP_PLACES=cores",
                "nvidia-smi dmon -s pucvmte -o T > nvdmon_job-${SLURM_JOB_ID}.log &",
                "abacus > abacus.stdout 2> abacus.stderr",
                "",
            ]
        )
    )


def prepare(root: Path) -> None:
    cell, atoms, base_mapping = build_atoms()
    root.mkdir(parents=True, exist_ok=True)
    structure_dir = root / "structure"
    structure_dir.mkdir(exist_ok=True)
    write_stru(structure_dir / "STRU_initial_compact", cell, atoms, "compact")
    (structure_dir / "mapping_initial.json").write_text(
        json.dumps(
            {
                **base_mapping,
                "compact_orbitals": orbital_mapping(atoms, "compact"),
                "precision_orbitals": orbital_mapping(atoms, "precision"),
            },
            indent=2,
        )
    )
    for nspin in [1, 2]:
        relax_dir = root / f"nspin{nspin}" / "relax"
        relax_dir.mkdir(parents=True, exist_ok=True)
        write_stru(relax_dir / "STRU", cell, atoms, "compact")
        write_input(relax_dir / "INPUT", nspin, "relax", COMPACT_ORB_DIR, out_hs=False)
        write_sbatch(relax_dir / "run_abacus.sbatch", f"ptco-r{nspin}")
        shutil.copy2(structure_dir / "mapping_initial.json", relax_dir / "mapping.json")


def prepare_final(root: Path, nspin: int) -> None:
    relax_dir = root / f"nspin{nspin}" / "relax"
    final_dir = root / f"nspin{nspin}" / "final_scf"
    source_stru = relax_dir / "OUT.ABACUS" / "STRU_ION_D"
    if not source_stru.exists():
        raise FileNotFoundError(f"Relaxed structure is missing: {source_stru}")
    final_dir.mkdir(parents=True, exist_ok=True)
    text = source_stru.read_text()
    for symbol in ["Pt", "C", "O"]:
        text = text.replace(SPECIES[symbol]["compact_orb"], SPECIES[symbol]["precision_orb"])
    (final_dir / "STRU").write_text(text)
    write_input(final_dir / "INPUT", nspin, "scf", PRECISION_ORB_DIR, out_hs=True)
    write_sbatch(final_dir / "run_abacus.sbatch", f"ptco-s{nspin}")
    shutil.copy2(relax_dir / "mapping.json", final_dir / "mapping.json")


def parse_stru_positions(path: Path) -> tuple[list[str], np.ndarray]:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    lat0 = 1.0
    if "LATTICE_CONSTANT" in lines:
        lat0 = float(lines[lines.index("LATTICE_CONSTANT") + 1].split()[0])
    if "LATTICE_VECTORS" in lines:
        cell = np.array(
            [[float(x) for x in lines[lines.index("LATTICE_VECTORS") + 1 + i].split()[:3]] for i in range(3)],
            dtype=float,
        )
        cell_angstrom = cell * lat0 / 1.889726125457828
    else:
        cell_angstrom = None
    idx = lines.index("ATOMIC_POSITIONS")
    coord_type = lines[idx + 1].split()[0]
    if coord_type not in {"Cartesian_angstrom", "Direct"}:
        raise ValueError(f"Only Cartesian_angstrom or Direct STRU is supported, got {lines[idx + 1]}")
    symbols = []
    positions = []
    cursor = idx + 2
    while cursor < len(lines):
        symbol = lines[cursor].split()[0]
        mag = lines[cursor + 1]
        natom = int(lines[cursor + 2].split()[0])
        del mag
        cursor += 3
        for _ in range(natom):
            parts = lines[cursor].split()
            symbols.append(symbol)
            position = np.array([float(parts[0]), float(parts[1]), float(parts[2])], dtype=float)
            if coord_type == "Direct":
                if cell_angstrom is None:
                    raise ValueError(f"Direct STRU requires LATTICE_VECTORS: {path}")
                position = position @ cell_angstrom
            positions.append(position.tolist())
            cursor += 1
    return symbols, np.array(positions, dtype=float)


def read_fermi(out_dir: Path) -> float:
    values = cohp.rao.read_etraj_fromlog(str(out_dir / "running_scf.log"), term="fermi")
    return float(values[-1])


def write_cohp_for_pair(out_dir: Path, cohp_dir: Path, label: str, pt_orbs: list[int], c_orbs: list[int], nspin: int) -> dict:
    cohp_dir.mkdir(parents=True, exist_ok=True)
    spins = ["sum"] if nspin == 1 else ["up", "down", "sum"]
    metrics = {}
    for spin in spins:
        prefix = cohp_dir / f"{label}_{spin}"
        e, y = cohp.run_outdir(
            out_dir=out_dir,
            atomI_orbs=pt_orbs,
            atomJ_orbs=c_orbs,
            testmethod="COHP",
            de=0.05,
            smooth=True,
            smooth_nstddev=4,
            shift_toefermi=True,
            invert_COHP=True,
            emin=-12,
            emax=8,
            width=None,
            output_prefix=str(prefix),
            spin=spin,
        )
        ef = read_fermi(out_dir)
        occupied = e <= ef
        icohp = float(np.trapezoid(y[occupied], e[occupied])) if np.any(occupied) else float("nan")
        metrics[spin] = {
            "icohp": icohp,
            "minus_icohp": -icohp,
            "min_cohp": float(np.min(y)),
            "max_cohp": float(np.max(y)),
            "points": int(len(e)),
        }
    return metrics


def analyze(root: Path, nspin: int) -> None:
    final_dir = root / f"nspin{nspin}" / "final_scf"
    out_dir = final_dir / "OUT.ABACUS"
    cohp_dir = root / f"nspin{nspin}" / "cohp"
    mapping = json.loads((final_dir / "mapping.json").read_text())
    precision_atoms = mapping["precision_orbitals"]["atoms"]
    top_pt = precision_atoms[mapping["top_pt_atom_index"]]
    carbon = precision_atoms[mapping["carbon_atom_index"]]
    pt_all = top_pt["shells"]["all"]
    c_all = carbon["shells"]["all"]
    metrics = {
        "top_Pt_C_total": write_cohp_for_pair(out_dir, cohp_dir, "top_Pt_C_total", pt_all, c_all, nspin),
        "top_Pt_d_C_p": write_cohp_for_pair(out_dir, cohp_dir, "top_Pt_d_C_p", top_pt["shells"]["d"], carbon["shells"]["p"], nspin),
        "top_Pt_d_C_s": write_cohp_for_pair(out_dir, cohp_dir, "top_Pt_d_C_s", top_pt["shells"]["d"], carbon["shells"]["s"], nspin),
    }

    stru_for_geom = out_dir / "STRU.cif"
    if (final_dir / "STRU").exists():
        symbols, positions = parse_stru_positions(final_dir / "STRU")
    else:
        symbols, positions = [], np.empty((0, 3))
    pt_c = float(np.linalg.norm(positions[mapping["top_pt_atom_index"]] - positions[mapping["carbon_atom_index"]])) if len(positions) else float("nan")
    c_o = float(np.linalg.norm(positions[mapping["carbon_atom_index"]] - positions[mapping["oxygen_atom_index"]])) if len(positions) else float("nan")
    del stru_for_geom, symbols

    report = {
        "nspin": nspin,
        "out_dir": str(out_dir),
        "efermi_ev": read_fermi(out_dir),
        "top_pt_c_distance_angstrom": pt_c,
        "c_o_distance_angstrom": c_o,
        "metrics": metrics,
        "interpretation": [
            "Use negative ICOHP in the native COHP sign convention as stronger occupied Pt-C bonding.",
            "Compare nspin=1 and nspin=2 total Pt-C -ICOHP to judge whether spin polarization changes the Pt-C bond.",
            "Compare Pt-d/C-p and Pt-d/C-s components to separate back-donation-like and donation-like contributions qualitatively.",
        ],
    }
    (cohp_dir / "summary.json").write_text(json.dumps(report, indent=2))
    md = [
        f"# Pt(111)-CO nspin={nspin} COHP summary",
        "",
        f"- E_Fermi: {report['efermi_ev']:.6f} eV",
        f"- Top Pt-C distance: {pt_c:.4f} Angstrom",
        f"- C-O distance: {c_o:.4f} Angstrom",
        "",
        "| Channel | Spin | ICOHP | -ICOHP | min COHP | max COHP |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for channel, spin_data in metrics.items():
        for spin, data in spin_data.items():
            md.append(
                f"| {channel} | {spin} | {data['icohp']:.6f} | {data['minus_icohp']:.6f} | "
                f"{data['min_cohp']:.6f} | {data['max_cohp']:.6f} |"
            )
    md += [
        "",
        "Native COHP sign is retained in the numeric ICOHP columns; plots use -COHP.",
        "More positive -ICOHP indicates stronger occupied bonding contribution within this ABACUS-NAO COHP definition.",
    ]
    (cohp_dir / "summary.md").write_text("\n".join(md) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(REPO_ROOT / "runs" / "pt111_co_top"))
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("prepare")
    pf = sub.add_parser("prepare-final")
    pf.add_argument("--nspin", type=int, choices=[1, 2], required=True)
    an = sub.add_parser("analyze")
    an.add_argument("--nspin", type=int, choices=[1, 2], required=True)
    args = parser.parse_args()
    root = Path(args.root)
    if args.cmd == "prepare":
        prepare(root)
    elif args.cmd == "prepare-final":
        prepare_final(root, args.nspin)
    elif args.cmd == "analyze":
        analyze(root, args.nspin)


if __name__ == "__main__":
    main()
