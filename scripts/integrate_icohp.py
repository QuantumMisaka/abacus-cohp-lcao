#!/usr/bin/env python3
"""Integrate occupied ICOHP from a two-column COHP curve."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
import cohp  # noqa: E402


def load_two_column_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"{path} must contain at least two columns: energy and COHP")
    return np.asarray(data[:, 0], dtype=float), np.asarray(data[:, 1], dtype=float)


def metadata_path_for_curve(path: Path) -> Path:
    if path.stem.endswith("_EminusEf"):
        return path.with_name(path.stem.removesuffix("_EminusEf") + ".meta.json")
    return path.with_suffix(".meta.json")


def load_metadata(path: Path) -> dict:
    metadata_path = metadata_path_for_curve(path)
    if not metadata_path.exists():
        return {}
    return json.loads(metadata_path.read_text())


def infer_efermi(path: Path, explicit_efermi: float | None) -> float:
    if explicit_efermi is not None:
        return float(explicit_efermi)
    if path.stem.endswith("_EminusEf"):
        return 0.0
    metadata = load_metadata(path)
    if "efermi_ev" in metadata:
        return float(metadata["efermi_ev"])
    return 0.0


def exact_metadata_icohp(path: Path) -> dict | None:
    icohp = load_metadata(path).get("icohp")
    if not isinstance(icohp, dict):
        return None
    if icohp.get("integration") != "occupied_discrete_state_sum":
        return None
    if not {"efermi_ev", "icohp", "minus_icohp"}.issubset(icohp):
        return None
    result = dict(icohp)
    result["source"] = "metadata_discrete_state_sum"
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Integrate occupied ICOHP/-ICOHP from a two-column COHP curve.",
    )
    parser.add_argument("curve", type=Path, help="Two-column COHP or -COHP curve file")
    parser.add_argument(
        "--input-convention",
        default="cohp",
        choices=["cohp", "minus-cohp"],
        help="Whether the second column is native COHP or plotted -COHP",
    )
    parser.add_argument("--efermi", type=float, help="Fermi energy in the curve energy units/eV")
    parser.add_argument(
        "--integrate-curve",
        action="store_true",
        help=(
            "Ignore an exact discrete-state ICOHP in adjacent metadata and instead "
            "trapezoid-integrate the broadened curve as a diagnostic approximation"
        ),
    )
    parser.add_argument("--output-json", type=Path, help="Write the ICOHP summary to this JSON file")
    args = parser.parse_args(argv)

    result = None
    if not args.integrate_curve and args.efermi is None and args.input_convention == "cohp":
        result = exact_metadata_icohp(args.curve)
    if result is None:
        energy, values = load_two_column_curve(args.curve)
        efermi = infer_efermi(args.curve, args.efermi)
        result = cohp.integrate_icohp(
            energy=energy,
            values=values,
            efermi=efermi,
            input_convention=args.input_convention,
        )
        result["source"] = "curve_trapezoid_approximation"
    result["input"] = str(args.curve)

    if args.output_json:
        args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print(f"E_Fermi = {result['efermi_ev']:.6f} eV")
    print(f"ICOHP = {result['icohp']:.6f} eV")
    print(f"-ICOHP = {result['minus_icohp']:.6f} eV")
    print(f"Integration source = {result['source']}")
    if args.output_json:
        print(f"json: {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
