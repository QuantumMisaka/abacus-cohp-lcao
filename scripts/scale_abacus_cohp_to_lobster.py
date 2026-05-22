#!/usr/bin/env python3
"""Scale ABACUS LCAO-COHP curves onto an empirical LOBSTER-like magnitude."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


SCALE_BASIS = "empirical_lobster_minus_icohp_over_abacus_minus_icohp"
HELP_EPILOG = """\
Examples:
  List available empirical presets:
    python scripts/scale_abacus_cohp_to_lobster.py --list-presets

  Scale a raw COHP file whose energy axis is absolute eV:
    python scripts/scale_abacus_cohp_to_lobster.py si_si_COHP.dat \\
      --preset Si-Si \\
      --efermi 7.111283804 \\
      --output-prefix si_si_lobster_like

  Scale a file that is already plotted as -COHP on an E-E_F axis:
    python scripts/scale_abacus_cohp_to_lobster.py pair_minus_cohp.dat \\
      --preset Pt-C \\
      --input-convention minus-cohp \\
      --efermi 0.0

Output:
  OUTPUT.dat columns are:
    energy_ev source_value minus_cohp lobster_like_minus_cohp
  OUTPUT.json records the preset, scale factor, input convention, integrated
  ABACUS -ICOHP, and scaled LOBSTER-like -ICOHP.

Scientific boundary:
  This is an empirical readability scale based on completed benchmark ratios.
  It does not prove ABACUS NAO-COHP and LOBSTER pCOHP are the same observable.
"""


@dataclass(frozen=True)
class PresetScale:
    scale_factor: float
    abacus_minus_icohp: float
    lobster_minus_icohp: float
    source: str
    note: str


@dataclass(frozen=True)
class ScaleResult:
    energy: np.ndarray
    source_cohp: np.ndarray
    minus_cohp: np.ndarray
    lobster_like_minus_cohp: np.ndarray
    preset: str
    scale_factor: float
    abacus_minus_icohp: float
    lobster_like_minus_icohp: float


PRESET_SCALES: dict[str, PresetScale] = {
    "Si-Si": PresetScale(
        scale_factor=74.5722164489454,
        abacus_minus_icohp=0.06198461330654695,
        lobster_minus_icohp=4.62233,
        source="Si2 ABACUS LTS 3.10.x vs VASP+LOBSTER Bunge",
        note="Full-window Si-Si -ICOHP ratio from the Si2 validation comparison.",
    ),
    "C-C-sp": PresetScale(
        scale_factor=112.15722194365841,
        abacus_minus_icohp=0.0854997996011099,
        lobster_minus_icohp=9.58942,
        source="Diamond ABACUS APNS precision C(sp)-C(sp) vs VASP+LOBSTER Bunge",
        note="Preferred diamond C-C preset when the ABACUS channel uses sp-like orbitals.",
    ),
    "C-C-all": PresetScale(
        scale_factor=166.88012210768915,
        abacus_minus_icohp=0.057462925355554736,
        lobster_minus_icohp=9.58942,
        source="Diamond ABACUS APNS precision C(all)-C(all) vs VASP+LOBSTER Bunge",
        note="Diagnostic diamond C-C preset for all ABACUS NAOs, including polarization orbitals.",
    ),
    "Pt-C": PresetScale(
        scale_factor=41.30589337731911,
        abacus_minus_icohp=0.11649766187219837,
        lobster_minus_icohp=4.81204,
        source="Pt(111)-CO top Pt-C ABACUS nspin=2 vs VASP+LOBSTER pbeVaspFit2015",
        note="Full-window top-site Pt-C total-channel -ICOHP ratio.",
    ),
    "C-O-Pt111": PresetScale(
        scale_factor=37.6713072780416,
        abacus_minus_icohp=0.48054637091269803,
        lobster_minus_icohp=18.10281,
        source="Pt(111)-CO C-O ABACUS nspin=2 vs VASP+LOBSTER pbeVaspFit2015",
        note="Full-window CO internal bond ratio from the Pt(111)-CO comparison.",
    ),
    "Ni-C": PresetScale(
        scale_factor=57.26049696941833,
        abacus_minus_icohp=0.05721186810078918,
        lobster_minus_icohp=3.27598,
        source="Ni(100)-CO precision-basis Ni-C ABACUS vs VASP+LOBSTER pbeVaspFit2015",
        note="Precision-basis Ni-C total-channel -ICOHP ratio.",
    ),
    "Ni-d_C-p": PresetScale(
        scale_factor=69.06790603856828,
        abacus_minus_icohp=0.04743129172282502,
        lobster_minus_icohp=3.27598,
        source="Ni(100)-CO precision-basis Ni-d/C-p ABACUS vs VASP+LOBSTER pbeVaspFit2015",
        note="Maps the explicit ABACUS Ni-d/C-p channel to the LOBSTER Ni-C total pair scale.",
    ),
    "C-O": PresetScale(
        scale_factor=30.565025365498112,
        abacus_minus_icohp=0.54119176418784,
        lobster_minus_icohp=16.54154,
        source="Ni(100)-CO precision-basis C-O ABACUS vs VASP+LOBSTER pbeVaspFit2015",
        note="Default C-O preset from the spin-polarized Ni(100)-CO comparison.",
    ),
    "C-p_O-p": PresetScale(
        scale_factor=48.31533231929231,
        abacus_minus_icohp=0.3423662677239822,
        lobster_minus_icohp=16.54154,
        source="Ni(100)-CO precision-basis C-p/O-p ABACUS vs VASP+LOBSTER pbeVaspFit2015",
        note="Maps the explicit ABACUS C-p/O-p channel to the LOBSTER C-O total pair scale.",
    ),
}


def load_two_column_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"{path} must contain at least two columns: energy and COHP")
    return np.asarray(data[:, 0], dtype=float), np.asarray(data[:, 1], dtype=float)


def integrate_occupied(energy: np.ndarray, curve: np.ndarray, efermi: float) -> float:
    occupied = energy <= efermi
    if np.count_nonzero(occupied) < 2:
        return float("nan")
    return float(np.trapz(curve[occupied], energy[occupied]))


def scale_curve(
    energy: np.ndarray,
    values: np.ndarray,
    preset: str,
    input_convention: str = "cohp",
    efermi: float = 0.0,
) -> ScaleResult:
    if preset not in PRESET_SCALES:
        valid = ", ".join(sorted(PRESET_SCALES))
        raise ValueError(f"Unknown preset {preset!r}. Valid presets: {valid}")
    if input_convention not in {"cohp", "minus-cohp"}:
        raise ValueError("input_convention must be 'cohp' or 'minus-cohp'")

    energy = np.asarray(energy, dtype=float)
    values = np.asarray(values, dtype=float)
    if energy.shape != values.shape:
        raise ValueError("energy and values must have the same shape")

    minus_cohp = -values if input_convention == "cohp" else values.copy()
    scale_factor = PRESET_SCALES[preset].scale_factor
    lobster_like = minus_cohp * scale_factor
    return ScaleResult(
        energy=energy,
        source_cohp=values,
        minus_cohp=minus_cohp,
        lobster_like_minus_cohp=lobster_like,
        preset=preset,
        scale_factor=scale_factor,
        abacus_minus_icohp=integrate_occupied(energy, minus_cohp, efermi),
        lobster_like_minus_icohp=integrate_occupied(energy, lobster_like, efermi),
    )


def write_scaled_outputs(
    input_path: Path,
    output_prefix: Path,
    preset: str,
    input_convention: str = "cohp",
    efermi: float = 0.0,
) -> dict:
    energy, values = load_two_column_curve(input_path)
    result = scale_curve(
        energy=energy,
        values=values,
        preset=preset,
        input_convention=input_convention,
        efermi=efermi,
    )

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    data_path = output_prefix.with_suffix(".dat")
    json_path = output_prefix.with_suffix(".json")
    np.savetxt(
        data_path,
        np.column_stack(
            [
                result.energy,
                result.source_cohp,
                result.minus_cohp,
                result.lobster_like_minus_cohp,
            ]
        ),
        header="energy_ev source_value minus_cohp lobster_like_minus_cohp",
    )

    preset_info = PRESET_SCALES[preset]
    summary = {
        "input": str(input_path),
        "output_dat": str(data_path),
        "preset": preset,
        "scale_basis": SCALE_BASIS,
        "scale_factor": result.scale_factor,
        "input_convention": input_convention,
        "efermi_ev": efermi,
        "abacus_minus_icohp": result.abacus_minus_icohp,
        "lobster_like_minus_icohp": result.lobster_like_minus_icohp,
        "preset_reference": asdict(preset_info),
        "interpretation_warning": (
            "This is an empirical readability scale, not a proof that ABACUS "
            "NAO-COHP and LOBSTER pCOHP are the same observable."
        ),
    }
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Scale a two-column ABACUS LCAO-COHP .dat curve by an empirical "
            "LOBSTER/ABACUS -ICOHP ratio so the plotted magnitude is easier "
            "to compare with LOBSTER-style reports."
        ),
        epilog=HELP_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", type=Path, nargs="?", help="Two-column ABACUS COHP .dat file.")
    parser.add_argument(
        "--preset",
        choices=sorted(PRESET_SCALES),
        help="Empirical LOBSTER/ABACUS -ICOHP scale preset.",
    )
    parser.add_argument(
        "--input-convention",
        choices=["cohp", "minus-cohp"],
        default="cohp",
        help="Whether the second input column is raw COHP or already -COHP.",
    )
    parser.add_argument(
        "--efermi",
        type=float,
        default=0.0,
        help="Fermi energy in the same energy coordinate as the input file. Use 0 for E-E_F data.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        help="Output prefix for .dat and .json. Defaults to INPUT stem plus _lobster_like.",
    )
    parser.add_argument("--list-presets", action="store_true", help="Print preset scale table and exit.")
    args = parser.parse_args()

    if args.list_presets:
        for name, preset in sorted(PRESET_SCALES.items()):
            print(
                f"{name}\t{preset.scale_factor:.8g}\t"
                f"LOBSTER {preset.lobster_minus_icohp:.6g} / ABACUS {preset.abacus_minus_icohp:.6g}\t"
                f"{preset.source}"
            )
        return

    if args.input is None:
        parser.error("input is required unless --list-presets is used")
    if args.preset is None:
        parser.error("--preset is required unless --list-presets is used")

    output_prefix = args.output_prefix or args.input.with_name(args.input.stem + "_lobster_like")
    summary = write_scaled_outputs(
        input_path=args.input,
        output_prefix=output_prefix,
        preset=args.preset,
        input_convention=args.input_convention,
        efermi=args.efermi,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
