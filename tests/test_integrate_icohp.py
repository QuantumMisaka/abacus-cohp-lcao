from pathlib import Path
import json
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import cohp  # noqa: E402
import read_abacus_out as rao  # noqa: E402


def test_run_outdir_rejects_incomplete_projected_cohp_modes(tmp_path: Path):
    with pytest.raises(NotImplementedError, match="pCOHP/pCOOP are disabled"):
        cohp.run_outdir(tmp_path, [0], [1], testmethod="pCOHP")


def test_atom_pair_weight_uses_ev_and_both_hermitian_directions():
    hamiltonian = np.array([[0.0, 0.5], [0.5, 0.0]], dtype=np.complex128)
    eigenvectors = np.ones((2, 1), dtype=np.complex128) / np.sqrt(2.0)
    energies, weights = cohp.cal_COHPvalsIJ_e(
        Hks=[hamiltonian],
        Sks=[np.eye(2)],
        Eks=[np.array([-0.25])],
        Cks=[eigenvectors],
        wk=[1.0],
        atomI_orbs=[0],
        atomJ_orbs=[1],
    )
    np.testing.assert_allclose(energies, [-0.25 * rao.RY_TO_EV])
    np.testing.assert_allclose(weights, [0.5 * rao.RY_TO_EV])


def test_discrete_icohp_and_broadening_preserve_state_weights():
    energy = np.array([-0.5, 0.5])
    weights = np.array([-0.4, 0.7])
    result = cohp.integrate_discrete_icohp(energy, weights, efermi=0.0)
    np.testing.assert_allclose(result["minus_icohp"], 0.4)
    for de in (0.2, 0.05):
        grid, spectrum = cohp.broaden_discrete_spectrum(energy, weights, de=de, sigma=0.15)
        np.testing.assert_allclose(np.trapezoid(spectrum, grid), weights.sum(), atol=2e-4)


def test_integrate_icohp_reports_native_and_minus_conventions():
    energy = np.array([-1.0, 0.0, 1.0])
    cohp_values = np.array([-0.1, -0.2, 0.3])

    result = cohp.integrate_icohp(energy, cohp_values, efermi=0.0, input_convention="cohp")
    minus_result = cohp.integrate_icohp(energy, -cohp_values, efermi=0.0, input_convention="minus-cohp")

    assert result["points_occupied"] == 2
    assert result["uses_plot_smoothing"] is False
    np.testing.assert_allclose(result["icohp"], -0.15)
    np.testing.assert_allclose(result["minus_icohp"], 0.15)
    np.testing.assert_allclose(minus_result["minus_icohp"], 0.15)


def test_integrate_icohp_cli_reads_metadata_efermi_and_writes_json(tmp_path: Path):
    curve = tmp_path / "pair.dat"
    curve.write_text("# Energy(eV) COHP\n4.0 -0.1\n5.0 -0.2\n6.0 0.3\n")
    (tmp_path / "pair.meta.json").write_text(json.dumps({"efermi_ev": 5.0}) + "\n")
    output = tmp_path / "pair_icohp.json"

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "integrate_icohp.py"),
            str(curve),
            "--output-json",
            str(output),
        ],
        check=True,
        text=True,
        capture_output=True,
    )

    summary = json.loads(output.read_text())
    np.testing.assert_allclose(summary["minus_icohp"], 0.15)
    assert summary["efermi_ev"] == 5.0
    assert summary["input_convention"] == "cohp"
    assert "-ICOHP = 0.150000 eV" in result.stdout


def test_integrate_icohp_cli_prefers_exact_discrete_metadata(tmp_path: Path):
    curve = tmp_path / "pair.dat"
    curve.write_text("# Energy(eV) COHP\n-1.0 -0.1\n0.0 -0.2\n1.0 0.3\n")
    metadata = {
        "efermi_ev": 0.0,
        "icohp": {
            "efermi_ev": 0.0,
            "icohp": -4.5,
            "minus_icohp": 4.5,
            "integration": "occupied_discrete_state_sum",
        },
    }
    (tmp_path / "pair.meta.json").write_text(json.dumps(metadata) + "\n")

    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "integrate_icohp.py"), str(curve)],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "-ICOHP = 4.500000 eV" in result.stdout
    assert "metadata_discrete_state_sum" in result.stdout


def test_integrate_icohp_cli_defaults_shifted_curve_to_zero_fermi(tmp_path: Path):
    curve = tmp_path / "pair_EminusEf.dat"
    curve.write_text("# Energy-E_Fermi(eV) COHP\n-1.0 -0.1\n0.0 -0.2\n1.0 0.3\n")

    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "integrate_icohp.py"), str(curve)],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "E_Fermi = 0.000000 eV" in result.stdout
    assert "-ICOHP = 0.150000 eV" in result.stdout


def test_scale_cli_rejects_invalid_historical_preset_without_opt_in(tmp_path: Path):
    curve = tmp_path / "pair.dat"
    curve.write_text("-1.0 -0.1\n0.0 -0.2\n")

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "scale_abacus_cohp_to_lobster.py"),
            str(curve),
            "--preset",
            "Si-Si",
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode != 0
    assert "invalid historical presets" in result.stderr
