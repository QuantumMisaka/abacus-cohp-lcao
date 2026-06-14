from pathlib import Path
import json
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import cohp  # noqa: E402


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
