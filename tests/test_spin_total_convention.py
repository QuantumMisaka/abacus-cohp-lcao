from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import cohp  # noqa: E402


def write_spin_selection_fixture(out_dir: Path, nspin: int) -> None:
    out_dir.mkdir()
    nk = 1
    nmatrix = nk if nspin == 1 else 2 * nk
    (out_dir / "kpoints").write_text("KPOINTS DIRECT COORDINATES\n1 0.0 0.0 0.0 1.0\n")
    for ik in range(nmatrix):
        (out_dir / f"data-{ik}-H").write_text("\n")
        (out_dir / f"data-{ik}-S").write_text("\n")
        (out_dir / f"WFC_NAO_K{ik + 1}.txt").write_text("\n")


def test_matrix_file_selection_applies_pyatb_spin_total_convention(tmp_path: Path):
    nspin1 = tmp_path / "nspin1"
    write_spin_selection_fixture(nspin1, nspin=1)

    selected = cohp._matrix_file_selection(nspin1, spin="sum")

    assert len(selected) == 1
    np.testing.assert_allclose(selected[0][-1], 2.0)

    nspin2 = tmp_path / "nspin2"
    write_spin_selection_fixture(nspin2, nspin=2)

    sum_selected = cohp._matrix_file_selection(nspin2, spin="sum")
    up_selected = cohp._matrix_file_selection(nspin2, spin="up")
    down_selected = cohp._matrix_file_selection(nspin2, spin="down")

    assert [item[0] for item in sum_selected] == [0, 1]
    assert [item[0] for item in up_selected] == [0]
    assert [item[0] for item in down_selected] == [1]
    np.testing.assert_allclose([item[-1] for item in sum_selected], [1.0, 1.0])
    np.testing.assert_allclose(up_selected[0][-1], 1.0)
    np.testing.assert_allclose(down_selected[0][-1], 1.0)


def test_initialize_from_outdir_applies_pyatb_spin_total_convention(monkeypatch, tmp_path: Path):
    out_dir = tmp_path / "OUT.ABACUS"
    write_spin_selection_fixture(out_dir, nspin=1)

    monkeypatch.setattr(cohp.rao, "read_mat_hs", lambda _path: np.eye(2))
    monkeypatch.setattr(
        cohp.rao,
        "read_lowf",
        lambda _path: (np.eye(2), np.zeros(3), np.array([0.0, 1.0]), np.array([1.0, 1.0])),
    )

    *_unused, wks, _efermi, _atom_i, _atom_j = cohp.initialize_from_outdir(
        out_dir=out_dir,
        atomI_orbs=[0],
        atomJ_orbs=[1],
        spin="sum",
    )

    np.testing.assert_allclose(wks, [2.0])


def test_cohp_output_metadata_records_spin_total_convention(tmp_path: Path):
    metadata = cohp._write_cohp_outputs(
        tmp_path / "pair",
        np.array([0.0]),
        np.array([1.0]),
        0.0,
        nspin=1,
        spin_degeneracy_factor=2.0,
    )

    assert metadata["nspin"] == 1
    assert metadata["spin_degeneracy_factor"] == 2.0
    assert metadata["spin_total_convention"] == cohp.SPIN_TOTAL_CONVENTION
