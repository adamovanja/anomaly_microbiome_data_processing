"""Tests for the antibiotic metadata of the T1D subcohort."""

import numpy as np
import pandas as pd

from src import process_vatanen19_abx as proc_abx

ID_COL = "ID, E=Espoo, Finland, T=Tartu, Estonia"


def _write_t1d_supp_material(path):
    """Supplementary table with one exposed (E1) and one unexposed (E2) infant."""
    pd.DataFrame(
        {
            ID_COL: ["E1", "E2"],
            "First antibiotic drug between 3 and 6 months": ["Amoxicillin", np.nan],
            "First antibiotic drug between 6 and 12 months": ["Cefaclor", np.nan],
            "Breastfeeding duration": [3, 6],
        }
    ).to_excel(path, sheet_name="Data", index=False)


def test_add_abx_info_t1d_flags_infants_with_a_course_as_ever_exposed(tmp_path):
    path = tmp_path / "mmc2.xlsx"
    _write_t1d_supp_material(path)
    md_both = pd.DataFrame(
        {
            "host_id": ["E1", "E1", "E2", "K1", "K2"],
            "study_subcohort": ["t1d", "t1d", "t1d", "karelia", "karelia"],
            "abx_ever": [np.nan, np.nan, np.nan, True, False],
        },
        index=["s1", "s2", "s3", "s4", "s5"],
    )

    out = proc_abx.add_abx_info_t1d(path, md_both)

    # T1D infants are flagged by their recorded courses, others keep their flag
    assert out["abx_ever"].to_dict() == {
        "s1": True,
        "s2": True,
        "s3": False,
        "s4": True,
        "s5": False,
    }
    assert "abx_ever_t1d" not in out.columns
