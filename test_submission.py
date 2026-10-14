import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cafa_baseline.submission import blend, cap_predictions, format_rows, write_submission


def approx(a, b, tol=1e-9):
    return abs(a - b) < tol


def test_blend_is_weighted_average():
    a = {"P": {"GO:1": 1.0}}
    b = {"P": {"GO:1": 0.0, "GO:2": 1.0}}
    out = blend([a, b])  # equal weights
    assert approx(out["P"]["GO:1"], 0.5)
    assert approx(out["P"]["GO:2"], 0.5)


def test_blend_respects_weights():
    a = {"P": {"GO:1": 1.0}}
    b = {"P": {"GO:1": 0.0}}
    out = blend([a, b], weights=[3.0, 1.0])
    assert approx(out["P"]["GO:1"], 0.75)


def test_cap_keeps_top_k_above_min_score():
    preds = {"P": {"GO:1": 0.9, "GO:2": 0.5, "GO:3": 0.0005}}
    out = cap_predictions(preds, top_k=1, min_score=0.001)
    assert out["P"] == {"GO:1": 0.9}


def test_format_rows_sorted_and_rounded():
    preds = {"B": {"GO:2": 0.4}, "A": {"GO:1": 0.12345, "GO:9": 0.00001}}
    rows = format_rows(preds, top_k=500, min_score=0.001, decimals=3)
    # proteins sorted; GO:9 dropped (below min_score); score rounded to 3 dp
    assert rows == [("A", "GO:1", "0.123"), ("B", "GO:2", "0.400")]


def test_write_submission_tsv_no_header_roundtrip():
    preds = {"P1": {"GO:0000001": 0.5, "GO:0000002": 1.0}}
    path = os.path.join(tempfile.mkdtemp(), "submission.tsv")
    n = write_submission(preds, path)
    assert n == 2
    with open(path, encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    assert lines[0] == "P1\tGO:0000002\t1.000"  # highest score first
    assert lines[1] == "P1\tGO:0000001\t0.500"


def test_write_submission_optional_header_and_delimiter():
    preds = {"P1": {"GO:0000001": 0.5}}
    path = os.path.join(tempfile.mkdtemp(), "submission.csv")
    write_submission(preds, path, delimiter=",", header=("Protein", "Term", "Score"))
    with open(path, encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    assert lines[0] == "Protein,Term,Score"
    assert lines[1] == "P1,GO:0000001,0.500"


if __name__ == "__main__":
    test_blend_is_weighted_average()
    test_blend_respects_weights()
    test_cap_keeps_top_k_above_min_score()
    test_format_rows_sorted_and_rounded()
    test_write_submission_tsv_no_header_roundtrip()
    test_write_submission_optional_header_and_delimiter()
    print("test_submission OK")
