import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cafa_baseline.transfer import transfer_scores

ANCESTORS = {
    "GO:0000001": set(),
    "GO:0000002": {"GO:0000001"},
    "GO:0000003": {"GO:0000002", "GO:0000001"},
}
TARGET_ANNOT = {"T1": {"GO:0000002"}, "T2": {"GO:0000003"}}


def approx(a, b, tol=1e-9):
    return abs(a - b) < tol


def test_bitscore_weighted_transfer_without_propagation():
    hits = {"Q": [("T1", 100.0), ("T2", 50.0)]}
    out = transfer_scores(hits, TARGET_ANNOT, ancestors=None, propagate=False)
    assert approx(out["Q"]["GO:0000002"], 100.0 / 150.0)
    assert approx(out["Q"]["GO:0000003"], 50.0 / 150.0)


def test_transfer_with_propagation():
    hits = [("Q", "T1", 100.0), ("Q", "T2", 50.0)]  # row form
    out = transfer_scores(hits, TARGET_ANNOT, ancestors=ANCESTORS, propagate=True)
    assert approx(out["Q"]["GO:0000003"], 50.0 / 150.0)
    # GO:0000002 = max(own 0.667, propagated-from-GO3 0.333) = 0.667
    assert approx(out["Q"]["GO:0000002"], 100.0 / 150.0)
    assert approx(out["Q"]["GO:0000001"], 100.0 / 150.0)


def test_same_term_on_multiple_hits_sums_weights():
    annot = {"T1": {"GO:0000002"}, "T2": {"GO:0000002"}}
    hits = {"Q": [("T1", 100.0), ("T2", 50.0)]}
    out = transfer_scores(hits, annot, propagate=False)
    assert approx(out["Q"]["GO:0000002"], 1.0)


def test_zero_total_bitscore_is_safe():
    out = transfer_scores({"Q": [("T1", 0.0)]}, TARGET_ANNOT, propagate=False)
    assert out["Q"] == {}


if __name__ == "__main__":
    test_bitscore_weighted_transfer_without_propagation()
    test_transfer_with_propagation()
    test_same_term_on_multiple_hits_sums_weights()
    test_zero_total_bitscore_is_safe()
    print("test_transfer OK")
