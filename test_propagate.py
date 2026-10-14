import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cafa_baseline.propagate import propagate_scores

ANCESTORS = {
    "GO:0000001": set(),
    "GO:0000002": {"GO:0000001"},
    "GO:0000003": {"GO:0000002", "GO:0000001"},
}


def test_single_term_propagates_to_all_ancestors():
    out = propagate_scores({"P1": {"GO:0000003": 0.9}}, ANCESTORS)
    assert out["P1"] == {"GO:0000003": 0.9, "GO:0000002": 0.9, "GO:0000001": 0.9}


def test_parent_takes_max_of_descendants():
    out = propagate_scores({"P1": {"GO:0000003": 0.5, "GO:0000002": 0.8}}, ANCESTORS)
    assert out["P1"]["GO:0000003"] == 0.5
    assert out["P1"]["GO:0000002"] == 0.8  # max(0.8, 0.5)
    assert out["P1"]["GO:0000001"] == 0.8  # max(0.5, 0.8)


def test_input_is_not_mutated():
    preds = {"P1": {"GO:0000003": 0.9}}
    propagate_scores(preds, ANCESTORS)
    assert preds == {"P1": {"GO:0000003": 0.9}}


if __name__ == "__main__":
    test_single_term_propagates_to_all_ancestors()
    test_parent_takes_max_of_descendants()
    test_input_is_not_mutated()
    print("test_propagate OK")
