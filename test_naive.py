import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cafa_baseline.naive import naive_frequencies, predict_naive

ANCESTORS = {
    "GO:0000001": set(),
    "GO:0000002": {"GO:0000001"},
    "GO:0000003": {"GO:0000002", "GO:0000001"},
}


def test_frequencies_include_propagated_parents():
    train = {"A": ["GO:0000003"], "B": ["GO:0000002"]}
    freqs = naive_frequencies(train, ANCESTORS)
    # A -> {3,2,1}; B -> {2,1}; over 2 proteins
    assert freqs["GO:0000001"] == 1.0
    assert freqs["GO:0000002"] == 1.0
    assert freqs["GO:0000003"] == 0.5


def test_empty_training_set():
    assert naive_frequencies({}, ANCESTORS) == {}


def test_predict_assigns_same_frequencies_to_each_protein():
    freqs = {"GO:0000001": 0.5}
    preds = predict_naive(["X", "Y"], freqs)
    assert preds["X"] == {"GO:0000001": 0.5}
    assert preds["Y"] == {"GO:0000001": 0.5}
    # must be independent copies
    preds["X"]["GO:0000001"] = 0.1
    assert preds["Y"]["GO:0000001"] == 0.5


if __name__ == "__main__":
    test_frequencies_include_propagated_parents()
    test_empty_training_set()
    test_predict_assigns_same_frequencies_to_each_protein()
    print("test_naive OK")
