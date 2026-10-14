"""The CAFA "Naive" baseline.

For every test protein, predict each GO term with a score equal to its
background frequency in the training set (the fraction of training proteins
annotated with that term, after true-path propagation).  This is the standard
floor that every CAFA method is compared against.
"""

from typing import Dict, Iterable, Set


def naive_frequencies(
    train_annotations: Dict[str, Iterable[str]],
    ancestors: Dict[str, Set[str]],
) -> Dict[str, float]:
    """Return ``term -> background frequency`` over the training set.

    ``train_annotations`` maps each training protein to the GO terms it is
    annotated with (leaf or any level); each protein's terms are propagated to
    their ancestors before counting, per the true-path rule.
    """
    n = len(train_annotations)
    if n == 0:
        return {}
    counts: Dict[str, int] = {}
    for terms in train_annotations.values():
        full: Set[str] = set()
        for term in terms:
            full.add(term)
            full |= ancestors.get(term, set())
        for t in full:
            counts[t] = counts.get(t, 0) + 1
    return {t: c / n for t, c in counts.items()}


def predict_naive(
    test_proteins: Iterable[str],
    frequencies: Dict[str, float],
) -> Dict[str, Dict[str, float]]:
    """Assign the same background ``frequencies`` to every test protein."""
    return {p: dict(frequencies) for p in test_proteins}
