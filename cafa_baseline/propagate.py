"""True-path propagation of predicted GO-term scores.

CAFA scoring requires that if a protein is predicted to have a term, it is
also (at least as strongly) predicted to have every ancestor of that term.
We implement the standard rule: a parent's score is the maximum over the
scores of itself and all of its descendants.
"""

from typing import Dict, Set


def propagate_scores(
    predictions: Dict[str, Dict[str, float]],
    ancestors: Dict[str, Set[str]],
    include_self: bool = True,
) -> Dict[str, Dict[str, float]]:
    """Propagate scores up the GO DAG using the max rule.

    ``predictions`` maps ``protein -> {term: score}``.  ``ancestors`` is the
    mapping produced by :func:`cafa_baseline.ontology.build_ancestors`.
    Returns a new dictionary; the input is not modified.
    """
    out: Dict[str, Dict[str, float]] = {}
    for protein, term_scores in predictions.items():
        acc: Dict[str, float] = {}
        for term, score in term_scores.items():
            targets = set(ancestors.get(term, ()))
            if include_self:
                targets.add(term)
            for t in targets:
                if t not in acc or score > acc[t]:
                    acc[t] = score
        out[protein] = acc
    return out
