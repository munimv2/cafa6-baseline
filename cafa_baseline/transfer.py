"""BLAST/DIAMOND annotation-transfer baseline.

This is the "BLAST" baseline popularised by DeepGOPlus (Kulmanov & Hoehndorf,
Bioinformatics 2020): for a query protein, the score of a GO term is the sum of
bit-scores of the hits annotated with that term, divided by the sum of
bit-scores of all hits.  Scores therefore lie in ``(0, 1]`` and are then
propagated up the GO DAG.

You run DIAMOND/BLAST yourself (CPU, no GPU) to produce the hits table; this
module only turns that table plus the training annotations into scores, so it
has no third-party dependencies and is fully unit-tested on synthetic data.
"""

from typing import Dict, Iterable, List, Set, Tuple, Union

from .propagate import propagate_scores

Hits = Union[
    Dict[str, Iterable[Tuple[str, float]]],
    Iterable[Tuple[str, str, float]],
]


def _as_query_map(hits: Hits) -> Dict[str, List[Tuple[str, float]]]:
    if isinstance(hits, dict):
        return {q: list(hitlist) for q, hitlist in hits.items()}
    by_query: Dict[str, List[Tuple[str, float]]] = {}
    for query, target, bitscore in hits:
        by_query.setdefault(query, []).append((target, bitscore))
    return by_query


def transfer_scores(
    hits: Hits,
    target_annotations: Dict[str, Iterable[str]],
    ancestors: Dict[str, Set[str]] = None,
    propagate: bool = True,
) -> Dict[str, Dict[str, float]]:
    """Score GO terms for each query from its alignment hits.

    ``hits`` is either ``{query: [(target, bitscore), ...]}`` or an iterable of
    ``(query, target, bitscore)`` rows.  ``target_annotations`` maps each
    (training) target protein to its GO terms.  If ``propagate`` is True and
    ``ancestors`` is given, scores are propagated up the GO DAG.
    """
    by_query = _as_query_map(hits)
    scored: Dict[str, Dict[str, float]] = {}
    for query, hitlist in by_query.items():
        total = sum(b for _, b in hitlist if b > 0)
        term_score: Dict[str, float] = {}
        if total <= 0:
            scored[query] = term_score
            continue
        for target, bitscore in hitlist:
            if bitscore <= 0:
                continue
            weight = bitscore / total
            for term in set(target_annotations.get(target, ())):
                term_score[term] = term_score.get(term, 0.0) + weight
        scored[query] = term_score
    if propagate and ancestors is not None:
        scored = propagate_scores(scored, ancestors, include_self=True)
    return scored
