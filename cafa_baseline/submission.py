"""Blend predictors and write a CAFA / Kaggle submission file.

The CAFA 6 (and CAFA 5) Kaggle submission is a delimited file of
``protein_id  GO_term  score`` rows with scores in ``[0, 1]``; proteins absent
from the file are scored as all-zero.  The exact delimiter and whether a header
row is required should be confirmed against the competition's
``sample_submission`` file -- so the writer here is configurable and defaults
to tab-separated with no header, which is the CAFA convention.
"""

from typing import Dict, List, Sequence, Tuple


def blend(
    prediction_sets: Sequence[Dict[str, Dict[str, float]]],
    weights: Sequence[float] = None,
) -> Dict[str, Dict[str, float]]:
    """Weighted average of several ``protein -> {term: score}`` predictors.

    Weights are normalised to sum to 1, so the blended scores stay in the same
    range as the inputs."""
    if not prediction_sets:
        return {}
    if weights is None:
        weights = [1.0] * len(prediction_sets)
    if len(weights) != len(prediction_sets):
        raise ValueError("weights and prediction_sets must be the same length")
    total_w = sum(weights) or 1.0
    out: Dict[str, Dict[str, float]] = {}
    for preds, w in zip(prediction_sets, weights):
        frac = w / total_w
        for protein, term_scores in preds.items():
            dest = out.setdefault(protein, {})
            for term, score in term_scores.items():
                dest[term] = dest.get(term, 0.0) + frac * score
    return out


def cap_predictions(
    predictions: Dict[str, Dict[str, float]],
    top_k: int = 500,
    min_score: float = 0.001,
) -> Dict[str, Dict[str, float]]:
    """Keep only the ``top_k`` highest-scoring terms per protein that are at
    least ``min_score``.  Returns a new dictionary."""
    out: Dict[str, Dict[str, float]] = {}
    for protein, term_scores in predictions.items():
        items = [(t, s) for t, s in term_scores.items() if s >= min_score]
        items.sort(key=lambda ts: (-ts[1], ts[0]))
        if top_k is not None:
            items = items[:top_k]
        out[protein] = dict(items)
    return out


def format_rows(
    predictions: Dict[str, Dict[str, float]],
    top_k: int = 500,
    min_score: float = 0.001,
    decimals: int = 3,
) -> List[Tuple[str, str, str]]:
    """Return sorted ``(protein, term, formatted_score)`` rows ready to write."""
    rows: List[Tuple[str, str, str]] = []
    for protein in sorted(predictions):
        items = [(t, s) for t, s in predictions[protein].items() if s >= min_score]
        items.sort(key=lambda ts: (-ts[1], ts[0]))
        if top_k is not None:
            items = items[:top_k]
        for term, score in items:
            score = min(1.0, max(0.0, score))
            value = round(score, decimals)
            if value <= 0:
                continue
            rows.append((protein, term, ("%.*f" % (decimals, value))))
    return rows


def write_submission(
    predictions: Dict[str, Dict[str, float]],
    path: str,
    top_k: int = 500,
    min_score: float = 0.001,
    decimals: int = 3,
    delimiter: str = "\t",
    header: Sequence[str] = None,
) -> int:
    """Write a CAFA submission file and return the number of rows written.

    Set ``header`` (e.g. ``("Protein Id", "GO Term Id", "Prediction")``) only if
    the competition's ``sample_submission`` includes one; the CAFA default is
    no header.
    """
    rows = format_rows(predictions, top_k, min_score, decimals)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        if header:
            fh.write(delimiter.join(header) + "\n")
        for protein, term, value in rows:
            fh.write(delimiter.join((protein, term, value)) + "\n")
    return len(rows)
