"""cafa_baseline -- a small, dependency-free toolkit for CPU-only protein
function prediction baselines (CAFA-style Gene Ontology term prediction).

The package implements the reusable pieces every CAFA entry needs, with no
third-party dependencies and no GPU required:

* ``ontology``   -- parse a Gene Ontology ``.obo`` file and build the ancestor
                    closure (over ``is_a`` and ``part_of``).
* ``propagate``  -- apply the true-path rule, propagating a term's score up to
                    all of its ancestors (max).
* ``naive``      -- the classic CAFA "Naive" baseline: per-term background
                    frequency over a training annotation set.
* ``transfer``   -- the DeepGOPlus-style BLAST/DIAMOND annotation-transfer
                    baseline from a precomputed hits table.
* ``submission`` -- blend predictors and write a CAFA/Kaggle submission file.

Built by MV2 of Munim, Inc. The code in this repository was written with AI
assistance; see the README for how to disclose that when you submit.
"""

from .ontology import GOTerm, parse_obo, build_ancestors
from .propagate import propagate_scores
from .naive import naive_frequencies, predict_naive
from .transfer import transfer_scores
from .submission import blend, cap_predictions, format_rows, write_submission

__version__ = "0.1.0"

__all__ = [
    "GOTerm",
    "parse_obo",
    "build_ancestors",
    "propagate_scores",
    "naive_frequencies",
    "predict_naive",
    "transfer_scores",
    "blend",
    "cap_predictions",
    "format_rows",
    "write_submission",
    "__version__",
]
