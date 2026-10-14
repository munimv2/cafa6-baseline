"""End-to-end demo of cafa6-baseline on tiny synthetic data -- no downloads,
no GPU. Run with:  python examples/run_example.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cafa_baseline import (
    parse_obo,
    build_ancestors,
    naive_frequencies,
    predict_naive,
    transfer_scores,
    blend,
    write_submission,
)

HERE = os.path.dirname(os.path.abspath(__file__))

# 1. Ontology
terms = parse_obo(open(os.path.join(HERE, "mini.obo")).read())
ancestors = build_ancestors(terms)

# 2. Training annotations (protein -> leaf GO terms) for two known proteins
train_annotations = {
    "TRAIN1": {"GO:0044237"},  # cellular metabolic process
    "TRAIN2": {"GO:0009987"},  # cellular process
}

# 3. Naive baseline over all test proteins
test_proteins = ["TEST_A", "TEST_B"]
freqs = naive_frequencies(train_annotations, ancestors)
naive = predict_naive(test_proteins, freqs)

# 4. BLAST/DIAMOND annotation transfer from a (synthetic) hits table:
#    columns are (query, target, bitscore) -- as produced by `diamond blastp`.
hits = [
    ("TEST_A", "TRAIN1", 220.0),
    ("TEST_A", "TRAIN2", 80.0),
    ("TEST_B", "TRAIN2", 150.0),
]
transfer = transfer_scores(hits, train_annotations, ancestors)

# 5. Blend, then write a CAFA submission file
final = blend([naive, transfer], weights=[0.2, 0.8])
out_path = os.path.join(HERE, "example_submission.tsv")
n = write_submission(final, out_path, top_k=500)

print("background frequencies:", freqs)
print("transfer scores for TEST_A:", transfer["TEST_A"])
print("wrote %d prediction rows to %s" % (n, out_path))
