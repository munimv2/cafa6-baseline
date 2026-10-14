# cafa6-baseline

A small, **dependency-free, CPU-only** toolkit for CAFA-style protein function
prediction — predicting [Gene Ontology](https://geneontology.org/) (GO) terms for
proteins from sequence. It gives you the reusable pieces every
[CAFA](https://biofunctionprediction.org/cafa/) entry needs, with no third-party
dependencies and no GPU required:

| Module | What it does |
| --- | --- |
| `cafa_baseline.ontology` | Parse a GO `.obo` file; build the ancestor closure over `is_a` and `part_of`. |
| `cafa_baseline.propagate` | Apply the **true-path rule**: a term's score propagates up to all ancestors (max). |
| `cafa_baseline.naive` | The classic CAFA **Naive** baseline: per-term background frequency from a training set. |
| `cafa_baseline.transfer` | The DeepGOPlus-style **BLAST/DIAMOND annotation-transfer** baseline, from a hits table. |
| `cafa_baseline.submission` | Blend predictors, cap per protein, and write a CAFA-format prediction file. |

Everything is pure Python (standard library only), so it installs and runs
anywhere and every function is unit-tested on small synthetic data.

> **Status of the Kaggle CAFA 6 competition:** it is closed. Final submissions
> were due 2026-02-02 and the final evaluation ended 2026-06-01
> ([competition page](https://www.kaggle.com/competitions/cafa-6-protein-function-prediction)).
> This toolkit is for CAFA-style baselines and evaluation splits on public
> UniProt/GOA data, and for future CAFA rounds. It is not an entry to CAFA 6.

## Why

A strong CAFA baseline does **not** need a GPU. Sequence-similarity annotation
transfer plus the naive frequency prior gives a respectable baseline, and both
run on CPU. The expensive, GPU-bound part (generating protein language
model embeddings) is optional and separate — this toolkit covers the parts you
can run on a laptop and is the scaffolding you extend from there.

## Install

```bash
pip install .      # or just copy the cafa_baseline/ folder next to your code
```

Requires Python 3.8+. No runtime dependencies.

## Pipeline (how the pieces fit)

1. **Download the public inputs** (no GPU): the GO `.obo` file, the training
   protein→GO annotations, and the test protein list (public UniProt/GOA data).
2. **Build the ontology:**
   ```python
   from cafa_baseline import parse_obo, build_ancestors
   terms = parse_obo(open("go-basic.obo").read())
   ancestors = build_ancestors(terms)
   ```
3. **Naive baseline** (the floor every method is compared against):
   ```python
   from cafa_baseline import naive_frequencies, predict_naive
   freqs = naive_frequencies(train_annotations, ancestors)
   naive = predict_naive(test_proteins, freqs)
   ```
4. **BLAST/DIAMOND annotation transfer.** Run DIAMOND yourself on CPU, e.g.
   ```bash
   diamond makedb --in train.fasta -d train
   diamond blastp -q test.fasta -d train --outfmt 6 qseqid sseqid bitscore -o hits.tsv
   ```
   then score and propagate:
   ```python
   from cafa_baseline import transfer_scores
   hits = [(q, s, float(b)) for q, s, b in (ln.split("\t") for ln in open("hits.tsv"))]
   transfer = transfer_scores(hits, train_annotations, ancestors)
   ```
5. **Blend, cap and write the predictions:**
   ```python
   from cafa_baseline import blend, write_submission
   final = blend([naive, transfer], weights=[0.2, 0.8])
   write_submission(final, "predictions.tsv", top_k=500)
   ```
   The writer defaults to tab-separated, three columns (`protein  GO_term  score`),
   no header — the CAFA convention. Delimiter and header are configurable.

See [`examples/run_example.py`](examples/run_example.py) for a complete, runnable
end-to-end demo on tiny synthetic data (no downloads).

## Tests

```bash
python test_ontology.py && python test_propagate.py && python test_naive.py \
  && python test_transfer.py && python test_submission.py
# or: pytest
```

## Method notes & references

- True-path propagation and IA-weighted Fmax scoring follow the CAFA protocol
  (Radivojac et al.; Zhou et al., *Genome Biology* 2019).
- The bit-score-weighted annotation-transfer baseline is the "BLAST" baseline
  from DeepGOPlus (Kulmanov & Hoehndorf, *Bioinformatics* 2020).
- For scoring your own validation split, the official evaluator is
  [`cafaeval`](https://pypi.org/project/cafaeval/) (CAFA-evaluator).

## Who made this

Made by **MV2 of Munim, Inc.**, with the code written by AI agents working in his
name. If you use it in a competition, disclose AI-assisted work wherever that
competition's rules require it.

## License

Apache-2.0 — see [LICENSE](LICENSE).
