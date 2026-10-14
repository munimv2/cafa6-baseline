import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cafa_baseline.ontology import parse_obo, build_ancestors

MINI = """
format-version: 1.2

[Term]
id: GO:0000001
name: root
namespace: biological_process

[Term]
id: GO:0000002
name: child
namespace: biological_process
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: grandchild
namespace: biological_process
is_a: GO:0000002 ! child
relationship: part_of GO:0000001 ! root

[Term]
id: GO:0000009
name: obsolete thing
namespace: biological_process
is_a: GO:0000001 ! root
is_obsolete: true

[Typedef]
id: part_of
name: part of
"""


def test_parses_terms_and_skips_obsolete_and_typedef():
    terms = parse_obo(MINI)
    assert set(terms) == {"GO:0000001", "GO:0000002", "GO:0000003"}
    assert terms["GO:0000002"].namespace == "biological_process"
    assert terms["GO:0000002"].name == "child"
    assert terms["GO:0000002"].parents == {"GO:0000001"}
    # is_a plus part_of both count as parents
    assert terms["GO:0000003"].parents == {"GO:0000002", "GO:0000001"}


def test_ancestor_closure():
    terms = parse_obo(MINI)
    anc = build_ancestors(terms)
    assert anc["GO:0000001"] == set()
    assert anc["GO:0000002"] == {"GO:0000001"}
    assert anc["GO:0000003"] == {"GO:0000002", "GO:0000001"}


def test_parents_to_obsolete_are_dropped():
    terms = parse_obo(MINI)
    for t in terms.values():
        assert "GO:0000009" not in t.parents


if __name__ == "__main__":
    test_parses_terms_and_skips_obsolete_and_typedef()
    test_ancestor_closure()
    test_parents_to_obsolete_are_dropped()
    print("test_ontology OK")
