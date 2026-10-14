"""Minimal, dependency-free parsing of the Gene Ontology (``.obo``) and
computation of the ancestor closure used for true-path propagation.

Only the fields a CAFA-style function-prediction baseline needs are read:
``id``, ``name``, ``namespace``, ``is_a`` and ``part_of`` relationships, and
``is_obsolete``.  The parser is deliberately small so it has no third-party
dependencies and is easy to audit; it is not a general-purpose OBO library.
"""

from typing import Dict, Set


class GOTerm:
    __slots__ = ("id", "name", "namespace", "parents", "is_obsolete")

    def __init__(self, id, name=None, namespace=None, parents=None, is_obsolete=False):
        self.id = id
        self.name = name
        self.namespace = namespace
        self.parents = set(parents or ())
        self.is_obsolete = is_obsolete

    def __repr__(self):
        return (
            "GOTerm(%r, namespace=%r, parents=%r)"
            % (self.id, self.namespace, sorted(self.parents))
        )


def parse_obo(text: str) -> Dict[str, GOTerm]:
    """Parse OBO-format text into a mapping ``term id -> GOTerm``.

    Obsolete terms are dropped, and parent references to unknown/obsolete
    terms are removed so the resulting graph is self-consistent.
    """
    terms: Dict[str, GOTerm] = {}
    in_term = False
    cur = None

    def flush():
        nonlocal cur
        if cur is not None and cur.id and not cur.is_obsolete:
            terms[cur.id] = cur
        cur = None

    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("[") and line.endswith("]"):
            flush()
            in_term = line == "[Term]"
            cur = GOTerm(id=None) if in_term else None
            continue
        if not in_term or cur is None or not line or line.startswith("!"):
            continue
        if ": " not in line:
            continue
        key, _, val = line.partition(": ")
        # strip trailing " ! human readable name" comments
        val = val.split(" ! ", 1)[0].strip()
        if key == "id":
            cur.id = val
        elif key == "name":
            cur.name = val
        elif key == "namespace":
            cur.namespace = val
        elif key == "is_a":
            cur.parents.add(val)
        elif key == "relationship":
            parts = val.split()
            if len(parts) >= 2 and parts[0] == "part_of":
                cur.parents.add(parts[1])
        elif key == "is_obsolete":
            cur.is_obsolete = val.lower() == "true"
    flush()

    for t in terms.values():
        t.parents = {p for p in t.parents if p in terms}
    return terms


def build_ancestors(terms: Dict[str, GOTerm]) -> Dict[str, Set[str]]:
    """Return ``term id -> set of all ancestor ids`` (not including the term
    itself).  Safe against cycles even though GO is a DAG."""
    cache: Dict[str, Set[str]] = {}

    def ancestors_of(tid: str) -> Set[str]:
        if tid in cache:
            return cache[tid]
        cache[tid] = set()  # guard against revisiting during recursion
        result: Set[str] = set()
        for parent in terms[tid].parents:
            if parent in terms:
                result.add(parent)
                result |= ancestors_of(parent)
        cache[tid] = result
        return result

    return {tid: ancestors_of(tid) for tid in terms}
