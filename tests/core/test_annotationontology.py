"""Unit tests for AnnotationOntology factories.

Covers from_prd_input (introduced alongside the bulk-reconstruction
endpoint in modelseed-api) plus a regression test for the
get_feature lookup bug in MSBuilder.build_from_annotaton_ontology.
"""

from __future__ import annotations

import pytest

from modelseedpy.core.annotationontology import (
    AnnotationOntology,
    AnnotationOntologyFeature,
)


# A trivial translator: maps a small set of namespaced terms to a
# pre-known list of MSRXN ids. Anything else returns []. This lets
# the tests assert on msrxns without dragging in cb_annotation_ontology_api
# data files.
def _fake_translator(term):
    table = {
        "KO:K00001": ["MSRXN:rxn00001", "MSRXN:rxn00002"],
        "EC:1.1.1.1": ["MSRXN:rxn00001"],
        "SSO:SSO_alcohol_dehydrogenase": ["MSRXN:rxn00001"],
    }
    return table.get(term, [])


def _empty_translator(term):
    return []


# ─────────────────────────────────────────────────────────────────────
# from_prd_input
# ─────────────────────────────────────────────────────────────────────


def test_from_prd_input_builds_minimal_genome(tmp_path):
    annotations = {
        "gene1": {"KO": [{"term": "K00001", "score": 0.9}]},
    }
    ao = AnnotationOntology.from_prd_input(
        "test-genome",
        annotations,
        data_dir=str(tmp_path),
        translator=_fake_translator,
    )
    assert ao.genome_ref == "test-genome"
    assert "gene1" in ao.genes
    assert ao.feature_types["gene1"] == "gene"
    # One event per ontology type
    assert len(ao.events) == 1
    assert ao.events[0].ontology.id == "KO"
    # Term registered with translated reactions
    assert "K00001" in ao.terms
    assert ao.terms["K00001"].msrxns == {"rxn00001", "rxn00002"}


def test_from_prd_input_multiple_genes_multiple_ontologies(tmp_path):
    annotations = {
        "geneA": {
            "KO": [{"term": "K00001", "score": 0.8}],
            "EC": [{"term": "1.1.1.1", "score": 0.6}],
        },
        "geneB": {
            "SSO": [{"term": "SSO_alcohol_dehydrogenase", "score": 1.0}],
        },
    }
    ao = AnnotationOntology.from_prd_input(
        "g",
        annotations,
        data_dir=str(tmp_path),
        translator=_fake_translator,
    )
    # Three events: KO, EC, SSO
    assert {e.ontology.id for e in ao.events} == {"KO", "EC", "SSO"}
    # Both features registered
    assert {"geneA", "geneB"} <= set(ao.genes.keys())
    # Each ontology event linked to the right features
    ko = next(e for e in ao.events if e.ontology.id == "KO")
    sso = next(e for e in ao.events if e.ontology.id == "SSO")
    assert "geneA" in ko.features and "geneB" not in ko.features
    assert "geneB" in sso.features and "geneA" not in sso.features


def test_from_prd_input_records_score_as_probability(tmp_path):
    ao = AnnotationOntology.from_prd_input(
        "g",
        {"gene1": {"KO": [{"term": "K00001", "score": 0.42}]}},
        data_dir=str(tmp_path),
        translator=_fake_translator,
    )
    feature = ao.genes["gene1"]
    event_id = list(feature.event_terms.keys())[0]
    evidence = feature.event_terms[event_id]["K00001"]
    assert evidence.probability == 0.42
    assert evidence.scores == {"probability": 0.42}


def test_from_prd_input_default_score_is_one(tmp_path):
    ao = AnnotationOntology.from_prd_input(
        "g",
        {"gene1": {"KO": [{"term": "K00001"}]}},
        data_dir=str(tmp_path),
        translator=_fake_translator,
    )
    feature = ao.genes["gene1"]
    event_id = list(feature.event_terms.keys())[0]
    evidence = feature.event_terms[event_id]["K00001"]
    assert evidence.probability == 1.0


def test_from_prd_input_unmapped_term_retained_with_empty_msrxns(tmp_path):
    """When the translator returns [], the term must still be retained
    in the AnnotationOntology with an empty msrxns set. This matches the
    PRD requirement that unmapped genes never silently disappear."""
    ao = AnnotationOntology.from_prd_input(
        "g",
        {"gene1": {"KO": [{"term": "K99999", "score": 0.5}]}},
        data_dir=str(tmp_path),
        translator=_empty_translator,
    )
    assert "gene1" in ao.genes
    assert "K99999" in ao.terms
    assert ao.terms["K99999"].msrxns == set()


def test_from_prd_input_passes_namespaced_term_to_translator(tmp_path):
    seen = []
    def remember(term):
        seen.append(term)
        return []
    AnnotationOntology.from_prd_input(
        "g",
        {"gene1": {"KO": [{"term": "K00001"}], "EC": [{"term": "1.1.1.1"}]}},
        data_dir=str(tmp_path),
        translator=remember,
    )
    assert "KO:K00001" in seen
    assert "EC:1.1.1.1" in seen


def test_from_prd_input_empty_input_produces_empty_ontology(tmp_path):
    ao = AnnotationOntology.from_prd_input(
        "g", {}, data_dir=str(tmp_path), translator=_fake_translator,
    )
    assert ao.genome_ref == "g"
    assert len(ao.genes) == 0
    assert len(ao.terms) == 0
    assert len(ao.events) == 0


# ─────────────────────────────────────────────────────────────────────
# Regression test: msbuilder.py:789 anno_ont.get_feature lookup
# ─────────────────────────────────────────────────────────────────────
# AnnotationOntology never had a `get_feature` method; the feature is
# keyed in `genes` or `cdss`. The MSBuilder pre-fix code would AttributeError
# the moment it tried to attach probability to a built reaction. The fix
# uses `genes.get(id) or cdss.get(id)` directly.


def test_anno_ont_feature_accessible_via_genes_dict(tmp_path):
    ao = AnnotationOntology.from_prd_input(
        "g",
        {"gene1": {"KO": [{"term": "K00001", "score": 0.5}]}},
        data_dir=str(tmp_path),
        translator=_fake_translator,
    )
    # The fix uses this exact accessor chain.
    looked_up = ao.genes.get("gene1") or ao.cdss.get("gene1")
    assert isinstance(looked_up, AnnotationOntologyFeature)
    assert looked_up.id == "gene1"


def test_anno_ont_has_no_get_feature_method():
    """Lock the bug shape so a future re-introduction is obvious.

    If someone adds a `get_feature` method later they should also update
    msbuilder.py:789's accessor to use it and remove this test. As of
    now (fix commit), the only safe path is direct dict lookup.
    """
    ao = AnnotationOntology(genome_ref="g", data_dir="/tmp")
    assert not hasattr(ao, "get_feature"), (
        "AnnotationOntology gained a get_feature method; reconcile "
        "msbuilder.py:789 to use it (and remove this guard)."
    )
