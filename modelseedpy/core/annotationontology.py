# -*- coding: utf-8 -*-
import logging
import re
import time
import json
import sys
import pandas as pd
import cobra
from cobra import DictList
from modelseedpy.core.msgenome import MSGenome

# from builtins import None

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.WARNING
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO

# Class structure
# AnnotationOntology -> Features/Events/Terms/Ontologies
#    AnnotationOntologyOntology -> Events/Terms
#    AnnotationOntologyEvent -> Features/Ontology
#        AnnotationOntologyFeature -> Term+Event->Evidence
#            AnnotationOntologyTerm -> Ontology/Events/Featurs
#            AnnotationOntologyEvidence -> --

allowable_score_types = [
    "probability",
    "evalue",
    "bitscore",
    "identity",
    "qalignstart",
    "qalignstop",
    "salignstart",
    "salignstop",
    "kmerhits",
    "tmscore",
    "rmsd",
    "hmmscore",
    "score"
]

def convert_to_search_role(role):
    role = role.lower()
    role = re.sub("\s","",role)
    role = re.sub("[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]*","",role)
    role = re.sub("\#.*$","",role)
    role = re.sub("\(ec:*\)","",role)
    role = re.sub("[\(\)\[\],-]","",role)
    return role

def split_role(role):
    return re.split("\s*;\s+|\s+[\@\/]\s+",role)

class AnnotationOntologyEvidence:
    def __init__(self, parent, event, term, probability=1, scores={}, ref_entity=None, entity_type=None):
        self.parent = parent
        self.event = event
        self.term = term
        self.probability = probability
        self.ref_entity = ref_entity
        self.entity_type = entity_type
        self.scores = scores
        #for item in self.scores:
            #if item not in allowable_score_types:
                #logger.warning(item + " not an allowable score type!")

    def to_data(self):
        output = {
            "event":self.event.method,
            "term":self.term.id,
            "ontology":self.term.ontology.id,
            "probability":self.probability
        }
        if self.ref_entity:
            output["ref_entity"] = self.ref_entity
        if self.entity_type:
            output["entity_type"] = self.entity_type
        if self.scores:
            output["scores"] = self.scores
        return output


class AnnotationOntologyTerm:
    def __init__(self, parent, term_id, ontology):
        self.id = term_id
        self.parent = parent
        self.ontology = ontology
        self.ontology.add_term(self)
        self.parent.add_term(self)
        self.msrxns = set()
        self.events = {}
        self.features = {}

    def add_msrxns(self, rxn_ids):
        for rxn_id in rxn_ids:
            if rxn_id[0:6] == "MSRXN:":
                rxn_id = rxn_id[6:]
                self.msrxns.update([rxn_id])

    def add_event(self, event):
        self.events[event.id] = event

    def add_feature(self, feature):
        self.features[feature.id] = feature


class AnnotationOntologyOntology:
    def __init__(self, parent, ontology_id):
        self.id = ontology_id
        self.parent = parent
        self.events = {}
        self.terms = {}

    def add_event(self, event):
        self.events[event.id] = event

    def add_term(self, term):
        self.terms[term.id] = term


class AnnotationOntologyFeature:
    def __init__(self, parent, feature_id, type=None):
        self.id = feature_id
        self.parent = parent
        parent.add_feature(self)
        self.type = type
        self.event_terms = {}
        self.term_events = {}

    def add_event_term(self, event, term, scores={}, ref_entity=None, entity_type=None,probability=1):
        if event.id not in self.event_terms:
            self.event_terms[event.id] = {}
        self.event_terms[event.id][term.id] = AnnotationOntologyEvidence(
            self,event,term,probability=probability,scores=scores,ref_entity=ref_entity,entity_type=entity_type
        )
        if term.id not in self.term_events:
            self.term_events[term.id] = {}
        self.term_events[term.id][event.id] = self.event_terms[event.id][term.id]

    def get_associated_terms(
        self,
        prioritized_event_list=None,
        ontologies=None,
        merge_all=False,
        translate_to_rast=False,
    ):
        output = {}
        for term_id in self.term_events:
            term = self.parent.terms[term_id]
            if not ontologies or term.ontology.id in ontologies:
                if merge_all or not prioritized_event_list:
                    for event_id in self.term_events[term_id]:
                        if (
                            not prioritized_event_list
                            or event_id in prioritized_event_list
                        ):
                            if term not in output:
                                output[term] = []
                            output[term].append(
                                self.term_events[term_id][event_id].to_data()
                            )
                else:
                    for event_id in prioritized_event_list:
                        if event_id in self.term_events[term_id]:
                            rxns = self.parent.terms[term_id].msrxns
                            if len(rxns) > 0:
                                if term not in output:
                                    output[term] = []
                                output[term].append(
                                    self.term_events[term_id][event_id].to_data()
                                )
                                break
        return output

    def get_associated_reactions(
        self, prioritized_event_list=None, ontologies=None, merge_all=False
    ):
        output = {}
        for term_id in self.term_events:
            if not ontologies or self.parent.terms[term_id].ontology.id in ontologies:
                if merge_all or not prioritized_event_list:
                    for event_id in self.term_events[term_id]:
                        if (
                            not prioritized_event_list
                            or event_id in prioritized_event_list
                        ):
                            rxns = self.parent.terms[term_id].msrxns
                            for rxn_id in rxns:
                                if rxn_id not in output:
                                    output[rxn_id] = []
                                output[rxn_id].append(
                                    self.term_events[term_id][event_id].to_data()
                                )
                else:
                    for event_id in prioritized_event_list:
                        if event_id in self.term_events[term_id]:
                            rxns = self.parent.terms[term_id].msrxns
                            for rxn_id in rxns:
                                if rxn_id not in output:
                                    output[rxn_id] = []
                                output[rxn_id].append(
                                    self.term_events[term_id][event_id].to_data()
                                )
                            if len(rxns) > 0:
                                break
        return output


class AnnotationOntologyEvent:
    def __init__(
        self,
        parent,
        event_id,
        ontology_id,
        method,
        method_version=None,
        description=None,
        timestamp=None,
    ):
        self.id = event_id
        self.parent = parent
        # Linking ontology
        self.ontology = self.parent.add_ontology(ontology_id)
        self.ontology.add_event(self)
        if not description:
            self.description = ""  # TODO
        else:
            self.description = description
        self.method = method
        self.method_version = method_version
        self.timestamp = timestamp
        self.features = {}

    @staticmethod
    def from_data(data, parent):
        if "method_version" not in data:
            data["method_version"] = None
        if "description" not in data:
            data["description"] = None
        if "timestamp" not in data:
            data["timestamp"] = None
        self = AnnotationOntologyEvent(
            parent,
            data["event_id"],
            data["ontology_id"],
            data["method"],
            data["method_version"],
            data["description"],
            data["timestamp"],
        )
        if "ontology_terms" in data:
            for feature_id in data["ontology_terms"]:
                feature = self.parent.add_feature(feature_id)
                self.add_feature(feature)
                for item in data["ontology_terms"][feature_id]:
                    term = self.parent.add_term(item["term"], self.ontology)
                    scores = {}
                    ref_entity = None
                    entity_type = None
                    if "evidence" in item:
                        if "scores" in item["evidence"]:
                            scores = item["evidence"]["scores"]
                        if "reference" in item["evidence"]:
                            ref_entity = item["evidence"]["reference"][1]
                            entity_type = item["evidence"]["reference"][0]
                    probability = 1/len(data["ontology_terms"][feature_id])
                    feature.add_event_term(self, term, scores, ref_entity, entity_type,probability)
                    if "modelseed_ids" in item:
                        term.add_msrxns(item["modelseed_ids"])
        return self

    def add_feature(self, feature):
        self.features[feature.id] = feature

    def to_data(self):
        data = {
            "event_id": self.event_id,
            "description": self.event_id,
            "ontology_id": self.ontology_id,
            "method": self.method,
            "method_version": self.method_version,
            "timestamp": self.timestamp,
            "ontology_terms": {},
        }
        for feature in self.features:
            data["ontology_terms"][feature] = {"term": None}  # TODO


class AnnotationOntology:
    mdlutls = {}

    @staticmethod
    def from_kbase_data(data, genome_ref=None, data_dir=None):
        self = AnnotationOntology(genome_ref, data_dir)
        if "feature_types" in data:
            self.feature_types = data["feature_types"]
        if "events" in data:
            for event in data["events"]:
                self.events += [AnnotationOntologyEvent.from_data(event, self)]
        return self

    @staticmethod
    def from_prd_input(
        genome_id,
        annotations,
        data_dir,
        translator,
        method="PRD",
        method_version="1.0",
        timestamp=None,
    ):
        """Build an AnnotationOntology from a probabilistic-annotation input dict.

        Mirrors ``from_kbase_data`` for the upstream-translated path: this
        factory is the equivalent entry point for inputs whose ontology
        terms have NOT yet been translated to ModelSEED reaction IDs.
        The supplied ``translator`` callable performs that lookup so the
        factory itself stays decoupled from any specific translation
        backend (KBUtilLib, a fake for tests, etc.).

        Parameters
        ----------
        genome_id:
            Identifier for the genome whose annotations are being loaded.
            Stored on the returned object as ``genome_ref``.
        annotations:
            Mapping ``{gene_id: {ontology_type: [{term, score}, ...]}}``.
            ``ontology_type`` examples: "SSO", "EC", "KO". ``score`` should
            be in [0.0, 1.0] (evidence/probability). A gene may carry
            multiple terms across multiple ontology types.
        data_dir:
            Path to ModelSEED ontology data directory (same value the
            existing ``__init__`` and ``get_term_name`` expect; used by
            downstream consumers, not by this factory).
        translator:
            Callable ``(namespaced_term: str) -> list[str]`` returning
            MSRXN ids for a given namespaced term (e.g. ``"KO:K00001"``).
            ``KBUtilLib.KBAnnotationUtils.translate_term_to_modelseed`` is
            the canonical implementation. Return an empty list when no
            reactions are known for the term; the term is still recorded
            (matches the "retain unmapped" requirement) with an empty
            ``msrxns`` set.
        method, method_version, timestamp:
            Provenance fields applied to the synthesized
            :class:`AnnotationOntologyEvent` records (one event per
            ontology type).
        """
        self = AnnotationOntology(genome_id, data_dir)
        # One synthesized event per ontology type, so downstream
        # priority-list logic can target individual ontologies.
        events_by_ontology = {}
        for gene_id, ont_dict in annotations.items():
            feature = self.add_feature(gene_id)
            self.feature_types[gene_id] = "gene"
            for ontology_type, term_list in ont_dict.items():
                if ontology_type not in events_by_ontology:
                    event = AnnotationOntologyEvent(
                        self,
                        event_id=method + ":" + ontology_type,
                        ontology_id=ontology_type,
                        method=method,
                        method_version=method_version,
                        timestamp=timestamp,
                    )
                    self.events += [event]
                    events_by_ontology[ontology_type] = event
                event = events_by_ontology[ontology_type]
                event.add_feature(feature)
                for term_entry in term_list:
                    term_id = term_entry["term"]
                    score = float(term_entry.get("score", 1.0))
                    term = self.add_term(term_id, event.ontology)
                    feature.add_event_term(
                        event,
                        term,
                        scores={"probability": score},
                        probability=score,
                    )
                    namespaced = ontology_type + ":" + term_id
                    msrxn_ids = translator(namespaced)
                    if msrxn_ids:
                        term.add_msrxns(msrxn_ids)
        return self

    def __init__(self, genome_ref, data_dir):
        self.genome_ref = genome_ref
        self.events = DictList()
        self.terms = {}
        self.ontologies = {}
        self.genes = {}
        self.cdss = {}
        self.data_dir = data_dir
        self.noncodings = {}
        self.feature_types = {}
        self.term_names = {}
        self.info = None

    def get_term_name(self, term):
        if term.ontology.id not in self.term_names:
            self.term_names[term.ontology.id] = {}
            if term.ontology.id in [
                "SSO",
                "AntiSmash",
                "EC",
                "TC",
                "META",
                "RO",
                "KO",
                "GO",
            ]:
                with open(
                    self.data_dir + "/" + term.ontology.id + "_dictionary.json"
                ) as json_file:
                    ontology = json.load(json_file)
                    for item in ontology["term_hash"]:
                        self.term_names[term.ontology.id][item] = ontology["term_hash"][
                            item
                        ]["name"]
        if term.id not in self.term_names[term.ontology.id]:
            return "Unknown"
        return self.term_names[term.ontology.id][term.id]

    def get_gene_term_hash(
        self,
        prioritized_event_list=None,
        ontologies=None,
        merge_all=False,
        feature_type=None,
        translate_to_rast=True,
    ):
        output = {}
        feature_hash = self.genes
        if len(self.genes) == 0 or (feature_type == "cds" and len(self.cdss) > 0):
            feature_hash = self.cdss
        for feature_id in feature_hash:
            if not feature_type or feature_type == self.feature_types[feature_id]:
                feature = feature_hash[feature_id]
                if feature not in output:
                    output[feature] = {}
                output[feature] = feature.get_associated_terms(
                    prioritized_event_list, ontologies, merge_all, translate_to_rast
                )
        return output

    def get_reaction_gene_hash(
        self,
        prioritized_event_list=None,
        ontologies=None,
        merge_all=False,
        cds_features=False,
        feature_type=None
    ):
        output = {}
        feature_hash = self.genes
        if len(self.genes) == 0 or (cds_features and len(self.cdss) == 0):
            feature_hash = self.cdss
        for feature_id in feature_hash:
            if not feature_type or feature_type == self.feature_types[feature_id]:
                reactions = feature_hash[feature_id].get_associated_reactions(
                    prioritized_event_list, ontologies, merge_all
                )
                for rxn_id in reactions:
                    if rxn_id not in output:
                        output[rxn_id] = {}
                    if feature_id not in output[rxn_id]:
                        output[rxn_id][feature_id] = {"probability": 0, "evidence": []}
                    for item in reactions[rxn_id]:
                        output[rxn_id][feature_id]["evidence"].append(item)
        for rxn_id in output:
            total_prob = 0
            for feature_id in output[rxn_id]:
                sub_total_prob = 0
                for evidence in output[rxn_id][feature_id]["evidence"]:
                    sub_total_prob += evidence["probability"]
                output[rxn_id][feature_id]["probability"] = sub_total_prob
                total_prob += sub_total_prob
            for feature_id in output[rxn_id]:
                output[rxn_id][feature_id]["probability"] = (
                    output[rxn_id][feature_id]["probability"] / total_prob
                )
        return output

    def add_term(self, term_or_id, ontology=None):
        if not isinstance(term_or_id, AnnotationOntologyTerm):
            if term_or_id in self.terms:
                return self.terms[term_or_id]
            else:
                return AnnotationOntologyTerm(self, term_or_id, ontology)
        if term_or_id.id in self.terms:
            logger.critical("Term with id " + term_or_id.id + " already in annotation!")
            return self.terms[term_or_id.id]
        else:
            self.terms[term_or_id.id] = term_or_id

    def add_ontology(self, ontology_or_id):
        if not isinstance(ontology_or_id, AnnotationOntologyOntology):
            if ontology_or_id in self.ontologies:
                return self.ontologies[ontology_or_id]
            else:
                return AnnotationOntologyOntology(self, ontology_or_id)
        if ontology_or_id.id in self.ontologies:
            logger.critical(
                "Ontology with id " + ontology_or_id.id + " already in annotation!"
            )
            return self.ontologies[ontology_or_id.id]
        else:
            self.ontologies[ontology_or_id.id] = ontology_or_id

    def get_feature_hash(self, feature_id):
        feature_hash = self.genes
        if feature_id in self.feature_types:
            if self.feature_types[feature_id] == "cds":
                feature_hash = self.cdss
            elif self.feature_types[feature_id] == "noncoding":
                feature_hash = self.noncodings
        return feature_hash

    def add_feature(self, feature_or_id):
        feature_hash = None
        if not isinstance(feature_or_id, AnnotationOntologyFeature):
            feature_hash = self.get_feature_hash(feature_or_id)
            if feature_or_id in feature_hash:
                return feature_hash[feature_or_id]
            else:
                feature_or_id = AnnotationOntologyFeature(self, feature_or_id)
        if not feature_hash:
            feature_hash = self.get_feature_hash(feature_or_id.id)
        if feature_or_id.id not in feature_hash:
            feature_hash[feature_or_id.id] = feature_or_id
        return feature_hash[feature_or_id.id]

    def get_msgenome(self,prioritized_event_list=None,ontologies=None,merge_all=False,feature_type=None,translate_to_rast=True):
        newgenome = MSGenome.from_annotation_ontology(
            self, prioritized_event_list, ontologies, merge_all,feature_type, translate_to_rast
        )
        newgenome.annoont = self
        return newgenome
    
    def get_events_from_priority_list(self,priority_list):
        event_list = []
        for item in priority_list:
            selected_merge = None
            for event in self.events:
                if item == "all":
                    if event.id not in event_list:
                        event_list.append(event.id)
                elif item == "RAST":
                    if len(event.method) > 4 and event.method[0:4] == "RAST" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "Prokka":
                    if len(event.method) > 6 and event.method[0:6] == "Prokka" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "DRAM":
                    if len(event.method) > 4 and event.method[0:4] == "DRAM" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "GLM4EC":
                    if len(event.method) > 6 and event.method[0:6] == "GLM4EC" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "PDB":
                    if event.method == "KBAnnotationApps.PDBAnnotation" and event.ontology.id == "EC" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "SNEKMER":
                    if len(event.method) > 7 and event.method[0:7] == "Snekmer" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "Import":
                    if len(event.method) > 6 and event.method[0:6] == "Import" and event.id not in event_list:
                        event_list.append(event.id)
                elif item == "Merge":
                    if len(event.method) > 5 and event.method[0:5] == "Merge" and event.id not in event_list:
                        selected_merge = event.id
                elif item.lower() in event.description.lower() or item.lower() in event.id.lower():
                    event_list.append(event.id)
            if selected_merge:
                event_list.append(selected_merge)
        return event_list