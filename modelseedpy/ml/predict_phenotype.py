# -*- coding: utf-8 -*-
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def get_functional_roles(genome, ontology_term):
    roles = set()
    for feature in genome.features:
        if ontology_term in feature.ontology_terms:
            roles.update(feature.ontology_terms[ontology_term])
    return roles


def get_list_functional_roles_from_kbase(genome_ref, ws_client):
    list_functional_roles = []
    genome_object_data = ws_client.get_objects2({"objects": [{"ref": genome_ref}]})[
        "data"
    ][0]["data"]

    # determine where functional roles are kept
    keys_location = genome_object_data.keys()
    if "features" in keys_location:
        location_of_functional_roles = genome_object_data["features"]
    elif "non_coding_features" in keys_location:
        location_of_functional_roles = genome_object_data["non_coding_features"]
    elif "cdss" in keys_location:
        location_of_functional_roles = genome_object_data["cdss"]
    else:
        raise ValueError(
            "The functional roles are not under 'features', 'non_coding_features', or 'cdss'."
        )

    # either the functional roles are under function or functions
    keys_function = location_of_functional_roles[0].keys()
    function_str = "function" if "function" in keys_function else "functions"
    for functional_role in location_of_functional_roles:
        try:
            role_to_insert = functional_role[function_str][0]
            if " @ " in role_to_insert:
                list_functional_roles.extend(role_to_insert.split(" @ "))
            elif " / " in role_to_insert:
                list_functional_roles.extend(role_to_insert.split(" / "))
            elif "; " in role_to_insert:
                list_functional_roles.extend(role_to_insert.split("; "))
            elif "hypothetical protein" in role_to_insert:
                pass
            else:
                list_functional_roles.append(role_to_insert)
        except KeyError as e:
            logger.error(e)
            # print("this is funcitonal role")
            # print(functional_role)
            # print("this is list_functional_roles")
            # print(list_functional_roles)

            # print("apparently some function list just don't have functions...")
            # ^^ this makes no sense...
            pass

    return list_functional_roles


def create_indicator_matrix_from_genomes(genomes, ontology_term, master_role_list=None):
    ref_to_role = {}
    for genome in genomes:
        ref_to_role[genome.id] = get_functional_roles(genome, ontology_term)

    return create_indicator_matrix(ref_to_role, master_role_list)


def _create_sorted_master_role_list(ref_to_role):
    master_role_set = set()
    for i, feature in ref_to_role.items():
        master_role_set.update(feature)

    master_role_list = sorted(list(master_role_set))
    return master_role_list


def create_indicator_matrix(ref_to_role, master_role_list=None):
    if master_role_list is None:
        master_role_list = _create_sorted_master_role_list(ref_to_role)

    ref_to_indication = {}
    # make indicator rows for each
    for genome_id, features in ref_to_role.items():
        matching_index = [
            i for i, role in enumerate(master_role_list) if role in set(features)
        ]
        indicators = np.zeros(len(master_role_list))
        try:
            indicators[np.array(matching_index)] = 1
        except IndexError:
            raise IndexError(
                "The genomes or genomeSet that you have submitted wasn’t annotated using the \
                RAST annotation pipeline. Please annotate the genomes via ‘Annotate Microbial Genome’ app \
                (https://narrative.kbase.us/#appcatalog/app/RAST_SDK/reannotate_microbial_genome/release)or \
                genomeSets via Annotate Multiple Microbial Genomes’ app \
                (https://narrative.kbase.us/#appcatalog/app/RAST_SDK/reannotate_microbial_genomes/release) and \
                resubmit the RAST annotated genome/genomeSets into the Predict Phenotype app. ("
            )
        ref_to_indication[genome_id] = indicators.astype(int)

    indicator_matrix = (
        pd.DataFrame.from_dict(
            data=ref_to_indication, orient="index", columns=master_role_list
        )
        .reset_index()
        .rename(columns={"index": "Genome Reference"})
    )
    return indicator_matrix, master_role_list
