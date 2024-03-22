# -*- coding: utf-8 -*-

from __future__ import absolute_import
from collections import deque, defaultdict
import os
import glob
import sys
import warnings
from itertools import combinations
import re
import pandas as pd
import numpy as np
import cobra
import networkx as nx

from modelseedpy.community import commhelper
from modelseedpy import MSModelUtil

warnings.filterwarnings("ignore")


def _create_graph_with_internal_reaction(organismsdata):
    """
    This function creates a NetworkX DiGraph object which consists of
    reactions and metabolites happening inside the organisms in a community.
    This makes use of the reaction information i.e., irreversible and
    reversible, which is obtained from another script fetch_reactions.

    Parameters
    ----------
    organismsdata : dict
        Dictionary containing the reaction information about organisms

    Returns
    -------
    G : NetworkX DiGraph Object
        Bipartite graph consisting of internal reactions in organisms
    """
    G = nx.DiGraph()
    for modelname in organismsdata:
        G.add_nodes_from(organismsdata[modelname]['irreversible_rxn_no'], bipartite=1)
        G.add_nodes_from(organismsdata[modelname]['reversible_rxn_no'], bipartite=1)
        G.add_nodes_from(organismsdata[modelname]['reversible_back_rxn_no'], bipartite=1)
        irrev_lhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['irreversible_lhs_nodes'] for item in sublist]))
        irrev_rhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['irreversible_rhs_nodes'] for item in sublist]))
        rev_lhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['reversible_lhs_nodes'] for item in sublist]))
        rev_rhs_nodes = list(set(
            [item for sublist in organismsdata[modelname]['reversible_rhs_nodes'] for item in sublist]))
        G.add_nodes_from(irrev_lhs_nodes, bipartite=0)
        G.add_nodes_from(irrev_rhs_nodes, bipartite=0)
        G.add_nodes_from(rev_lhs_nodes, bipartite=0)
        G.add_nodes_from(rev_rhs_nodes, bipartite=0)
        for irrevidx in range(len(organismsdata[modelname]['irreversible_rxn_no'])):
            for lhsmetidx in range(len(organismsdata[modelname]['irreversible_lhs_nodes'][irrevidx])):
                G.add_edges_from([(organismsdata[modelname]['irreversible_lhs_nodes'][irrevidx][lhsmetidx],
                                   organismsdata[modelname]['irreversible_rxn_no'][irrevidx])])
            for rhsmetidx in range(len(organismsdata[modelname]['irreversible_rhs_nodes'][irrevidx])):
                G.add_edges_from([(organismsdata[modelname]['irreversible_rxn_no'][irrevidx],
                                   organismsdata[modelname]['irreversible_rhs_nodes'][irrevidx][rhsmetidx])])
        for revidx in range(len(organismsdata[modelname]['reversible_rxn_no'])):
            for lhsmetidxrev in range(len(organismsdata[modelname]['reversible_lhs_nodes'][revidx])):
                G.add_edges_from([(organismsdata[modelname]['reversible_lhs_nodes'][revidx][lhsmetidxrev],
                                   organismsdata[modelname]['reversible_rxn_no'][revidx])])
                G.add_edges_from([(organismsdata[modelname]['reversible_back_rxn_no'][revidx],
                                   organismsdata[modelname]['reversible_lhs_nodes'][revidx][lhsmetidxrev])])
            for rhsmetidxrev in range(len(organismsdata[modelname]['reversible_rhs_nodes'][revidx])):
                G.add_edges_from([(organismsdata[modelname]['reversible_rxn_no'][revidx],
                                   organismsdata[modelname]['reversible_rhs_nodes'][revidx][rhsmetidxrev])])
                G.add_edges_from([(organismsdata[modelname]['reversible_rhs_nodes'][revidx][rhsmetidxrev],
                                   organismsdata[modelname]['reversible_back_rxn_no'][revidx])])
    return G


def _create_graph_with_exchange_reactions(G, orgs, namemap):
    """
    This function first identifies the common exchange metabolites
    and the non-common exchange metabolites and adds them to the
    DiGraph object generated above.

    Parameters
    ----------
    G : NetworkX DiGraph Object
        Bipartite graph of reaction network from organisms
    orgs : dict
        Dictionary consisting of irreversible, reversible and exchange
        reactions pertaining to the organisms. If more than one organism
        is used, this dictionary consists of information about all the
        organisms.
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model

    Returns
    -------
    G : NetworkX DiGraph Object
        Bipartite graph consisting of internal and exchange reactions in organisms
    namemap : dict
        Dictionary mapping the adhoc exchange reaction names to reaction names in
        the model
    """
    metabolite_exchanged = []
    for orgnames in orgs:
        exc_met = orgs[orgnames]['exchange_metab_nodes']
        metabolite_exchanged.append(exc_met)
    # Common exchange metabolites in different organisms
    common_exchange_metabolite = list(set.intersection(*list(map(set, metabolite_exchanged))))
    common_exchange_metabolite.sort()
    #  Adding the common exchange metabolites to the graph
    for orgnames in orgs:
        renamed_exc_met = [f"{orgnames} {comexcmet}" for comexcmet in common_exchange_metabolite]
        number_exc_met = list(range(0, len(common_exchange_metabolite)))
        mod_exc_rxn_number = [f'Org_{orgnames} ER{str(num + 1)}' for num in number_exc_met]
        mod_exc_rev_rxn_number = [f'Org_{orgnames} ERR{str(num + 1)}' for num in number_exc_met]
        G.add_nodes_from(mod_exc_rxn_number, bipartite=1)
        G.add_nodes_from(mod_exc_rev_rxn_number, bipartite=1)
        G.add_nodes_from(common_exchange_metabolite, bipartite=0)
        G.add_nodes_from(renamed_exc_met, bipartite=0)
        for k in range(len(renamed_exc_met)):
            namemap[mod_exc_rxn_number[k]] = common_exchange_metabolite[k]
            namemap[mod_exc_rev_rxn_number[k]] = common_exchange_metabolite[k]
            G.add_edges_from([(renamed_exc_met[k], mod_exc_rxn_number[k])])
            G.add_edges_from([(mod_exc_rxn_number[k], common_exchange_metabolite[k])])
            G.add_edges_from([(common_exchange_metabolite[k], mod_exc_rev_rxn_number[k])])
            G.add_edges_from([(mod_exc_rev_rxn_number[k], renamed_exc_met[k])])
    #  Adding the uncommon exchange metabolites to the graph
    for orgnames in orgs:
        metitems = orgs[orgnames]['exchange_metab_nodes']
        non_common_exc_met = list(set(metitems) - set(common_exchange_metabolite))
        non_common_exc_met.sort()
        renamed_non_common_exc_met = [f"{orgnames} {s}" for s in non_common_exc_met]
        number_non_common_exc_met = list(range(0, len(non_common_exc_met)))
        mod_non_common_exc_rxn_number = [f"Org_{orgnames} NCER{str(num + 1)}" for num in number_non_common_exc_met]
        mod_non_common_exc_rev_rxn_number = [f"Org_{orgnames} NCERR{str(num + 1)}"
                                             for num in number_non_common_exc_met]
        G.add_nodes_from(mod_non_common_exc_rxn_number, bipartite=1)
        G.add_nodes_from(mod_non_common_exc_rev_rxn_number, bipartite=1)
        G.add_nodes_from(non_common_exc_met, bipartite=0)
        G.add_nodes_from(renamed_non_common_exc_met, bipartite=0)
        for k in range(len(renamed_non_common_exc_met)):
            namemap[mod_non_common_exc_rxn_number[k]] = non_common_exc_met[k]
            namemap[mod_non_common_exc_rev_rxn_number[k]] = non_common_exc_met[k]
            G.add_edges_from([(renamed_non_common_exc_met[k], mod_non_common_exc_rxn_number[k])])
            G.add_edges_from([(mod_non_common_exc_rxn_number[k], non_common_exc_met[k])])
            G.add_edges_from([(non_common_exc_met[k], mod_non_common_exc_rev_rxn_number[k])])
            G.add_edges_from([(mod_non_common_exc_rev_rxn_number[k], renamed_non_common_exc_met[k])])
    return G, namemap

def create_graph(file_names, no_of_orgs):
    """
    This function creates bipartite graph of the organisms based on the
    path provided and the number of organsisms. For instance, if a folder
    has 3 model files, and the number of organisms is 2, 3 (3C2) different
    bipartite graphs are created. The graph objects and the dictionary
    are saved as gpickle and pickle files respectively.

    Parameters
    ----------
    file_names : list
        List containing the file names of models
    no_of_orgs : int
        Number of organisms to be used for creating the DiGraph.

    Returns
    -------
    H : NetworkX DiGraph Object
        Bipartite graph consisting of internal and exchange reactions in organisms
    full_name_map : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model
    """

    H=[]
    organisms_reaction_data, partial_name_map = segregate_reactions_from_models(file_names)
    if organisms_reaction_data:
        organisms_names = list(organisms_reaction_data.keys())
        all_possible_combis = list(combinations(list(range(len(organisms_names))), int(no_of_orgs)))
        if int(no_of_orgs)>1 and sorted(organisms_names)[0][0]=='0':
            all_possible_combis = all_possible_combis[:len(organisms_names)-1]
        if all_possible_combis:
            for ncom in range(len(all_possible_combis)):
                file_name = ''
                current_combination = {}
                for numincom in range(len(all_possible_combis[ncom])):
                    current_combination[organisms_names[all_possible_combis[ncom][numincom]]] = \
                        organisms_reaction_data[organisms_names[all_possible_combis[ncom][numincom]]]
                    file_name = file_name + organisms_names[all_possible_combis[ncom][numincom]] + '_'
                H.append(_create_graph_with_internal_reaction(current_combination))
                temp, full_name_map = _create_graph_with_exchange_reactions(
                    H[ncom], current_combination, partial_name_map)
                H[ncom]=temp
                print(len(H), H[ncom])
                print('Number of edges in graph', len(H[ncom].edges()))
                print('Number of nodes in graph', len(H[ncom].nodes()))

                # Uncomment the following code to save the graph files externally in your machine
                # Note: Graph files can occupy a large space for large datasets
                '''
                if os.access(path_name_with_models, os.W_OK):
                    with open(file_name + 'namemap' + '.pickle', 'wb') as filetodump:
                        dump(full_name_map, filetodump)
                    nx.write_gpickle(H[ncom], file_name + '.gpickle')
                    print('Graph and namemap saved for file(s) in', path_name_with_models)
                '''
        else:
            print(
                'Number of organisms for creating a consortium graph is more than the models given')
            print('Program will now exit')
            sys.exit()
    else:
        print("Cannot create graph")
        sys.exit()
    return H, full_name_map


def forward_pass(graph_object, media):
    """
    This function carries out the Guided Breadth First Search on a directed
    bipartite graph starting from the entries in seed metabolite set.

    Parameters
    ----------
    graph_object : NetworkX DiGraph Object
        Bipartite graph of the metabolic network

    seedmet : set
        Set of seed metabolites including the source

    Returns
    -------
    lower_bound_metabolite : defaultdict
        Minimum number of steps required to reach a metabolite
    status_dict : defaultdict
        Dictionary pertaining to the status of every reaction - whether it
        has been visited or not
    scope : set
        Set of metabolites that can be produced from the given set of
        seed metabolites

    Notes
    -----
    Starting with the set of seed metabolites S, the algorithm first finds
    all the reactions from the set R, whose precursor metabolites are in S.
    Such reactions are marked visited and added to the visited reaction set.
    Metabolites produced by these reactions are checked. The reactions where
    these metabolites participate are then checked for the presence of all its
    predecessors and are added to the queue. This traversal continues in a
    breadth-first manner and stops when there are no further reactions to
    be visited.
    """
    pred = graph_object.predecessors
    succ = graph_object.successors
    lower_bound_metabolite = {cpd: [0] for cpd in media}
    lower_bound_reaction = defaultdict(list)
    status_dict = defaultdict(str)
    # Using a deque since deques have O(1) speed for appendleft() and popleft()
    # while lists have O(n) performance for inserting and popping.
    queue = deque([])
    # All seed metabolites are always present, hence require 0 steps
    stage = 1
    mediaMets = list(media.keys())
    scope = list(media.keys())
    starting_rxn_node = []
    # First stage where starting_rxn_node list contains all the reactions
    # which require only the seed metabolites as input
    for starting_met_nodes in mediaMets:
        # Essential when analysing mutiple networks with same seed metabolite
        # set, although would be redundant in case of single network
        if starting_met_nodes in graph_object:
            for startingrxns in succ(starting_met_nodes):
                if set(pred(startingrxns)).issubset(mediaMets):
                    if startingrxns not in starting_rxn_node:
                        starting_rxn_node.append(startingrxns)
                    for metsprod in succ(startingrxns):
                        scope.add(metsprod)
                        if stage not in lower_bound_metabolite[metsprod]:
                            lower_bound_metabolite[metsprod].append(stage)
                    if stage not in lower_bound_reaction[startingrxns]:
                        lower_bound_reaction[startingrxns].append(stage)
    for rxn in starting_rxn_node:
        for metabs in succ(rxn):
            for nextrxn in succ(metabs):
                if set(pred(nextrxn)).issubset(scope):
                    if nextrxn not in queue:  queue.append(nextrxn)
        status_dict[rxn] = 'V'
    while queue:
        stage += 1
        for parentrxn in list(queue):
            if status_dict[parentrxn] == '':
                if stage not in lower_bound_reaction[parentrxn]:
                    lower_bound_reaction[parentrxn].append(stage)
                for mets in succ(parentrxn):
                    scope.add(mets)
                    if stage not in lower_bound_metabolite[mets]:
                        lower_bound_metabolite[mets].append(stage)
                    for progeny in succ(mets):
                        if set(pred(progeny)).issubset(scope):
                            if status_dict[progeny] != 'V':
                                if progeny not in queue:  queue.append(progeny)
                status_dict[parentrxn] = 'V'
            elif status_dict[parentrxn] == 'V':
                for mets in succ(parentrxn):
                    if stage not in lower_bound_metabolite[mets]:  lower_bound_metabolite[mets].append(stage)
            queue.popleft()
    return lower_bound_metabolite, status_dict, scope

def find_different_reaction_types(stoi_matrix, model, current_model_name):
    """
    This function finds the exchange, irreversible and the reversible reactions
    from the model.

    Parameters
    ----------
    stoi_matrix : numpy array
        full path name where the model files are
    model : COBRA model object
        COBRA model object created from SBML models
    current_model_name : str
        Name which is to be prefixed against every
        reaction/metabolite (to differentiate the entries in multiple organisms,
        when a community model is built)
    Returns
    -------
    exchange_met_ids : list
        Metabolite identifiers of exchange metabolites
    irrev_lhs_nodes : list
        Metabolite identifiers of reactants of irreversible reactions
    irrev_rhs_nodes : list
        Metabolite identifiers of products of irreversible reactions
    rev_lhs_nodes : list
        Metabolite identifiers of reactants of reversible reactions
    rev_rhs_nodes : list
        Metabolite identifiers of products of reversible reactions
    exchange_rxn_ids : list
        Reaction identifers of exchange reactions
    irrev_rxn_ids : list
        Reaction identifiers of irreversible reactions
    rev_rxn_ids : list
        Reaction identifiers of reversible reactions

    """

    xdim = np.shape(stoi_matrix)
    reactants_of_reaction, total_metabolites_in_reaction, products_of_reaction = [], [], []
    number_of_reactants_in_reaction, total_number_of_metabs_in_reaction = [], []
    number_of_products_in_reaction, exchange_reaction_idx = [], []
    reaction_identifiers, reaction_in_model, metabolite_identifiers = [], [], []
    for metab in model.metabolites:
        metabolite_identifiers.append(metab.id)
    for rxns in model.reactions:
        reaction_identifiers.append(rxns.id)
        reaction_in_model.append(rxns.reaction)
    for rxnidx in range(xdim[0]):
        reactants_of_reaction.append(np.where(stoi_matrix[rxnidx] == -1))
        total_metabolites_in_reaction.append(np.where(stoi_matrix[rxnidx] != 0))
        products_of_reaction.append(np.where(stoi_matrix[rxnidx] == 1))
        number_of_reactants_in_reaction.append(len(reactants_of_reaction[rxnidx][0]))
        total_number_of_metabs_in_reaction.append(len(total_metabolites_in_reaction[rxnidx][0]))
        number_of_products_in_reaction.append(len(products_of_reaction[rxnidx][0]))

        # Case 1 - Presence of bulk metabolites in the medium

        if reaction_in_model[rxnidx][-1] == 'b':  # Assuming the bulk metabolites end in 'b'
            if number_of_reactants_in_reaction[rxnidx] == 1 and number_of_products_in_reaction[rxnidx] == 1:
                exchange_reaction_idx.append(rxnidx)
        # Case 2 - Presence of exchange metabolites
        elif number_of_reactants_in_reaction[rxnidx] == 1 and total_number_of_metabs_in_reaction[rxnidx] == 1:
            exchange_reaction_idx.append(rxnidx)
        elif number_of_products_in_reaction[rxnidx] == 1 and total_number_of_metabs_in_reaction[rxnidx] == 1:
            exchange_reaction_idx.append(rxnidx)
    exchange_met_ids, exchange_met_index, exchange_rxn_ids = [], [], []
    for excentry in exchange_reaction_idx:
        exchange_rxn_ids.append(reaction_identifiers[excentry])
        if reaction_in_model[excentry][-1] == 'b':
            exchange_met_ids.append(metabolite_identifiers[np.nonzero(stoi_matrix[excentry])[0][0]])
        else:  exchange_met_index.append(np.nonzero(stoi_matrix[excentry])[0].tolist()[0])
    if exchange_met_index:
        for metind in exchange_met_index:
            exchange_met_ids.append(metabolite_identifiers[metind])
    all_rxn_idx = list(range(len(reaction_in_model)))
    internal_rxns = list(set(all_rxn_idx) ^ set(exchange_reaction_idx))
    reversible_rxns, irreversible_rxns, rxns_lowerbound, rxns_upperbound = [], [], [], []
    for rxns in model.reactions:
        rxns_lowerbound.append(rxns.lower_bound)  ;  rxns_upperbound.append(rxns.upper_bound)
    for idxint in internal_rxns:
        if rxns_lowerbound[idxint] < 0 and rxns_upperbound[idxint] >= 0:  reversible_rxns.append(idxint)
        elif rxns_lowerbound[idxint] >= 0 and rxns_upperbound[idxint] >= 0:  irreversible_rxns.append(idxint)
    #  Irreversible reaction nodes
    irrev_lhs_temporary, irrev_rhs_temporary, irrev_lhs_nodes, irrev_rhs_nodes, irrev_rxn_ids = [], [], [], [], []
    for irridx in irreversible_rxns:
        irrev_rxn_ids.append(reaction_identifiers[irridx])
        irrev_lhs_temporary.append(np.where(stoi_matrix[irridx] < 0)[0].tolist())
        irrev_rhs_temporary.append(np.where(stoi_matrix[irridx] > 0)[0].tolist())
    for lhsirridx in range(len(irrev_lhs_temporary)):
        temp_metab_list_lhs = []
        for met_idx_lhs in irrev_lhs_temporary[lhsirridx]:
            met_namech_lhs = f"{current_model_name} {metabolite_identifiers[met_idx_lhs]}"
            temp_metab_list_lhs.append(met_namech_lhs)
        irrev_lhs_nodes.append(temp_metab_list_lhs)
    for rhsirridx in range(len(irrev_rhs_temporary)):
        temp_metab_list_rhs = []
        for met_idx_rhs in irrev_rhs_temporary[rhsirridx]:
            met_namech_rhs = f"{current_model_name} {metabolite_identifiers[met_idx_rhs]}"
            temp_metab_list_rhs.append(met_namech_rhs)
        irrev_rhs_nodes.append(temp_metab_list_rhs)

    #  Reversible reaction nodes
    rev_lhs_temporary, rev_rhs_temporary, rev_lhs_nodes, rev_rhs_nodes, rev_rxn_ids = [], [], [], [], []
    for rridx in reversible_rxns:
        rev_rxn_ids.append(reaction_identifiers[rridx])
        rev_lhs_temporary.append(np.where(stoi_matrix[rridx] < 0)[0].tolist())
        rev_rhs_temporary.append(np.where(stoi_matrix[rridx] > 0)[0].tolist())
    for lhsrevidx in range(len(rev_lhs_temporary)):
        temp_metab_list_lhs_rev = []
        for met_idx_lhs in rev_lhs_temporary[lhsrevidx]:
            met_namech_lhs = "%s %s" % (current_model_name, metabolite_identifiers[met_idx_lhs])
            temp_metab_list_lhs_rev.append(met_namech_lhs)
        rev_lhs_nodes.append(temp_metab_list_lhs_rev)
    for rhsrevidx in range(len(rev_rhs_temporary)):
        temp_metab_list_rhs_rev = []
        for met_idx_rhs in rev_rhs_temporary[rhsrevidx]:
            met_namech_rhs = "%s %s" % (current_model_name, metabolite_identifiers[met_idx_rhs])
            temp_metab_list_rhs_rev.append(met_namech_rhs)
        rev_rhs_nodes.append(temp_metab_list_rhs_rev)
    return (exchange_met_ids, irrev_lhs_nodes, irrev_rhs_nodes, rev_lhs_nodes,
            rev_rhs_nodes, exchange_rxn_ids, irrev_rxn_ids, rev_rxn_ids)


def segregate_reactions_from_models(models):
    """
    This function gets the data pertaining to the reactions and the
    metabolites from the models of multiple organisms.
    This requires as input the pathname where the '.xml' files are located.
    From this path, this function reads all the files using the functions
    in the COBRA toolbox and generates the stoichiometric model for these
    SBML models.

    Parameters
    ----------
    models : list
        List of model objects

    Returns
    -------
    all_organisms_info : dict
        Dictionary of all model data (reaction information about all the
        organisms)
    namemap : dict
        Dictionary mapping the adhoc reaction names to reaction names in
        the model

    """
    all_organisms_info = {}
    namemap = {}
    for model in models:
        stoi = cobra.util.array.create_stoichiometric_matrix(model)
        current_organisms_info = {}
        rxns_in_model, mets_in_model = [], []
        for metab in model.metabolites:
            mets_in_model.append(metab.id)
        for reac in model.reactions:
            rxns_in_model.append(reac.id)
        stoi_matrix = stoi.T
        (exchange_nodes, irrev_lhs_nodes, irrev_rhs_nodes, rev_lhs_nodes, rev_rhs_nodes,
         exc_name, irrev_rxn_name, rev_rxn_name
         ) = find_different_reaction_types(stoi_matrix, model, model.id)
        current_organisms_info[model.id] = {
            'exchange_metab_nodes': exchange_nodes, 'irreversible_lhs_nodes': irrev_lhs_nodes,
            'irreversible_rhs_nodes': irrev_rhs_nodes, 'reversible_lhs_nodes': rev_lhs_nodes,
            'reversible_rhs_nodes': rev_rhs_nodes, 'exch_rxn_name': exc_name, 'irrev_rxn_name': irrev_rxn_name,
            'rev_rxn_name': rev_rxn_name}

        irrev_rxn_number = []
        for num in range(len(irrev_lhs_nodes)):
            modified_name_irrev = f'Org_{model.id} IR' + str(num + 1)
            irrev_rxn_number.append(modified_name_irrev)
            namemap[modified_name_irrev] = irrev_rxn_name[num]

        rev_rxn_number = []
        for num in range(len(rev_lhs_nodes)):
            modified_name_rev = f'Org_{model.id} RR' + str(num + 1)
            rev_rxn_number.append(modified_name_rev)
            namemap[modified_name_rev] = rev_rxn_name[num]

        rev_back_rxn_number = []
        for num in range(len(rev_lhs_nodes)):
            modified_name_back_rev = f'Org_{model.id} RevBR' + str(num + 1)
            rev_back_rxn_number.append(modified_name_back_rev)
            namemap[modified_name_back_rev] = rev_rxn_name[num]

        current_organisms_info[model.id]['reversible_rxn_no'] = rev_rxn_number
        current_organisms_info[model.id]['irreversible_rxn_no'] = irrev_rxn_number
        current_organisms_info[model.id]['total_nodes'] = len(
            exchange_nodes) + len(irrev_lhs_nodes) + len(rev_lhs_nodes)
        current_organisms_info[model.id]['model_rxns'] = rxns_in_model
        current_organisms_info[model.id]['reversible_back_rxn_no'] = rev_back_rxn_number
        current_organisms_info[model.id]['metabolites'] = mets_in_model
        all_organisms_info.update(current_organisms_info)
    return all_organisms_info, namemap

def find_relievedrxns(model, org_info, org_info_pert):
    relieved = {i: list(set(org_info_pert[i]) - set(org_info[i])) for i in org_info_pert}
    detailed_rel_rxns, rel_rxns_name = {}, {}

    for i in model:
        j = i.id
        detailed_rel_rxns[j] = []
        rel_rxns_name[j] = []
        if len(relieved[j]):
            rxn_ids = []
            for r in i.reactions:
                rxn_ids.append(r.id)
            for rel in relieved[j]:
                rel_rxn = i.reactions[rxn_ids.index(rel)].reaction
                detailed_rel_rxns[j].append(rel_rxn)
                rel_rxns_name[j].append(i.reactions[rxn_ids.index(rel)].name)

    return relieved, detailed_rel_rxns, rel_rxns_name

def find_stuckrxns(model, community, media, no_of_orgs):
    # Constructing graphs
    warnings.filterwarnings("ignore")
    G, full_name_map = create_graph(community, no_of_orgs)
    if not os.path.exists('results'):  os.makedirs('results')
    all_possible_combis = list(combinations(list(range(len(community))), int(no_of_orgs)))
    if no_of_orgs > 1 and sorted(community)[0][0] == '0':
        all_possible_combis = all_possible_combis[:len(community) - 1]
    org_info = {}
    scope = {}
    print('No. of graphs constructed: ', len(G))

    # This loop finds all the stuck reaction
    for i in range(len(all_possible_combis)):
        lbm, sd, s = forward_pass(G[i], media)
        for j in range(len(all_possible_combis[i])):
            stuck, rxnNode = [], []
            model1 = model[all_possible_combis[i][j]].id
            visited = list(sd.keys())
            for r in G[i].nodes:
                if r.find(model1) >= 0:  rxnNode.append(r)
            for rxn in rxnNode:
                if rxn in visited:  continue
                elif rxn.find('ERR') >= 0:  continue
                elif rxn.find('Org') >= 0:
                    if (rxn[len(model1) + 5] == 'I') or (rxn[len(model1) + 5] == 'R'):  stuck.append(rxn)
            org_info[model1] = stuck
            scope[model1] = s
    return org_info, scope, full_name_map

def decrypt_orginfo(org_info, namemap):
    """
    This function decrypts the rxn ids using the data in corresponding namemaps
    :param org_info:
    :param namemap:
    :return:
        org_info: An dictionary of decrypted rxn ids for each community
    """
    for i in org_info:
        for j in range(len(org_info[i])):
            org_info[i][j] = namemap[org_info[i][j]]
    return org_info

def make_perturbed_community(rem_org, pert_models, pert_community):
    pert_model_ids = [i.id for i in pert_models]
    for i in rem_org:
        if i in pert_model_ids:
            pert_models.remove(pert_models[pert_model_ids.index(i)])
            pert_community.remove(pert_community[pert_model_ids.index(i)])
            pert_model_ids.remove(i)

    return pert_models, pert_community, pert_model_ids

def perform_task(media, model, transport_rxns, pert_community,
                 org_info_wo_trans_rxn, rem_org_list, n):
    org_info_pert, scope_pert, namemap_pert = find_stuckrxns(model, pert_community, media, len(pert_community))
    org_info_pert = decrypt_orginfo(org_info_pert, namemap_pert)
    org_info_pert_wo_trans_rxn = {i:list(set(org_info_pert[i]) - set(transport_rxns)) for i in org_info_pert}

    with open(f"results/Community_without_clus{str(n)}.csv", "w") as g:
        for m in org_info_pert_wo_trans_rxn:
            g.write(m + ',' + str(len(org_info_pert_wo_trans_rxn[m])) + '\n')
    stuck_com = stuck_pert_com = 0
    for i in org_info_wo_trans_rxn:
        if i not in rem_org_list:  stuck_com += len(org_info_wo_trans_rxn[i])
    for i in org_info_pert_wo_trans_rxn:
        stuck_pert_com += len(org_info_pert_wo_trans_rxn[i])
    msi = 1 - (stuck_com / stuck_pert_com)
    print(n, 'th cluster')
    return org_info_pert, org_info_pert_wo_trans_rxn, msi

def write_relieved_rxns(g, relieved, detailed_rel_rxns, rel_rxns_name):
    g.write('acceptor\trelieved reactions\n')
    for i in relieved:
        g.write(i + '\t')
        for j in list(set(relieved[i])):
            g.write(j + '\t\n\t')
        for d in list(set(rel_rxns_name[i])):
            g.write(d + '\t\n\t')
        for k in list(set(detailed_rel_rxns[i])):
            g.write(k + '\t\n')

def write_relieved_rxn_metadata(h, org_info_wo_trans_rxn, org_info_pert_wo_trans_rxn):
    nrelieved = {}
    for i in org_info_pert_wo_trans_rxn:
        nrelieved[i] = len(org_info_pert_wo_trans_rxn[i]) - len(org_info_wo_trans_rxn[i])
        if nrelieved[i]:
            h.write(i + ',' + str(len(org_info_wo_trans_rxn[i])) + ',' + str(
                len(org_info_pert_wo_trans_rxn[i])) + ',' + str(nrelieved[i]) + '\n')

def find_relieved_rxn(model, media_name, org_info_single, org_info_pair):
    """
    This function extracts and writes the relieved rxns into a tsv file
    :param model:
    :param media_name: name of the media used (identifer to know what media is used when analysis is done using multiple media)
    :param org_info_single: Dictionary containing stuck reactions of all microbes in the community
    :param org_info_pair: Dictionary containing stuck reactions of all microbes in the community
    :return: None
    """
    relieved = {}
    for org1 in model:
        for org2 in model:
            if org1.id + '_' + org2.id in org_info_pair.keys():
                relieved[org1.id + '_' + org2.id] = []
                temp = list(set(org_info_single[org1.id + '_' + org1.id]) - set(org_info_pair[org1.id + '_' + org2.id]))
                for j in temp:
                    relieved[org1.id + '_' + org2.id].append(j)
            else:  continue

    rel_rxns_name, detailed_rel_rxns = {}, {}
    for i in model:
        rxn_ids = [r.id for r in i.reactions]
        for j in model:
            org1 = i.id  ;  org2 = j.id
            if org1 + '_' + org2 in relieved.keys():
                detailed_rel_rxns[org1 + '_' + org2] = []
                rel_rxns_name[org1 + '_' + org2] = []
                for rel in relieved[org1 + '_' + org2]:
                    rel_rxn = i.reactions[rxn_ids.index(rel)].reaction
                    detailed_rel_rxns[org1 + '_' + org2].append(rel_rxn)
                    rel_rxns_name[org1 + '_' + org2].append(i.reactions[rxn_ids.index(rel)].name)

    relieved_rxn_output_file = f'results/relieved_rxns_{media_name}_w_excrxns.tsv'
    with open(relieved_rxn_output_file, 'w') as g:
        header = 'acceptor\tdonor\trelieved reactions\n'
        g.write(header)
        for i in model:
            for j in model:
                org1 = i.id  ;  org2 = j.id
                if org1 + '_' + org2 in relieved.keys():
                    g.write(org1 + '\t' + org2 + '\t')
                    rel_rxns = list(set(relieved[org1 + '_' + org2]))
                    det_rel_rxns = list(set(detailed_rel_rxns[org1 + '_' + org2]))
                    rel_rxn_nam = list(set(rel_rxns_name[org1 + '_' + org2]))
                    for x in rel_rxns:
                        g.write(x + '\t\n\t\t')
                    for d in rel_rxn_nam:
                        g.write(d + '\t\n\t\t')
                    for k in det_rel_rxns:
                        g.write(k + '\t\n')
    print('relieved reactions are written at:\n', relieved_rxn_output_file)

def find_stuck_rxns(models, community, media, comm_size):
    """
    Constructs graphs using MetQuest and finds all stuck reactions in the cellular compartment
    :param models: list of GEMs
    :param community: the community model
    :param seedmet_file: path to txt file containing seed metabolites
    :param comm_size: number of organisms in a community
    :return:
        org_info: Dictionary containing stuck reactions of all microbes in the community
        scope: Dictionary containing all the metabolites that can be produced by the microbes in the community
        namemap: Dictionaru containing all the decrypted rxn ids
    """
    warnings.filterwarnings("ignore")
    G, full_name_map = create_graph(community, comm_size)
    if not os.path.exists('results'):  os.makedirs('results')

    all_possible_combis = combinations(models, comm_size)
    org_info, scope, vis = {}, {}, {}
    print('No. of graphs constructed: ', len(G))

    # This loop finds all the stuck reaction
    for i in range(len(all_possible_combis)):
        lbm, sd, s = forward_pass(G[i], media)
        for j in range(len(all_possible_combis[i])):
            stuck, rxnNode = [], []
            model1 = models[all_possible_combis[i][j]].id
            visited = list(sd.keys())
            for r in G[i].nodes:
                if r.find(model1) >= 0:  rxnNode.append(r)
            for rxn in rxnNode:
                if rxn in visited or rxn.find('ERR') >= 0:  continue
                elif rxn.find('Org') >= 0:
                    if (rxn[len(model1) + 5] == 'I') or (rxn[len(model1) + 5] == 'R'):  stuck.append(rxn)
            model2 = models[all_possible_combis[i][j - 1]].id
            org_info[model1 + '_' + model2] = stuck
            scope[model1 + '_' + model2] = s
            vis[model1 + '_' + model2] = visited
    return org_info, scope, full_name_map, vis

def decrypt_org_info(org_info, namemap):
    """
    This function decrypts the rxn ids using the data in corresponding namemaps
    :param org_info:
    :param namemap:
    :return:
        org_info: An dictionary of decrypted rxn ids for each community
    """
    for i in org_info:
        for j in range(len(org_info[i])):
            org_info[i][j] = namemap[org_info[i][j]]
    return org_info

def pMSI(models, media):
    """
    Calculates MSI for CarveMe models
    Extracts and writes relieved reactions in every pair
    :param community: list of GSMM files
    :param sd_file: path to txt file containing seed metabolites
    :return: msi: Dictionary containing MSI values for every pair
    """
    # find all transport reactions
    community_model = commhelper.build_from_species_models(models)
    comm_util = MSModelUtil(community_model)
    # find stuck reactions
    org_info_single, scope_sin, namemap_sin, vis = find_stuck_rxns(models, community_model, media, 1)
    org_info_pair, scope_pair, namemap_pair, vis = find_stuck_rxns(models, models, media, 2)
    # decrypt the stuck reactions
    org_info_single = decrypt_org_info(org_info_single, namemap_sin)
    org_info_pair = decrypt_org_info(org_info_pair, namemap_pair)
    # Filter out the transport reactions from every stuck reaction list
    org_info_single_wo_trans_rxn, org_info_pair_wo_trans_rxn = {}, {}
    for i in org_info_single:
        org_info_single_wo_trans_rxn[i] = list(set(org_info_single[i]) - set(comm_util.transport_list()))
    for i in org_info_pair:
        org_info_pair_wo_trans_rxn[i] = list(set(org_info_pair[i]) - set(comm_util.transport_list()))
    # find all the relieved reactions in every pairs
    find_relieved_rxn(models, "relieved_rxns", org_info_single, org_info_pair)
    # calculate MSI for every pair
    msi = {}
    for org1 in models:
        stuck_A = len(org_info_single_wo_trans_rxn[org1.id + '_' + org1.id])
        for org2 in models:
            if org1.id + '_' + org2.id in org_info_pair_wo_trans_rxn.keys():
                stuck_AUB = len(org_info_pair_wo_trans_rxn[org1.id + '_' + org2.id])
                if stuck_A == 0:  msi[org1.id + '_' + org2.id] = 0
                else:  msi[org1.id + '_' + org2.id] = 1 - (stuck_AUB / stuck_A)
    return msi, community_model

def calculate_pairwiseMSI(models, media):
    """
    This function calculates pairwise-MSI for all given microbes.

    Creates a csv file containing the MSI values of all pairs.

    Creates an tsv file containing the list of reaction relieved
    in all acceptor microbes in the presence of corresponding donor microbes.

    :param path: path to all xml files
    :param sd_file: path to txt file containing seed metabolites
    """

    warnings.filterwarnings("ignore")
    msi, community_model = pMSI(models, media)
    msi_output_file = f"results/MSI_{os.path.basename(media).replace('.txt', '')}.csv"
    with open(msi_output_file, 'w') as f:
        header = 'organism,in_the_presence,msi_value\n'
        f.write(header)
        for org1, org2 in combinations(models, 2):
            if org1.id + '_' + org2.id in msi.keys():
                f.write(f"{org1.id},{org2.id},{str(msi[org1.id + '_' + org2.id])}\n")
    print('MSI values are written at:\n', msi_output_file)

def calculate_higherorderMSI(models, media, clusters = 'individual_clusters'):
    community_model = commhelper.build_from_species_models(models)
    comm_util = MSModelUtil(community_model)
    org_info, scope, namemap = find_stuckrxns(model, community, media, len(community))
    org_info = decrypt_orginfo(org_info, namemap)
    org_info_wo_trans_rxn = {i: list(set(org_info[i]) - set(comm_util.transport_list())) for i in org_info}

    with open(f"results/community_unperturbed.csv", 'w') as f:
        for i, diff in org_info_wo_trans_rxn.items():
            f.write(i + ',' + str(len(diff)) + '\n')

    if clusters == 'individual_clusters':
        rem_org_list1, rem_org_list2 = {}, {}
        for i, model in enumerate(models):
            rem_org_list1[i] = model.id  ;  rem_org_list2[i] = model.id
    else:
        cluster_data = pd.read_csv(clusters, sep=',')
        rem_org_list1 = cluster_data.set_index('Cluster').T.to_dict('list')
        for n in rem_org_list1:
            rem_org_list1[n] = [j for j in rem_org_list1[n] if pd.isna(j) is False]
        for n in rem_org_list1:
            rem_org_list1[n] = [cobra.io.read_sbml_model(i).id for i in rem_org_list1[n]]
            # rem_org_list1[n] = [model_ids[model_ids.index(i)] for i in rem_org_list1[n]]
        rem_org_list2 = rem_org_list1.copy()

    for nclus in rem_org_list2:
        rem_org_list2[nclus] = [x.replace('.xml', '') for x in rem_org_list2[nclus]]

    with open(f"results/higher_order_msi.csv", 'w') as f:
        for n in rem_org_list1:
            # os.chdir(path)
            # new_models = model.copy()
            # new_community = glob.glob('*.xml')
            # if not new_community:
            #     new_community = glob.glob('*.sbml')
            # new_community.sort()

            pert_models, pert_community, pert_model_ids = make_perturbed_community(rem_org_list1[n], new_models,
                                                                                   new_community)

            org_info_pert, org_info_pert_wo_trans_rxn, msi = perform_task(
                media, pert_models, transport_rxns, pert_community, org_info_wo_trans_rxn, rem_org_list2[n], n)
            for i in rem_org_list2[n]:
                f.write('Comm,clus_' + str(n) + '#' + i + ',' + str(msi) + '\n')

            if msi:
                relieved, detailed_rel_rxns, rel_rxns_name = find_relievedrxns(pert_models, org_info, org_info_pert)
                with open(f'results/clusterKO_/data_analysis/relieved_rxns_Comm--clus{n}.tsv', 'w') as g:
                    write_relieved_rxns(g, relieved, detailed_rel_rxns, rel_rxns_name)
                with open(f'results/clusterKO_/data_analysis/Comm--clus{n}.tsv', 'w') as h:
                    h.write('Comm--clus' + str(n) + '\n')
                    for i in rem_org_list2[n]:
                        h.write(i + '\n')
                    h.write('num of rxns relieved in the below orgs in the presence of clust' + str(n) + '\n')
                    h.write('org,unpert,clust_' + str(n) + 'KO,rxns relieved\n')
                    write_relieved_rxn_metadata(h, org_info_wo_trans_rxn, org_info_pert_wo_trans_rxn)
                print('Comm--clus' + str(n))

            new_models = model.copy()
            new_community = glob.glob('*.xml')
            if not new_community:
                new_community = glob.glob('*.sbml')
            new_community.sort()
            ko_models, ko_community, model_ids = make_perturbed_community(pert_model_ids, new_models, new_community)
            ko_org_list = [x for x in pert_model_ids]
            if len(ko_org_list) < len(model):
                org_info_pert, org_info_pert_wo_trans_rxn, msi = perform_task(
                    media, ko_models, transport_rxns, ko_community, org_info_wo_trans_rxn, ko_org_list, n)
                for i in ko_community:
                    f.write('clus_' + str(n) + '#' + i + ',Comm,' + str(msi) + '\n')

                if msi:
                    relieved, detailed_rel_rxns, rel_rxns_name = find_relievedrxns(ko_models, org_info, org_info_pert)
                    with open(f'results/clusterKO_/data_analysis/relieved_rxns_Comm--clus{n}.tsv', 'w') as g:
                        write_relieved_rxns(g, relieved, detailed_rel_rxns, rel_rxns_name)
                    with open(f'results/clusterKO_/data_analysis/Comm{n}--clus.tsv', 'w') as h:
                        h.write('clus' + str(n) + '--Comm\n')
                        for i in ko_org_list:
                            h.write(i + '\n')
                        h.write('num of rxns relieved in the below orgs in the presence of Comm')
                        h.write('org,unpert,commKO,rxns relieved\n')
                        write_relieved_rxn_metadata(h, org_info_wo_trans_rxn, org_info_pert_wo_trans_rxn)
                    print('clus' + str(n) + '--Comm')
