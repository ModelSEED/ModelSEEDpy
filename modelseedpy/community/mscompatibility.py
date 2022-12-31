from collections import OrderedDict, namedtuple
from cobra import Reaction
from cobra.io.json import save_json_model
from modelseedpy.core.fbahelper import FBAHelper
from typing import Iterable
from zipfile import ZipFile, ZIP_LZMA
from math import isclose
from icecream import ic
from pprint import pprint
from deepdiff import DeepDiff
import platform, logging, json, re, os #, lzma

ic.configureOutput(includeContext=True, contextAbsPath=False)

logging.basicConfig(filename="mscompatability.log", format='%(asctime)s %(message)s', filemode='w')
logger = logging.getLogger(__name__)


# open the parsed ModelSEED Database reactions and compounds content
with open(os.path.join(os.path.dirname(__file__), "..", "data", "compound_Xrefs.json"), 'r') as cpdXRefs:
    compounds_cross_references = json.load(cpdXRefs)
with open(os.path.join(os.path.dirname(__file__), "..", "data", "compoundNames.json"), 'r') as cpdNames:
    compoundNames = json.load(cpdNames)

# define generic helper functions
def _remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    return string

def _remove_suffix(string, suffix):
    if string.endswith(suffix):
        return string[:-len(suffix)]
    return string

def _print_changes(change):
    print('\n')
    if float(platform.python_version()[:3]) >= 3.8:
        pprint(change, sort_dicts=False)
    else:
        pprint(change)

def _define_vars(*variables):
    return [var or [] for var in variables]

def _check_cross_references(met, general_met, met_name): 
    met_refs = compounds_cross_references[compoundNames[met_name]]
    matches = []
    for db, content in compounds_cross_references[general_met].items():
        for cross_ref in content:
            if db in met_refs and cross_ref in met_refs[db]:
                matches.append(db)
    return matches


# define a results object
resultsTup = namedtuple("resultsTup", ("new_met_id", "unknown_met_id", "changed_mets", "changed_rxns"))
        
    
class MSCompatibility:
        
    @staticmethod
    def standardize(models, metabolites:bool=True, exchanges:bool=True, conflicts_file_name:str=None,
                    model_names:list=None, export_directory:str=None, view_unknown_mets:bool=False, printing:bool=True,
                    unknown_mets:Iterable=None, changed_mets:Iterable=None, changed_rxns:Iterable=None):
        unknown_mets, changed_mets, changed_rxns = _define_vars(unknown_mets, changed_mets, changed_rxns)
        new_models = []
        models = [models] if not isinstance(models, (list, tuple, set)) else models
        for org_model in models:  # Develop a singular model function and then an abstracted version for multiple models
            model = org_model.copy()  # model_util cannot be used for a circular import
            model_exchanges = [rxn for rxn in model.reactions if "EX_" in rxn.id]
            reactions = {}
            # standardize metabolites
            if not metabolites:
                # TODO develop a correction of reactions based upon their stoichiometry
                continue
            ## standardize exchanges
            if exchanges:
                if printing:
                    message = f"\n\n\nStandardize exchange reactions in {model.id}"
                    print(message, "\n", "="*len(message))
                ex_mets = [met.id for ex_rxn in model_exchanges for met in ex_rxn.metabolites]
                ### correct each metabolite in each exchange reaction
                for ex_rxn in model_exchanges:
                    for met in ex_rxn.metabolites:
                        model, met, reactions, results = MSCompatibility._correct_met(
                            model, met, reactions, True, printing)
                        new_rxn_id = 'EX_'+results.new_met_id
                        if all(['cpd' not in met.id, results.new_met_id not in ex_mets, results.unknown_met_id]):
                            unknown_mets.append(met.id)
                            logger.warning(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                        changed_rxns.extend(results.changed_rxns) ; changed_mets.extend(results.changed_mets)
                    ex_rxn.id = new_rxn_id
            else:
                for met in model.metabolites:
                    model, met, reactions, results = MSCompatibility._correct_met(model, met, reactions, True, printing)
                    if 'cpd' not in met.id:
                        unknown_mets.append(met.id)
                        logger.warning(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')

                if conflicts_file_name is not None:
                    model_names = model_names or [model.id for model in models]
                    MSCompatibility._export(
                        models, {'metabolite_changes':changed_mets, 'reaction_changes':changed_rxns},
                        conflicts_file_name, model_names, export_directory)
            new_models.append(model)
            MSCompatibility._validate_results(model, org_model, unknown_mets)
        models_id = ",".join([model.id for model in models])
        if len(changed_rxns) == len(changed_mets) == 0:
            print(f"The {'exchange ' if exchanges else ''} metabolite ID's of the model {models_id} "
                  f"are completely standardized to ModelSEED.")
        else:
            print(f'\n\n{len(changed_rxns)} reactions were substituted and '
                  f'{len(changed_mets)} metabolite IDs were redefined in {models_id} by standardize().')
        if view_unknown_mets:
            return new_models, unknown_mets
        return new_models
       
    @staticmethod
    # TODO - verify that this method appropriately aligns the exchange reactions of the two models
    def align_exchanges(models, conflicts_file_name:str=None, model_names:list=None, export_directory:str=None, printing:bool=True, extras=False): 
        unknown_mets, changed_mets, changed_rxns, unique_names, established_mets = [], [], [], [], []
        unique_mets, met_conflicts = OrderedDict(), OrderedDict()
        new_models = []
        for model_index, org_model in enumerate(models):
            model = org_model.copy()
            model_exchanges = [rxn for rxn in model.reactions if "EX_" in rxn.id]
            message = f"\n\n\nAlign exchange reactions in {model.id}"
            print(message, "\n", "="*len(message))
            model_metabolites = {met.id:met for met in model.metabolites}
            reactions = {}
            for ex_rxn in model_exchanges:
                for met in ex_rxn.metabolites:
                    met_name = re.sub('_\w\d$', '', met.name) 
                    if met.id not in unique_mets and met.id not in established_mets: 
                        if met_name not in unique_names:
                            # identify the unique metabolite
                            unique_mets[met.id] = {
                                f'model{model_index}_id': met.id,
                                f'model{model_index}_met': met
                                }
                            unique_names.append(met_name)
                        else:
                            # describe the metabolite conflict between the ID and name
                            former_id = list(unique_mets.keys())[unique_names.index(met_name)]
                            former_model_index = _remove_prefix(list(unique_mets[former_id].keys())[0].split('_')[0], 'model')
                            if met.name not in met_conflicts:
                                met_conflicts[met_name] = {
                                        f'model{former_model_index}_id': former_id,
                                        f'model{former_model_index}_met': unique_mets[former_id][f'model{former_model_index}_met'],
                                        f'model{model_index}_id': met.id,
                                        f'model{model_index}_met': met
                                    }
                            else:
                                met_conflicts[met_name].update({
                                        f'model{model_index}_id': met.id,
                                        f'model{model_index}_met': met
                                    })
                            model, met, reactions, results = MSCompatibility._correct_met(model, met, reactions, False, printing)
                    else:
                        former_name = unique_names[list(unique_mets.keys()).index(met.id)]
                        former_model_index = _remove_prefix(list(unique_mets[met.id].keys())[0].split('_')[0], 'model')
                        if met_name == former_name:
                            # remove the metabolite that is no longer unique
                            del unique_names[list(unique_mets.keys()).index(met.id)]
                            unique_mets.pop(met.id)    
                            established_mets.append(met.id)
                        else:
                            # describe the conflicting metabolite names
                            if met.id not in met_conflicts:
                                met_conflicts[met.id] = {
                                        f'model{former_model_index}_name': former_name,
                                        f'model{former_model_index}_met': unique_mets[met.id][f'model{former_model_index}_met'],
                                        f'model{model_index}_name': met.name,
                                        f'model{model_index}_met': met
                                    }
                            else:
                                if f'model{model_index}_name' not in met_conflicts[met.id]:
                                    met_conflicts[met.id].update({
                                            f'model{model_index}_name': met.name,
                                            f'model{model_index}_met': met
                                        })
                                else:
                                    iteration = 0
                                    while f'model{model_index}_{iteration}_name' in met_conflicts[met.id]:
                                        iteration += 1
                                        
                                    met_conflicts[met.id].update({
                                            f'model{model_index}_{iteration}_name': met.name,
                                            f'model{model_index}_{iteration}_met': met
                                        })
                            model, met, reactions, results = MSCompatibility._correct_met(model, met, reactions, False, printing)
                    
                # correct the reaction ID
                reaction = _remove_prefix(re.sub('(_\w\d$)', '', ex_rxn.id), 'EX_')
                if reaction in model_metabolites:
                    suffix = re.search('(_\w\d$)', reaction).group()
                    model, met, reactions, results = MSCompatibility._correct_met(
                        model, _remove_suffix(reaction, suffix), reactions, False, printing)
                    ex_rxn.id = 'EX_'+results.new_met_id+suffix
            new_models.append(model)
            MSCompatibility._validate_results(model, org_model, unknown_mets)

        if conflicts_file_name:
            export_met_conflicts = {}
            for met_id, content in met_conflicts.items():
                export_met_conflicts[met_id] = {}
                for key, val in content.items():
                    if "_met" not in key:
                        export_met_conflicts[met_id][key] = val
                    else:
                        export_met_conflicts[met_id][key.replace('_met','_formula')] = val.formula
                        
            model_names = model_names or [model.id for model in models]
            MSCompatibility._export(new_models, export_met_conflicts, conflicts_file_name, model_names, export_directory)

        print(f'\n\n{len(changed_rxns)} exchange reactions were substituted and '
              f'{len(changed_mets)} exchange metabolite IDs were redefined by align_exchanges().')
        if extras:
            return models, (unique_mets, unknown_mets, changed_mets, changed_rxns)
        return models

    @staticmethod
    # !!! This does not catch the errors, perhaps from faulty unknown_mets
    def _validate_results(model, org_model, unknown_mets, standardize=True):
        # ensure that all non-standard exchanges have been corrected
        model_exchanges = [rxn for rxn in model.reactions if "EX_" in rxn.id]
        if standardize:
            residual_nonstandard_mets = [met.id for ex_rxn in model_exchanges for met in ex_rxn.metabolites if "cpd" not in met.id]
            residuals = set(residual_nonstandard_mets)-set(unknown_mets)
            if residuals:
                logger.error(f"The {model.id} model has residual non-standard metabolites in its exchange reactions:"
                             f" {residuals}. Specify a True `printing` parameter to view why these metabolites were not corrected.")
        else:  # TODO develop a check for aligned_exchanges that ID's match between models
            pass

        # verify that no duplicate reactions were added to the model
        reactions = [rxn.name for rxn in model.variables]
        if len(reactions) != len(set(reactions)):
            duplicate_reactions = DeepDiff(set(reactions), reactions)
            logger.critical(f'CodeError: The model {org_model.id} contains {duplicate_reactions}'
                            f' that compromise the model.')

        # verify that the objective value is practically unchanged
        original_objective_value = org_model.slim_optimize()
        new_objective_value = model.slim_optimize()
        if not isclose(original_objective_value, new_objective_value, rel_tol=1e-6):
            logger.critical(f"The original objective value {original_objective_value} for {org_model.id}"
                            f" does not equal the new objective value {new_objective_value}.")

    @staticmethod
    def _export(models, conflicts, conflicts_file_name, model_names, export_directory):
        if export_directory is None:
            export_directory = os.getcwd()

        file_paths = []
        if conflicts_file_name:
            path = os.path.join(export_directory,conflicts_file_name)
            file_paths.append(os.path.relpath(path, export_directory))
            with open(path, 'w') as out:
                json.dump(conflicts, out, indent = 3)
        if model_names:
            for index, model in enumerate(models):
                path = os.path.join(export_directory,f'{model_names[index]}.json')
                file_paths.append(os.path.relpath(path, export_directory))
                save_json_model(model, path)
        with ZipFile(
            "_".join(model_names[:4]) + ".zip", "w", compression=ZIP_LZMA
        ) as zip:
            for file in file_paths:
                zip.write(file)
                os.remove(file)

    @staticmethod
    def _correct_met(model, met, reactions, standardize, printing):
        # identify a matching metabolite name in the ModelSEED Database
        base_name = ''.join(met.name.split('-')[1:]).capitalize()
        comp = re.compile("(_\w\d+$)")
        general_name = comp.sub('', met.name)
        general_met = comp.sub("", met.id)
        met_name = None
        for possible_name in [met.name, met.name.capitalize(), general_name, general_name.capitalize(), base_name]:
            if possible_name in compoundNames:
                met_name = possible_name
                break
        if not met_name:
            # TODO - add a search through cross-references to confirm that the metabolite is unknown, despite different names
            # general_met = re.sub("(_\w\d+$)", "", met.id)
            # matches = _check_cross_references(met, general_met, None)
            # if not matches:
            #     logger.warning(f"ModelSEEDError: The old metabolite {met.id} cross-references"
            #     f" ({compounds_cross_references[general_met]}) do not overlap with those"
            #     f" ({compounds_cross_references[compoundNames[met_name]]}) of the new metabolite {new_met_id}.")
            metabolite_desc = " | ".join([x for x in [met.id, met.name, base_name, general_name] if x != ""])
            logger.warning(f"ModelSEEDError: The metabolite ({metabolite_desc}) is not recognized by the ModelSEED Database")
            return model, met, reactions, resultsTup(met.id, met.id, [], [])
        compartment = comp.search(met.id).group()
        if _remove_suffix(met.id, compartment) == compoundNames[met_name]:
            return model, met, reactions, resultsTup(met.id, None, [], [])

        # correct the metabolite whose name matches a metabolite in the ModelSEED Database but whose ID deviates
        changed_mets, changed_rxns, matches = [], [], []
        new_met_id = original_id = met.id
        original_name = met.name 
        ## affirm the match with cross-references, where it is possible for ModelSEED compounds
        if 'cpd' in met.id:  # TODO correct an anomaly where a valid MSID is repeated with itself
            logger.warning(f"IDWarning: The original ID {met.id} is a ModelSEED ID, and "
                           f"may not be desirably changed to {new_met_id}.")
            if compounds_cross_references[general_met] != {}:
                matches = _check_cross_references(met, general_met, met_name)
                if not matches:
                    logger.warning(f"ModelSEEDError: The old metabolite {met.id} cross-references"
                    f" ({compounds_cross_references[general_met]}) do not overlap with those"
                    f" ({compounds_cross_references[compoundNames[met_name]]}) of the new metabolite {new_met_id}.")
        new_met_id = compoundNames[met_name]+compartment
        if new_met_id in model.metabolites:
            ## replace the undesirable compound in each of its reactions, since the new compound ID already exists
            for org_rxn in met.reactions:
                original_reaction = org_rxn.reaction
                reaction_dict = {}
                ### remove the old exchange reaction
                if org_rxn.id == 'EX_'+new_met_id:
                    change = {'original': {'reaction': original_reaction},
                              'new': "-- Deleted --",
                              'justification': f"A {new_met_id} exchange reaction already exists in model {model.id},"
                              f" thus this reaction ({org_rxn.id}) must be deleted to permit its replacement."}
                    if matches:
                        change['justification'] += f' The ID match was verified with the {matches} cross-reference(s).'
                    model.remove_reactions([org_rxn])
                    changed_rxns.append(change)
                    if printing:
                        _print_changes(change)
                ### a new reaction is created in all of its instances
                else:
                    # !!! duplicate reactions and hindrance to multiple changes on the same reaction remain challenges.
                    #### Copy the reaction stoichiometry while introducing the new metabolite
                    rxn = org_rxn if org_rxn.id not in reactions else reactions[org_rxn.id]
                    reaction_met_ids = {}
                    redundant_mets = False
                    for rxn_met in rxn.metabolites:
                        stoich = float(rxn.metabolites[rxn_met])
                        new_met = rxn_met.copy()
                        if rxn_met.id == met.id:
                            new_met.id, new_met.name = new_met_id, met_name
                            if new_met not in model.metabolites:
                                model.add_metabolites(new_met)
                        if new_met.id in reaction_met_ids:
                            redundant_mets = True
                            reaction_dict[reaction_met_ids[new_met.id]] += stoich
                        else:
                            reaction_met_ids[new_met.id] = new_met
                            reaction_dict[new_met] = stoich

                    #### replace the original reaction with a renewed reaction
                    new_reactants = sum([1 for val in reaction_dict.values() if val < 0])
                    new_products = len(reaction_dict) - new_reactants
                    if len(rxn.reactants) == new_reactants and len(rxn.products) == new_products or redundant_mets:
                        new_rxn = Reaction(id=rxn.id, name=rxn.name, subsystem=rxn.subsystem,
                                           lower_bound=rxn.lower_bound, upper_bound=rxn.upper_bound)
                        model.remove_reactions([rxn])
                        model.add_reaction(new_rxn)
                        new_rxn.add_metabolites(reaction_dict)
                        change = {'original': {'reaction': original_reaction},
                                  'new': {'reaction': new_rxn.reaction},
                                  'justification': f"The new {new_met_id} ID for {met.id} already exists in model"
                                                   f" {model.id}, so each reaction (here {rxn.id}) must be replaced."}
                        if matches:
                            change['justification'] += f' The ID match was verified with the {matches} cross-reference(s).'
                        changed_rxns.append(change)
                        if printing:
                            _print_changes(change)
                        reactions[org_rxn.id] = new_rxn
                    else:
                        ##### print the discrepancy in the reagents from the original reaction
                        rxn_diff = DeepDiff(FBAHelper.IDRxnMets(org_rxn), FBAHelper.IDRxnMets(reaction_dict))
                        logger.error(f"CodeError: The new reaction of {rxn.id} with"
                                     f" {new_reactants} reactants | {new_products} products"
                                     f" differs from the original reaction with "
                                     f"{len(rxn.reactants)} reactants | {len(rxn.products)} products,"
                                     f" {rxn_diff} and is therefore skipped.")
        else:
            ## rename the undesirable compound
            met.name = met_name+compartment
            met.id = new_met_id
            change = {'original': {'id': original_id, 'name': original_name},
                      'new': {'id': met.id, 'name': met.name},
                      'justification': f'The {original_id} and {met.id} distinction in {model.id} is incompatible; '
                                       f'hence, the {met.id} ID and {met.name} are used.'}
            if 'cpd' not in original_id:
                change['justification'] += f' The {original_id} ID is not a ModelSEED Database ID.'
            if standardize:
                change['justification'] += f' The {original_id} and {met.id} metabolites were matched via their name.'
            if matches:
                change['justification'] += f' The ID match was verified with the {matches} cross-reference(s).'
            changed_mets.append(change)
            if printing:
                _print_changes(change)
                    
        return model, met, reactions, resultsTup(new_met_id, None, changed_mets, changed_rxns)
