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


def remove_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    return string


def remove_suffix(string, suffix):
    if string.endswith(suffix):
        return string[:-len(suffix)]
    return string


def print_changes(change):
    if float(platform.python_version()[:3]) >= 3.8:
        pprint(change, sort_dicts=False)
    else:
        pprint(change)


def define_vars(variables):
    return [var or [] for var in variables]


resultsTup = namedtuple("resultsTup", ("new_met_id", "unknown_met_id", "changed_mets", "changed_rxns"))
        
    
class MSCompatibility:
        
    @staticmethod
    def standardize(models:Iterable, metabolites:bool=True, exchanges:bool=True, conflicts_file_name:str=None, 
                    model_names:list=None, export_directory:str=None, view_unknown_mets:bool=False, printing:bool=True,
                    unknown_met_ids:Iterable=None, changed_metabolites:Iterable=None, changed_reactions:Iterable=None):
        unknown_met_ids, changed_metabolites, changed_reactions = define_vars([unknown_met_ids, changed_metabolites, changed_reactions])
        new_models = []
        for org_model in models:
            model = org_model.copy()
            reactions = {}
            # standardize metabolites
            if metabolites:
                if exchanges:
                    if printing:
                        message = f"\n\n\nStandardize exchange reactions in {model.id}"
                        print(message, "\n", "="*len(message))
                    ex_mets = [met.id for ex_rxn in FBAHelper.exchange_reactions(model) for met in ex_rxn.metabolites]
                    for ex_rxn in FBAHelper.exchange_reactions(model):
                        for met in ex_rxn.metabolites:
                            model, met, reactions, results = MSCompatibility._fix_met(model, met, reactions, True, printing)
                            if results.unknown_met_id:
                                unknown_met_ids.append(results.unknown_met_id)
                            try:  # catching errors of repeated exchange IDs
                                ex_rxn.id = 'EX_'+met.id
                            except:
                                ex_rxn.id = 'EX_'+results.new_met_id
                                
                            if all(['cpd' not in met.id, results.new_met_id not in ex_mets,
                                    not results.changed_rxns+results.changed_mets]):
                                unknown_met_ids.append(met.id)
                                logger.warning(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                            changed_reactions.extend(results.changed_rxns)
                            changed_metabolites.extend(results.changed_mets)
                else:
                    for met in model.metabolites:
                        model, met, results = MSCompatibility._fix_met(model, met, True, printing)
                        if 'cpd' not in met.id: 
                            unknown_met_ids.append(met.id)
                            logger.warning(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                
                if conflicts_file_name is not None:
                    model_names = model_names or [model.id for model in models]
                    MSCompatibility._export(
                        models, {'metabolite_changes':changed_metabolites, 'reaction_changes':changed_reactions},
                        conflicts_file_name, model_names, export_directory)
            new_models.append(model)
            MSCompatibility._validate_results(model, org_model, unknown_met_ids)
        print(f'\n\n{len(changed_reactions)} reactions were substituted and {len(changed_metabolites)} metabolite IDs were redefined by standardize().')
        if view_unknown_mets:
            return new_models, unknown_met_ids
        return new_models
       
    @staticmethod
    # TODO - verify that this method appropriately aligns the exchange reactions of the two models
    def align_exchanges(models, conflicts_file_name:str=None, model_names:list=None, export_directory:str=None, printing:bool=True, extras=False): 
        unknown_met_ids, changed_metabolites, changed_reactions, unique_names, established_mets = [], [], [], [], []
        unique_mets, met_conflicts = OrderedDict(), OrderedDict()
        new_models = []
        for model_index, org_model in enumerate(models):
            model = org_model.copy()
            message = f"\n\n\nAlign exchange reactions in {model.id}"
            print(message, "\n", "="*len(message))
            model_metabolites = {met.id:met for met in model.metabolites}
            reactions = {}
            for ex_rxn in FBAHelper.exchange_reactions(model):
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
                            former_model_index = remove_prefix(list(unique_mets[former_id].keys())[0].split('_')[0], 'model')
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
                            model, met, reactions, results = MSCompatibility._fix_met(model, met, reactions, False, printing)
                    else:
                        former_name = unique_names[list(unique_mets.keys()).index(met.id)]
                        former_model_index = remove_prefix(list(unique_mets[met.id].keys())[0].split('_')[0], 'model')
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
                            model, met, reactions, results = MSCompatibility._fix_met(model, met, reactions, False, printing)
                    
                # correct the reaction ID
                reaction = remove_prefix(re.sub('(_\w\d$)', '', ex_rxn.id), 'EX_')
                if reaction in model_metabolites:
                    suffix = re.search('(_\w\d$)', reaction).group()
                    model, met, reactions, results = MSCompatibility._fix_met(model, remove_suffix(reaction, suffix), reactions, False, printing)
                    ex_rxn.id = 'EX_'+results.new_met_id+suffix
            new_models.append(model)
            MSCompatibility._validate_results(model, org_model, unknown_met_ids)

        if conflicts_file_name:
            export_met_conflicts = {}
            for met_id, content in met_conflicts.items():
                export_met_conflicts[met_id] = {}
                for key, val in content.items():
                    if '_met' not in key:
                        export_met_conflicts[met_id][key] = val
                    else:
                        export_met_conflicts[met_id][key.replace('_met','_formula')] = val.formula
                        
            model_names = model_names or [model.id for model in models]
            MSCompatibility._export(new_models, export_met_conflicts, conflicts_file_name, model_names, export_directory)

        print(f'\n\n{len(changed_reactions)} exchange reactions were substituted and '
              f'{len(changed_metabolites)} exchange metabolite IDs were redefined by align_exchanges().')
        if extras:
            return models, (unique_mets, unknown_met_ids, changed_metabolites, changed_reactions)
        return models

    @staticmethod
    # !!! This does not catch the errors, perhaps from faulty unknown_met_ids
    def _validate_results(model, org_model, unknown_met_ids, standardize=True):
        # ensure that all non-standard exchanges have been corrected
        if standardize:
            residual_nonstandard_mets = [met.id for ex_rxn in FBAHelper.exchange_reactions(model) for met in ex_rxn.metabolites if "cpd" not in met.id]
            residuals = set(residual_nonstandard_mets)-set(unknown_met_ids)
            if residuals:
                logger.error(f"The {model.id} model has residual non-standard metabolites in its exchange reactions: {residuals}."
                             " Specify a True `printing` parameter to view why these metabolites were not corrected.")
        else:  # TODO develop a check for aligned_exchanges between models
            pass

        # verify that no duplicate reactions were added to the model
        reactions = [rxn.name for rxn in model.variables]
        duplicate_reactions = DeepDiff(sorted(reactions), sorted(set(reactions)))
        if duplicate_reactions:
            logger.critical(f'CodeError: The model {model.id} contains {duplicate_reactions}'
                         f' that compromise the model.')

        # verify that the objective value is (essentially) identical to the original model
        original_objective_value = org_model.slim_optimize()
        new_objective_value = org_model.slim_optimize()
        if not isclose(original_objective_value, new_objective_value, rel_tol=1e-6):
            logger.critical(f"The original objective value {original_objective_value} "
                         f"does not equal the new objective value {new_objective_value}.")

    @staticmethod
    def _fix_met(model, met, reactions, standardize, printing):
        # correct the conflict
        base_name = ''.join(met.name.split('-')[1:]).capitalize()
        met_name = re.sub('_\w\d$', '', met.name)
        new_met_id = met.id
        for possible_name in [met.name, met.name.capitalize(), met_name, met_name.capitalize(), base_name]:
            if possible_name in compoundNames:
                model, met, reactions, new_met_id, changed_metabolites, changed_reactions = MSCompatibility._correct_met(
                    model, met, reactions, possible_name, standardize, printing)
                return model, met, reactions, resultsTup(new_met_id, None, changed_metabolites, changed_reactions)
        
        # TODO - add a search through cross-references to confirm that the metabolite is unknown, despite different names
        # general_met = re.sub("(_\w\d+$)", "", met.id)
        # matches = MSCompatibility._check_cross_references(met, general_met, None)
        # if not matches:
        #     logger.warning(f"ModelSEEDError: The old metabolite {met.id} cross-references"
        #     f" ({compounds_cross_references[general_met]}) do not overlap with those"
        #     f" ({compounds_cross_references[compoundNames[met_name]]}) of the new metabolite {new_met_id}.")
        logger.warning(f'ModelSEEDError: The metabolite ({" | ".join([x for x in [met.id, met.name, base_name, met_name] if x != ""])})'
                       " is not recognized by the ModelSEED Database")  
        return model, met, reactions, resultsTup(new_met_id, met.id, [], [])

    @staticmethod
    def _export(models, conflicts, conflicts_file_name, model_names, export_directory):
        if export_directory is None:
            export_directory = os.getcwd()
                    
        file_paths = []
        if conflicts_file_name is not None:
            path = os.path.join(export_directory,conflicts_file_name)
            file_paths.append(os.path.relpath(path, export_directory))
            with open(path, 'w') as out:
                json.dump(conflicts, out, indent = 3)
        if model_names:
            for index, model in enumerate(models):
                path = os.path.join(export_directory,f'{model_names[index]}.json')
                file_paths.append(os.path.relpath(path, export_directory))
                save_json_model(model, path)
        with ZipFile('_'.join(model_names[:4])+'.zip', 'w', compression = ZIP_LZMA) as zip:
            for file in file_paths:
                zip.write(file)
                os.remove(file)

    @staticmethod
    def _check_cross_references(met, general_met, met_name): 
        met_refs = compounds_cross_references[compoundNames[met_name]]
        matches = []
        for db, content in compounds_cross_references[general_met].items():
            for cross_ref in content:
                if db in met_refs and cross_ref in met_refs[db]:
                    matches.append(db)
        return matches

    @staticmethod
    def _correct_met(model, met, reactions, met_name, standardize, printing):
        changed_metabolites, changed_reactions = [], []
        original_id = new_met_id = met.id
        original_name = met.name 
        # affirm the match with cross-references, where it is possible for ModelSEED compounds
        general_met = re.sub("(_\w\d+$)", "", met.id)
        matches = []
        if 'cpd' in met.id and compounds_cross_references[general_met] != {}:
            matches = MSCompatibility._check_cross_references(met, general_met, met_name)
            if not matches:
                logger.warning(f"ModelSEEDError: The old metabolite {met.id} cross-references"
                f" ({compounds_cross_references[general_met]}) do not overlap with those"
                f" ({compounds_cross_references[compoundNames[met_name]]}) of the new metabolite {new_met_id}.")
        compartment = re.search('(_\w\d+$)', met.id).group()
        if remove_suffix(met.id, compartment) != compoundNames[met_name]:  # Assess whether the ID deviates from that associated with the name in the ModelSEED Database
            model_exchange_ids = [ex_rxn.id for ex_rxn in FBAHelper.exchange_reactions(model)]
            new_met_id = compoundNames[met_name]+compartment
            general_met = re.sub("(_\w\d+$)", "", new_met_id)
            if 'cpd' in met.id:
                logger.warning(f"IDWarning: The original ID {met.id} is a ModelSEED ID, and "
                               f"may not be desirably changed to {new_met_id}.")
            if new_met_id in model.metabolites:
                # replace the undesirable isomer in every instance, since it cannot be renamed
                for org_rxn in met.reactions:
                    original_reaction = org_rxn.reaction
                    reaction_dict = {} ; change = None
                    # remove duplicate exchange reaction
                    if 'EX_' in org_rxn.id and 'EX_'+new_met_id in model_exchange_ids:
                        change = {'original': {'reaction': original_reaction},
                                'new': "-- Deleted --",
                                'justification': f"A {new_met_id} exchange reaction already exists in model {model.id},"
                                f" thus this duplicative exchange reaction ({org_rxn.id}) is deleted."}
                        if matches:
                            change['justification'] += f' The ID match was verified with the {matches} cross-reference(s).'
                        model.remove_reactions([org_rxn])
                        changed_reactions.append(change)
                        if printing:
                            print('\n')
                            print_changes(change)
                    else:
                        rxn = org_rxn
                        if org_rxn.id in reactions:
                            rxn = reactions[org_rxn.id]
                        # !!! This clause creates a new reaction with each edit, which can both duplicate reactions and prevent multiple
                        # !!! changes to the same reaction.
                        reaction_met_ids = {}
                        redundant_mets = False
                        for rxn_met in rxn.metabolites:
                            stoich = float(rxn.metabolites[rxn_met])
                            new_met = rxn_met.copy()
                            # create a new metabolite for the replacement reaction
                            if rxn_met.id == met.id:
                                compartment = re.search('(_\w\d$)', rxn_met.id).group()
                                new_met.id = new_met_id ; new_met.name = met_name
                                # add metabolites that are not currently in the model (according to IDs/names)
                                if new_met not in model.metabolites:
                                    model.add_metabolites(new_met)
                            if new_met.id in reaction_met_ids:
                                redundant_mets = True
                                reaction_dict[reaction_met_ids[new_met.id]] += stoich
                            else:
                                reaction_dict[new_met] = stoich
                                reaction_met_ids[new_met.id] = new_met
                            
                        # reconstruct the reactions
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
                                    'justification': f"The new {new_met_id} ID for {met.id} already exists in model {model.id},"
                                    f" so each reaction (here {rxn.id}) must be updated."}
                            if matches:
                                change['justification'] += f' The ID match was verified with the {matches} cross-reference(s).'
                            changed_reactions.append(change)
                            if printing:
                                print('\n')
                                print_changes(change)
                            reactions[org_rxn.id] = new_rxn
                        else:
                            rxn_diff = DeepDiff(FBAHelper.IDRxnMets(org_rxn), FBAHelper.IDRxnMets(reaction_dict))
                            logger.error(f"CodeError: The new reaction of {rxn.id} with"
                                         f" {new_reactants} reactants | {new_products} products"
                                         f" differs from the original reaction with "
                                         f"{len(rxn.reactants)} reactants | {len(rxn.products)} products,"
                                         f" {rxn_diff} and is therefore skipped.")
            else:
                # rename the undesirable isomer
                met.name = met_name+compartment
                met.id = new_met_id
                change = {'original': {'id': original_id, 'name': original_name},
                        'new': {'id': met.id, 'name': met.name},
                        'justification': f'The {original_id} and {met.id} distinction in {model.id} is incompatible.'}
                if 'cpd' not in original_id:
                    change['justification'] += f' The {original_id} ID is not a ModelSEED Database ID.'
                if standardize:
                    change['justification'] += f' The {original_id} and {met.id} metabolites were matched via their name.'
                if matches:
                    change['justification'] += f' The ID match was verified with the {matches} cross-reference(s).'
                changed_metabolites.append(change)
                if printing:
                    print('\n')
                    print_changes(change)
                    
        return model, met, reactions, new_met_id, changed_metabolites, changed_reactions