from collections import OrderedDict
from cobra.core.metabolite import Metabolite
from cobra.io.json import save_json_model
from zipfile import ZipFile, ZIP_LZMA
from pprint import pprint
import logging, json, re, os #, lzma

logging.basicConfig(filename="mscompatability.log", format='%(asctime)s %(message)s', filemode='w')
logger = logging.getLogger(__name__)


# open the parsed ModelSEED Database reactions and compounds content
with open("compound_Xrefs.json" , 'w') as cpdXRefs_out:
    compounds_cross_references = json.load(cpdXRefs_out)
with open("compoundNames.json" , 'w') as cpdNames_out:
    compoundNames = json.load(cpdNames_out)
    
class MSCompatibility():
        
    @staticmethod
    def standardize(self, models, unknown_met_ids, changed_metabolites, changed_reactions, metabolites:bool=True, exchanges:bool=True, conflicts_file_name:str=None, model_names:list=None, export_directory:str=None, view_unknown_mets:bool=True, printing:bool=True):
        changed_ids_count = changed_rxn_count = 0
        for model_index, model in enumerate(models):
            with model:
                # standardize metabolites
                if metabolites:
                    if exchanges:
                        model_metabolites = [met.id for met in model.metabolites]
                        for ex_rxn in model.exchanges:
                            for met in ex_rxn.metabolites:
                                model, met, new_met_id, unknown_met_ids, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility._fix_met(
                                    model, model_index, met, unknown_met_ids, changed_metabolites, changed_reactions, True, printing)
                                changed_rxn_count += rxn_count; changed_ids_count += ids_count
                                try:
                                    ex_rxn.id = 'EX_'+met.id
                                except:
                                    ex_rxn.id = 'EX_'+new_met_id
                                if 'cpd' not in met.id and rxn_count+ids_count != 0 and new_met_id not in model_metabolites: 
                                    unknown_met_ids.append(met.id)
                                    logger.warning(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                    else:
                        for met in model.metabolites:
                            model, met, new_met_id, unknown_met_ids, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility._fix_met(
                                model, model_index, met, unknown_met_ids, changed_metabolites, changed_reactions, True, printing)
                            changed_rxn_count += rxn_count; changed_ids_count += ids_count
                            if 'cpd' not in met.id: 
                                unknown_met_ids.append(met.id)
                                logger.warning(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                    
                    if conflicts_file_name is not None:
                        MSCompatibility._export(
                            models, {'metabolite_changes':changed_metabolites, 'reaction_changes':changed_reactions},
                            conflicts_file_name, model_names, export_directory)
                            
                # standardize reactions
                # else:  # TODO explore standardizing the model reactions, since using modelreactions appears incorrect
                    # with open("reactionIDs.json" , 'w') as rxnID_out:
                    #     reactionIDs = json.load(rxnID_out)
                #     modelreactions_ids = {re.sub('(_\w\d$)', '', rxn['id']).removeprefix('R-'):rxn for rxn in model.modelreactions}
                #     with open(os.path.join(export_directory, 'modelreactions.json'), 'w') as out:
                #         json.dump(modelreactions_ids, out, indent = 3)
                #     model_metabolites = {met.id:met for met in model.metabolites}
                #     missed_reactions = 0
                #     for rxn in model.reactions:
                #         if 'EX_' in rxn.id:
                #             continue
                #         original_reaction = rxn.reaction
                #         rxn.add_metabolites({rxn_met:0 for rxn_met in rxn.metabolites}, combine = False)
                        
                #         if re.sub('(_\w\d$)', '', rxn.id) in modelreactions_ids:
                #             reaction_dict = _parse_modelReactionReagents(
                #                 modelreactions_ids[re.sub('(_\w\d$)', '', rxn.id)]['modelReactionReagents'], model_metabolites
                #                 )
                #         elif rxn.id in modelreactions_ids:
                #             reaction_dict = _parse_modelReactionReagents(
                #                 modelreactions_ids[rxn.id]['modelReactionReagents'], model_metabolites
                #                 )
                #         else:
                #             logger.warning(f'ModelSEEDError: The reaction ID {rxn.id} is not captured by the modelreactions.')
                        
                #         try:
                #             rxn.add_metabolites(reaction_dict, combine = False)  
                #         except:
                #             new_reaction_dict = {}
                #             for met, content in reaction_dict.items():
                #                 if isinstance(met, str):
                #                     met = re.sub('_\w\d', '', met)
                #                 else:
                #                     if re.sub('_\w\d', '', met.id) not in model.metabolites:
                #                         met.id = re.sub('_\w\d', '', met.id)
                #                 new_reaction_dict[met] = content
                #             reaction_dict = new_reaction_dict
                #         if rxn.id not in reactionIDs: 
                #             missed_reactions += 1
                #             # logger.warning(f'ModelSEEDError: The {rxn.id} | {rxn.name} reaction is not recognized by the ModelSEED Database')
                
                #         # describe the change
                #         if original_reaction != rxn.reaction:
                #             change = {
                #                     'original': {
                #                             'reaction': original_reaction
                #                         },
                #                     'new': {
                #                             'reaction': rxn.reaction
                #                         },
                #                     'explanation': f'The reaction {rxn.id} was reconstructed from the ModelSEED Database.'
                #                 }
                #             changed_reactions.append(change)
                        
                #     if export_directory is not None:
                #         with open(os.path.join(export_directory, 'standardized_reactions.txt'), 'w') as out:
                #             json.dump(changed_reactions, out, indent = 3)
                            
                #     total_reactions = 0
                #     for model in models:
                #         total_reactions += len(model.reactions)
                        
                #     logger.warning(f'\nModelSEEDError: {missed_reactions}/{total_reactions} reactions were not captured by the ModelSEED modelreaction IDs.')
            
                models[model_index] = model
        print(f'\n\n{changed_rxn_count} reactions were substituted and {changed_ids_count} metabolite IDs were redefined.')
        if view_unknown_mets:
            return models, unknown_met_ids
        return models
       
    @staticmethod
    def align_exchanges(models, standardize:bool=False, conflicts_file_name:str=None, model_names:list=None, export_directory:str=None, printing:bool=True, extras=False): 
        unknown_met_ids, changed_metabolites, changed_reactions, unique_names, established_mets = [], [], [], [], []
        changed_ids_count = changed_rxn_count = 0
        if standardize:
            standardize(models, unknown_met_ids, changed_metabolites, changed_reactions, 
                        conflicts_file_name='standardized_exchange_metabolites.json', model_names=model_names, printing=printing)
            
        unique_mets, met_conflicts = OrderedDict(), OrderedDict()
        for model_index, model in enumerate(models):
            with model:
                model_metabolites = {met.id:met for met in model.metabolites}
                for ex_rxn in model.exchanges:
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
                                former_model_index = list(unique_mets[former_id].keys())[0].split('_')[0].removeprefix('model')
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
                                model, met, new_met_id, unknown_met_ids, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility._fix_met(
                                    model, model_index, met, unknown_met_ids, changed_metabolites, changed_reactions, standardize, printing)
                                changed_rxn_count += rxn_count; changed_ids_count += ids_count
                        else:
                            former_name = unique_names[list(unique_mets.keys()).index(met.id)]
                            former_model_index = list(unique_mets[met.id].keys())[0].split('_')[0].removeprefix('model')
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
                                model, met, new_met_id, unknown_met_ids, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility._fix_met(
                                    model, model_index, met, unknown_met_ids, changed_metabolites, changed_reactions, standardize, printing)
                                changed_rxn_count += rxn_count; changed_ids_count += ids_count
                 
                        models[model_index] = model
                        
                    # correct the reaction ID
                    if re.sub('(_\w\d$)', '', ex_rxn.id).removeprefix('EX_') in model_metabolites:
                        suffix = re.search('(_\w\d$)', ex_rxn.id).group()
                        model, met, new_met_id, unknown_met_ids, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility._fix_met(
                            model, model_index, re.sub('(_\w\d$)', '', ex_rxn.id).removeprefix('EX_'), 
                            unknown_met_ids, changed_metabolites, changed_reactions, standardize, printing)
                        changed_rxn_count += rxn_count; changed_ids_count += ids_count
                        ex_rxn.id = 'EX_'+new_met_id+suffix

        if conflicts_file_name:
            export_met_conflicts = {}
            for met_id, content in met_conflicts.items():
                export_met_conflicts[met_id] = {}
                for key, val in content.items():
                    if '_met' not in key:
                        export_met_conflicts[met_id][key] = val
                    else:
                        export_met_conflicts[met_id][key.replace('_met','_formula')] = val.formula
                        
            MSCompatibility._export(models, export_met_conflicts, conflicts_file_name, model_names, export_directory)

        print(f'\n\n{changed_rxn_count} exchange reactions were substituted and {changed_ids_count} exchange metabolite IDs were redefined.')
        if extras:
            return models, (unique_mets, unknown_met_ids, changed_metabolites, changed_reactions)
        return models
    
    def _fix_met(model, model_index, met, compound_names, unknown_met_ids, changed_metabolites, changed_reactions, standardize, printing):
        # correct the conflict
        base_name = ''.join(met.name.split('-')[1:]).capitalize()
        met_name = re.sub('_\w\d$', '', met.name)
        new_met_id = met.id
        if met.name in compound_names:
            model, met, new_met_id, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility.__correct_met(
                model, model_index, met, met.name, changed_metabolites, changed_reactions, standardize, printing)
        elif met.name.capitalize() in compound_names:
            model, model_index, met, new_met_id, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility.__correct_met(
                model, met, met.name.capitalize(), changed_metabolites, changed_reactions, standardize, printing)
        elif met_name in compound_names:
            model, model_index, met, new_met_id, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility.__correct_met(
                model, met, met_name, changed_metabolites, changed_reactions, standardize, printing)
        elif met_name.capitalize() in compound_names:
            model, model_index, met, new_met_id, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility.__correct_met(
                model, met, met_name.capitalize(), changed_metabolites, changed_reactions, standardize, printing)
        elif base_name in compound_names and base_name != '':
            model, model_index, met, new_met_id, changed_mets, changed_rxns, rxn_count, ids_count = MSCompatibility.__correct_met(
                model, met, base_name, changed_metabolites, changed_reactions, standardize, printing)
        else:
            unknown_met_ids.append(met.id)
            rxn_count = ids_count = 0
            logger.warning(f'ModelSEEDError: The metabolite ({" | ".join([x for x in [met.id, met.name, base_name, met_name] if x != ""])}) is not recognized by the ModelSEED Database')    
        return model, met, new_met_id, unknown_met_ids, changed_mets, changed_rxns, rxn_count, ids_count
    
    def _export(self, models,
                conflicts,       # the conflicts dictionary that will be exported
                conflicts_file_name,   # the metabolite conflicts are stored and organized, where None does not the conflicts
                model_names,           # specifies the name of the exported model, where None does not export the models
                export_directory       # specifies the directory to which all of the content will be exported 
                ):
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
        
    def __correct_met(model, model_index, met, met_name, changed_metabolites, changed_reactions, standardize, printing):
        def check_cross_references(met, general_met):
            met_refs = compounds_cross_references[compoundNames[met_name]]
            match = False
            for db in compounds_cross_references[general_met]:
                for cross_ref in compounds_cross_references[general_met][db]:
                    if db in met_refs and cross_ref in met_refs[db]:
                        match = True
                        break
                if match:
                    break
            return match, db
        
        changed_rxn_count = changed_ids_count = 0
        original_id = new_met_id = met.id
        compartment = re.search('(_\w\d$)', met.id).group()
        if met.id.removesuffix(compartment) != compoundNames[met_name]:  # If the ID associated with the name deviates from that in the ModelSEED Database
            new_met_id = compoundNames[met_name]+compartment
            if new_met_id in met.model.metabolites:
                # replace the undesirable isomer in every instance, since it cannot be renamed
                for rxn in met.reactions:
                    double_reagent = False
                    original_reaction = rxn.reaction
                    removal_dict, reaction_dict = {}, {}
                    for rxn_met in rxn.reactants+rxn.products:  # The REACTANTS+PRODUCTS may resolve metabolites that are both, more than the METABOLITES attribute
                        match = False
                        stoich = float(rxn.metabolites[rxn_met])
                        compartment = re.search('(_\w\d$)', rxn_met.id).group()
                        new_met = rxn_met
                        if rxn_met.id == met.id:
                            if new_met_id in [old_met.id for old_met in rxn.metabolites]:
                                double_reagent = True
                                logger.warning(f'CodeError: The metabolite {new_met_id} replacement for {met.id} already exists in the reaction {rxn.id}, thus the reaction cannot be updated.')
                                break
                            
                            # affirm the match with cross-references, where it is possible for ModelSEED compounds
                            general_met = re.sub("(_\w\d$)", "", met.id)
                            if 'cpd' in met.id and compounds_cross_references[general_met] != {}:
                                match, db = check_cross_references(met, general_met)
                                if not match:
                                    logger.warning(f'ModelSEEDError: The old metabolite {met.id} cross-references ({compounds_cross_references[general_met]}) do not overlap with those ({compounds_cross_references[compoundNames[met_name]]}) of the new metabolite {new_met_id}.')
                                    
                            # remove duplicate exchange reaction
                            if 'EX_' in rxn.id and 'EX_'+new_met_id in [ex_rxn.id for ex_rxn in model.exchanges]:
                                change = {
                                        'original': {
                                                'reaction': original_reaction
                                            },
                                        'new': {
                                                'reaction': None
                                            },
                                        'justification': f'A {new_met_id} exchange reaction already exists in model {model_index}, thus this duplicative exchange reaction ({rxn.id}) is deleted.'
                                    }
                                if match:
                                    change['justification'] += f' The ID match was verified with {db} cross-references.'
                                model.remove_reactions([rxn.id])
                                changed_reactions.append(change)
                                if printing:
                                    print('\n')
                                    pprint(change, sort_dicts=False)
                                changed_rxn_count += 1
                                double_reagent = True
                                break
                                
                            # define the metabolite with the new name
                            new_met = Metabolite(
                                id = new_met_id, 
                                name = met_name,
                                formula = met.formula,
                                charge = met.charge,
                                compartment = met.compartment
                                )
                            
                        removal_dict[rxn_met] = 0
                        reaction_dict[new_met] = stoich
                        
                    # reconstruct the reactions
                    if double_reagent:
                        continue
                    new_reactants = 0
                    for key, val in reaction_dict.items():
                        new_reactants += 1 if val < 0 else 0
                    new_products = len(reaction_dict) - new_reactants
                    num_reactants, num_products = len(rxn.reactants), len(rxn.products)
                    if num_reactants == new_reactants and num_products == new_products:
                        rxn.add_metabolites(removal_dict, combine = False)
                        rxn.add_metabolites(reaction_dict, combine = False)  
                        change = {
                                'original': {
                                        'reaction': original_reaction
                                    },
                                'new': {
                                        'reaction': rxn.reaction
                                    },
                                'justification': f'The {new_met_id} replacement for {met.id} already exists in model {model_index}, so each reaction (here {rxn.id}) must be updated.'
                            }
                        if match:
                            change['justification'] += f' The ID match was verified with {db} cross-references.'
                        changed_reactions.append(change)
                        if printing:
                            print('\n')
                            pprint(change, sort_dicts=False)
                            
                        changed_rxn_count += 1
                    else:
                        logger.warning(f'CodeError: The reaction {reaction_dict} | {new_reactants} {new_products} possesses a different number of reagents than the original reaction {original_reaction} | {num_reactants} {num_products}, and is skipped.')
            else:
                # affirm the match with cross-references, where it is possible for ModelSEED compounds
                match = False
                general_met = re.sub("(_\w\d$)", "", met.id)
                if 'cpd' in met.id and compounds_cross_references[general_met] != {}:
                    match, db = check_cross_references(met, general_met)
                    if not match:
                        logger.warning(f'ModelSEEDError: The old metabolite {met.id} cross-references ({compounds_cross_references[general_met]}) do not overlap with those ({compounds_cross_references[compoundNames[met_name]]}) of the new metabolite {new_met_id}.')
                            
                # rename the undesirable isomer
                met.id = compoundNames[met_name]+compartment
                change = {
                        'original': {
                                'id': original_id,
                                'name': met.name
                            },
                        'new': {
                                'id': met.id,
                                'name': met_name+compartment
                            },
                        'justification': f'The {original_id} and {met.id} distinction in {model_index} is incompatible.'
                    }
                if 'cpd' not in original_id:
                    change['justification'] = f'The {original_id} ID is not a ModelSEED Database ID.'
                if standardize:
                    change['justification'] = f'The {original_id} and {met.id} metabolites were matched via their name.'
                if match:
                    change['justification'] += f' The ID match was verified with {db} cross-references.'
                changed_metabolites.append(change)
                if printing:
                    print('\n')
                    pprint(change, sort_dicts=False)
                changed_ids_count += 1

        return model, met, new_met_id, changed_metabolites, changed_reactions, changed_rxn_count, changed_ids_count
            