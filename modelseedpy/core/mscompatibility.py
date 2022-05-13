from collections import OrderedDict
from cobra.core.metabolite import Metabolite
from numpy import negative
from warnings import warn
from pprint import pprint
import json, re, os

class MSCompatibility():
    def __init__(self,
                 modelseed_db_path: str, # the local path to the ModelSEEDDatabase repository
                 printing = True         # specifies whether results are printed
                 ):       
        self.printing = printing
        
        # import and parse ModelSEED Database reactions and compounds
        with open(os.path.join(modelseed_db_path, 'Biochemistry', 'reactions.json'), 'r') as rxns:
            self.reactions = json.load(rxns)
            self.reaction_ids = OrderedDict()
            for rxn in self.reactions:
                self.reaction_ids[rxn['id']] = rxn['name']
                
        with open(os.path.join(modelseed_db_path, 'Biochemistry', 'compounds.json'), 'r') as rxns:
            self.compounds = json.load(rxns)
            self.compounds_id_indexed = self.compound_names = OrderedDict()
            for cpd in self.compounds:
                self.compounds_id_indexed[cpd['id']] = cpd
                if 'aliases' in cpd and cpd['aliases'] is not None:
                    names = [name.strip() for name in cpd['aliases'][0].split(';')]
                    names.append(cpd['name'])
                    for name in names:
                        if name not in self.compound_names:
                            self.compound_names[name] = cpd['id']
    
    # def _parse_modelReactionReagents(self, modelReactionReagents, model_metabolites):
    #     rxn_dict = {}
    #     for cpd in modelReactionReagents:
    #         met = re.search('(?<=id\/)(.+)', cpd['modelcompound_ref']).group()
    #         stoich = float(cpd['coefficient'])
    #         if met in model_metabolites:
    #             met = model_metabolites[met]
    #         elif re.sub('_\w\d', '', met) in model_metabolites:
    #             met = model_metabolites[re.sub('_\w\d', '', met)]
    #         else:
    #             KeyError(f'ModelSEEDError: The metabolite {met} in the reactions is not in the modelreactions.')
    #         rxn_dict[met] = stoich

    #     return rxn_dict
        
    def standardize_MSD(self, models, output_directory = None, metabolites = True, exchanges = True):
        self.unique_mets, self.met_conflicts = OrderedDict(), OrderedDict()
        self.unknown_met_ids, self.changed_metabolites, self.changed_reactions= [], [], []
        if exchanges:
            output = None
            if output_directory is not None:
                os.path.join(output_directory, 'exchanges_conflicts.txt')
            self.correct_exchanges(models, output)
        for model in models:
            model_index = models.index(model)
            
            # standardize metabolites
            if metabolites:
                for met in model.metabolites:
                    met = self._fix_met(met)
                    if 'cpd' not in met.id: 
                        self.unknown_met_ids.append(met.id)
                        warn(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
            
                if output_directory is not None:
                    with open(os.path.join(output_directory, 'standardized_metabolites.txt'), 'w') as out:
                        json.dump(self.changed_metabolites+self.changed_reactions, out, indent = 3)
                        
            # standardize reactions
        #     else:
        #         modelreactions_ids = {re.sub('(_\w\d$)', '', rxn['id']).removeprefix('R-'):rxn for rxn in model.modelreactions}
        #         with open(os.path.join(output_directory, 'modelreactions.json'), 'w') as out:
        #             json.dump(modelreactions_ids, out, indent = 3)
        #         model_metabolites = {met.id:met for met in model.metabolites}
        #         missed_reactions = 0
        #         for rxn in model.reactions:
        #             if 'EX_' in rxn.id:
        #                 continue
        #             original_reaction = rxn.reaction
        #             rxn.add_metabolites({rxn_met:0 for rxn_met in rxn.metabolites}, combine = False)
                    
        #             if re.sub('(_\w\d$)', '', rxn.id) in modelreactions_ids:
        #                 reaction_dict = self._parse_modelReactionReagents(
        #                     modelreactions_ids[re.sub('(_\w\d$)', '', rxn.id)]['modelReactionReagents'], model_metabolites
        #                     )
        #             elif rxn.id in modelreactions_ids:
        #                 reaction_dict = self._parse_modelReactionReagents(
        #                     modelreactions_ids[rxn.id]['modelReactionReagents'], model_metabolites
        #                     )
        #             else:
        #                 warn(f'ModelSEEDError: The reaction ID {rxn.id} is not captured by the modelreactions.')
                    
        #             try:
        #                 rxn.add_metabolites(reaction_dict, combine = False)  
        #             except:
        #                 new_reaction_dict = {}
        #                 for met, content in reaction_dict.items():
        #                     if isinstance(met, str):
        #                         met = re.sub('_\w\d', '', met)
        #                     else:
        #                         if re.sub('_\w\d', '', met.id) not in model.metabolites:
        #                             met.id = re.sub('_\w\d', '', met.id)
        #                     new_reaction_dict[met] = content
        #                 reaction_dict = new_reaction_dict
        #             if rxn.id not in self.reaction_ids: 
        #                 missed_reactions += 1
        #                 # warn(f'ModelSEEDError: The {rxn.id} | {rxn.name} reaction is not recognized by the ModelSEED Database')
            
        #             # describe the change
        #             if original_reaction != rxn.reaction:
        #                 change = {
        #                         'original': {
        #                                 'reaction': original_reaction
        #                             },
        #                         'new': {
        #                                 'reaction': rxn.reaction
        #                             },
        #                         'explanation': f'The reaction {rxn.id} was reconstructed from the ModelSEED Database.'
        #                     }
        #                 self.changed_reactions.append(change)
                    
        #         if output_directory is not None:
        #             with open(os.path.join(output_directory, 'standardized_reactions.txt'), 'w') as out:
        #                 json.dump(self.changed_reactions, out, indent = 3)
                        
        #         total_reactions = 0
        #         for model in models:
        #             total_reactions += len(model.reactions)
                    
        #         warn(f'\nModelSEEDError: {missed_reactions}/{total_reactions} reactions were not captured by the ModelSEED modelreaction IDs.')
        # models[model_index] = model

                
        return models
    
    def correct_exchanges(self, models,                      # the collection of cobrakbase models that will be compared
                          conflicts_path_file: str = None,   # the metabolite conflicts are stored and organized
                          standardize: bool = False,         # standardize the model names and reactions to the ModelSEED Database
                          ): 
        if standardize:
            models = self.standardize_MSD(models)
            
        unique_names, established_mets, self.unknown_met_ids, self.changed_metabolites, self.changed_reactions = [], [], [], [], []
        self.unique_mets, self.met_conflicts = OrderedDict(), OrderedDict()
        for model in models:
            model_index = models.index(model)
            model_metabolites = {met.id:met for met in model.metabolites}
            for ex_rxn in model.exchanges:
                for met in ex_rxn.metabolites:
                    met_name = re.sub('_\w\d', '', met.name) #!!! The MODELCOMPOUNDS annotations should be matched preferentially to the metabolite names
                    if met.id not in self.unique_mets and met.id not in established_mets: #!!! The edge case where a few models use one ID and a few use another ID must be clarified
                        if met_name not in unique_names:
                            # identify the unique metabolite
                            self.unique_mets[met.id] = {
                                f'model{model_index}_id': met.id,
                                f'model{model_index}_met': met
                                }
                            unique_names.append(met_name)
                        else:
                            # describe the metabolite conflict between the ID and name
                            former_id = list(self.unique_mets.keys())[unique_names.index(met_name)]
                            former_model_index = list(self.unique_mets[former_id].keys())[0].split('_')[0].removeprefix('model')
                            if met.name not in self.met_conflicts:
                                self.met_conflicts[met_name] = {
                                        f'model{former_model_index}_id': former_id,
                                        f'model{former_model_index}_met': self.unique_mets[former_id][f'model{former_model_index}_met'],
                                        f'model{model_index}_id': met.id,
                                        f'model{model_index}_met': met
                                    }
                            else:
                                self.met_conflicts[met_name].update({
                                        f'model{model_index}_id': met.id,
                                        f'model{model_index}_met': met
                                    })
                            met = self._fix_met(met)
                    else:
                        former_name = unique_names[list(self.unique_mets.keys()).index(met.id)]
                        former_model_index = list(self.unique_mets[met.id].keys())[0].split('_')[0].removeprefix('model')
                        if met_name == former_name:
                            # remove the metabolite that is no longer unique
                            del unique_names[list(self.unique_mets.keys()).index(met.id)]
                            self.unique_mets.pop(met.id)    
                            established_mets.append(met.id)
                        else:
                            # describe the conflicting metabolite names
                            if met.id not in self.met_conflicts:
                                self.met_conflicts[met.id] = {
                                        f'model{former_model_index}_name': former_name,
                                        f'model{former_model_index}_met': self.unique_mets[former_id][f'model{former_model_index}_met'],
                                        f'model{model_index}_name': met.name,
                                        f'model{model_index}_met': met
                                    }
                            else:
                                if f'model{model_index}_name' not in self.met_conflicts[met.id]:
                                    self.met_conflicts[met.id].update({
                                            f'model{model_index}_name': met.name,
                                            f'model{model_index}_met': met
                                        })
                                else:
                                    iteration = 0
                                    while f'model{model_index}_{iteration}_name' in self.met_conflicts[met.id]:
                                        iteration += 1
                                        
                                    self.met_conflicts[met.id].update({
                                            f'model{model_index}_{iteration}_name': met.name,
                                            f'model{model_index}_{iteration}_met': met
                                        })
                            met = self._fix_met(met)
             
                    models[model_index] = model
                    
                # correct the reaction ID
                if re.sub('(_\w\d$)', '', ex_rxn.id).removeprefix('EX_') in model_metabolites:
                    suffix = re.search('(_\w\d$)', ex_rxn.id).group()
                    rxn_met = self._fix_met(re.sub('(_\w\d$)', '', ex_rxn.id).removeprefix('EX_'))
                    ex_rxn.id = 'EX_'+rxn_met.id+suffix
                            
        if conflicts_path_file is not None:
            export_met_conflicts = {}
            for met_id, content in self.met_conflicts.items():
                export_met_conflicts[met_id] = {}
                for key, val in content.items():
                    if '_met' not in key:
                        export_met_conflicts[met_id][key] = val
                    else:
                        export_met_conflicts[met_id][key.replace('_met','_formula')] = val.formula
            with open(conflicts_path_file, 'w') as out:
                json.dump(export_met_conflicts, out, indent = 3)

        return models
    
    def _fix_met(self,met):
        # correct the conflict
        base_name = ''.join(met.name.split('-')[1:]).capitalize()
        met_name = re.sub('_\w\d', '', met.name)
        if met.name in self.compound_names:
            met = self.__correct_met(met, met.name)
        elif met_name in self.compound_names:
            met = self.__correct_met(met, met_name)
        elif base_name in self.compound_names:
            met = self.__correct_met(met, base_name)
        else:
            self.unknown_met_ids.append(met.id)
            warn(f'ModelSEEDError: The metabolite {met.id} | {met.name} | {base_name} | {met_name} is not recognized by the ModelSEED Database')    
        return met
        
    def __correct_met(self, met, met_name, standardize = False):
        original_id = re.sub('(_\w\d)', '', met.id)
        if original_id != self.compound_names[met_name]:  # If the ID associated with the name deviates from that in the ModelSEED Database
            if self.compound_names[met_name] in met.model.metabolites:
                # replace the undesirable isomer in every instance, since it cannot be renamed
                for rxn in met.reactions:
                    original_reaction = rxn.reaction
                    reaction_dict = {}
                    for rxn_met in rxn.reactants+rxn.products:  # The REACTANTS+PRODUCTS may help resolve metabolites that are both more than METABOLITES
                        stoich = float(rxn.metabolites[rxn_met])
                        new_met = rxn_met
                        if rxn_met.id == met.id:
                            new_met = Metabolite(
                                id = self.compound_names[met_name], 
                                name = met_name,
                                formula = met.formula,
                                charge = met.charge,
                                compartment = met.compartment
                                )
                            
                        reaction_dict[new_met] = stoich
                        rxn.add_metabolites({rxn_met:0}, combine = False)
                        
                    rxn.add_metabolites(reaction_dict, combine = False)  
                    change = {
                            'original': {
                                    'reaction': original_reaction
                                },
                            'new': {
                                    'reaction': rxn.reaction
                                },
                            'justification': f'The non-ModelSEED ID {met.id} in this reaction ({rxn.id}) must be replaced.'
                        }
                    self.changed_reactions.append(change)
                    if self.printing:
                        print('\n')
                        pprint(change, sort_dicts=False)
            else:
                # rename the undesirable isomer
                met.id = self.compound_names[met_name]
                change = {
                        'original': {
                                'id': original_id,
                                'name': met.name
                            },
                        'new': {
                                'id': met.id,
                                'name': met_name
                            },
                        'justification': f'The {original_id} and {met.id} distinction is incompatible.'
                    }
                if 'cpd' not in original_id:
                    change['justification'] = f'The {original_id} ID is not a ModelSEED Database ID.'
                if standardize:
                    change['justification'] = f'The {original_id} and {met.id} metabolites were matched via their name.'
                self.changed_metabolites.append(change)
                if self.printing:
                    print('\n')
                    pprint(change, sort_dicts=False)

        return met
            