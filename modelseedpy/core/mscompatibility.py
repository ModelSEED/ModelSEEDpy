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
        
    def standardize_MSD(self, models, output_path = None):
        self.unique_mets, self.met_conflicts = OrderedDict(), OrderedDict()
        self.unknown_met_ids, self.changed_metabolites, self.changed_reactions= [], [], []
        for model in models:
            model_index = models.index(model)
            
            # standardize metabolites
            for met in model.metabolites:
                base_name = ''.join(met.name.split('-')[1:]).capitalize()
                met_name = re.sub('_\w\d', '', met.name)
                if met.name in self.compound_names:
                    met = self._correct_met(met, met.name)
                elif met_name in self.compound_names:
                    met = self._correct_met(met, met_name)
                elif base_name in self.compound_names:
                    met = self._correct_met(met, base_name)
                else:
                    self.unknown_met_ids.append(met.id)
                    warn(f'ModelSEEDError: The metabolite {met.id} | {met.name} is not recognized by the ModelSEED Database')
                        
                # ensure that the metabolite ID is conventional in the ModelSEED Database
                if 'cpd' not in met.id: 
                    self.unknown_met_ids.append(met.id)
                    warn(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                        
            # standardize reactions
            missed_reactions = 0
            for rxn in model.reactions:
                if re.sub('(_\w\d)', '', rxn.id) in self.reaction_ids:
                    rxn.name = self.reaction_ids[re.sub('(_\w\d)', '', rxn.id)]
                    reaction_dict = {}
                    for met in rxn.metabolites:
                        new_met = self._correct_met(met, re.sub('(_\w\d)', '', met.name))
                        reaction_dict[new_met] = float(rxn.metabolites[met])
                    for met in rxn.metabolites:
                        rxn.add_metabolites({met:0}, combine = False)
                    rxn.add_metabolites(reaction_dict, combine = False)
                else:
                    missed_reactions += 1
                    # warn(f'ModelSEEDError: The {rxn.id} | {rxn.name} reaction is not recognized by the ModelSEED Database')
                    
            print(model_index, missed_reactions)
            models[model_index] = model
            
        if output_path is not None:
            with open(output_path, 'w') as out:
                json.dump(self.changed_metabolites, out, indent = 3)
                
        return models
    
    def correct_exchanges(self, models,                      # the collection of cobrakbase models that will be compared
                          conflicts_path_file: str = None,   # the metabolite conflicts are stored and organized
                          standardize: bool = False,         # standardize the model names and reactions to the ModelSEED Database
                          ): 
        def fix_met(met):
            # correct the conflict
            base_name = ''.join(met.name.split('-')[1:]).capitalize()
            if met.name in self.compound_names:
                met = self._correct_met(met, met.name)
            elif base_name in self.compound_names:
                met = self._correct_met(met, base_name)
            else:
                self.unknown_met_ids.append(met.id)
                warn(f'ModelSEEDError: The metabolite {met.id} | {met.name} | {base_name} is not recognized by the ModelSEED Database')    
            return met
        
        if standardize:
            models = self.standardize_MSD(models)
            
        unique_names, established_mets, self.unknown_met_ids, self.changed_metabolites, self.changed_reactions = [], [], [], [], []
        self.unique_mets, self.met_conflicts = OrderedDict(), OrderedDict()
        for model in models:
            model_index = models.index(model)
            for ex_rxn in model.exchanges:
                for met in ex_rxn.metabolites:
                    met_name = re.sub('_\w\d', '', met.name)
                    if met.id not in self.unique_mets and met.id not in established_mets:
                        if met_name not in unique_names:
                            self.unique_mets[met.id] = {
                                f'model{model_index}_id': met.id,
                                f'model{model_index}_met': met
                                }
                            unique_names.append(met_name)
                        else:
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
                            # correct the conflict
                            met = fix_met(met)
    
                    else:
                        former_name = unique_names[list(self.unique_mets.keys()).index(met.id)]
                        former_model_index = list(self.unique_mets[met.id].keys())[0].split('_')[0].removeprefix('model')
                        if met_name == former_name:
                            # remove the metabolite that is no longer unique
                            del unique_names[list(self.unique_mets.keys()).index(met.id)]
                            self.unique_mets.pop(met.id)    
                            established_mets.append(met.id)
                        else:
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
                            # correct the conflict
                            met = fix_met(met)
             
                    models[model_index] = model
                            
        if conflicts_path_file is not None:
            export_met_conflicts = {}
            for met_id, content in self.met_conflicts.items():
                print(content)
                export_met_conflicts[met_id] = {}
                for key, val in content.items():
                    if '_met' not in key:
                        export_met_conflicts[met_id][key] = val
                    else:
                        export_met_conflicts[met_id][key.replace('_met', '_met_formula')] = val.formula
            with open(conflicts_path_file, 'w') as out:
                json.dump(export_met_conflicts, out, indent = 3)

        return models
        
    def _correct_met(self, met, met_name, standardize = False):
        original_id = re.sub('(_\w\d)', '', met.id)
        if original_id != self.compound_names[met_name]:  # If the ID associated with the name deviates from that in the ModelSEED Database
            # replace the undesirable isomer
            if self.compound_names[met_name] in met.model.metabolites:
                for rxn in met.reactions:
                    original_reaction = rxn.reaction
                    for rxn_met in rxn.reactants+rxn.products:
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
                        
                        rxn.add_metabolites({rxn_met:0}, combine = False)
                        rxn.add_metabolites({new_met:stoich}, combine = False)  
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
                    change['justification'] = 'The {original_id} ID is not a ModelSEED Database ID.'
                if standardize:
                    change['justification'] = 'The {original_id} and {met.id} metabolites were matched via their name.'
                self.changed_metabolites.append(change)
                if self.printing:
                    print('\n')
                    pprint(change, sort_dicts=False)

                return met