from collections import OrderedDict
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
                    for name in names:
                        self.compound_names[name] = cpd['id']
            
    def _correct_met(self, met, metabolites_set, met_name):
        try:
            original_id = met.id
            if re.sub('(_\w\d)', '', original_id) != self.compound_names[met_name]:
                original_name = met.name
                met.id = self.compound_names[met_name]
                if self.printing:
                    print('\noriginal ID\t', original_id, '\t', original_name)
                    print('new ID\t\t', met.id, '\t', met.name)
        except: # the met.id already exists in the model
            for rxn in met.reactions:
                original_reaction = rxn.reaction
                if re.sub('(_\w\d)', '', original_reaction) != re.sub('(_\w\d)', '', rxn.reaction): #!!! What is the intention of this line?
                    rxn_mets = [rxn_met.id for rxn_met in rxn.metabolites]
                    if metabolites_set[met_name] not in rxn_mets:
                        rxn.add_metabolites({
                                    met: 0, metabolites_set[met_name]: rxn.metabolites[met]
                                }, combine = False)
                    else:   #!!! The designation of metabolites that are never corrected must be made, perhaps by returning None
                        warn(f'The {met.id} metabolite in the {rxn.id} reaction could not be corrected.')
                    if self.printing:
                        print('\noriginal met_rxn\t', original_reaction)
                        print('new met_rxn\t\t', rxn.reaction)
                    
        return met
        
    def standardize_MSD(self,model):
        for met in model.metabolites:
            if re.sub('(_\w\d)', '', met.id) in self.compounds_id_indexed:
                # standardize the metabolite names
                met.name = self.compounds_id_indexed[re.sub('(_\w\d)', '', met.id)]['name']
            else:
                warn(f'The {met.id} metabolite is not captured by the ModelSEED Database')

        for rxn in model.reactions:
            if re.sub('(_\w\d)', '', rxn.id) in self.reaction_ids:
                # standardize the reaction names
                rxn.name = self.reaction_ids[re.sub('(_\w\d)', '', rxn.id)]
                
                # standardize the reactions to the ModelSEED Database
                for index, rxn_id in enumerate(list(self.reaction_ids.keys())):
                    if rxn_id == re.sub('(_\w\d)', '', rxn.id):
                        reaction_dict = {}
                        for met in rxn.metabolites:
                            reaction_dict[met] = float(rxn.metabolites[met])
                        break
                
                for met in rxn.metabolites:
                    rxn.add_metabolites({met:0}, combine = False)
                rxn.add_metabolites(reaction_dict, combine = False)
            else:
                warn(f'The {rxn.id} reaction is not captured by the ModelSEED Database')
                for met in rxn.metabolites:
                    for name in self.compound_names:
                        if met.name == name:
                            met = self._correct_met(met, self.compound_names, name)
                            break
        
        return model
    
    def compare_models(self, model_1, model_2,       # arbitrary cobrakbase models
                       metabolites: bool = True,     # contrast metabolites (True) or reactions (False) between the models
                       standardize: bool = False     # standardize the model names and reactions to the ModelSEED Database
                       ): 
        misaligned = []
        if metabolites: # contrast metabolites 
            model2_ids = OrderedDict({met.id:re.sub('(_\w\d)', '', met.name) for met in model_2.metabolites})
            compared_met_counter = 0
            for met in model_1.metabolites:  #!!! develop dictionaries of ids and names to facilitate the comparison
                if met.id in model2_ids: 
                    met_name = re.sub('(_\w\d)', '', met.name)
                    if met_name != model2_ids[met.id]:
                        print(f'\nmisaligned {met.id} names\n', f'model1: {met_name}', f'model2: {model2_ids[met.id]}')
                        misaligned.append({
                                    'model_1': met_name,
                                    'model_2': model2_ids[met.id],
                                })
                        if self.printing:
                            print(met.name, model2_ids[met.id])
                        
                    compared_met_counter += 1
                
            if self.printing:
                print(f'\n\n{compared_met_counter} of the {len(model_1.metabolites)} model_1 metabolites and {len(model_2.metabolites)} model_2 metabolites are shared and were compared.')        
        else: # contrast reactions 
            model2_rxns = {rxn.id:rxn for rxn in model_2.reactions}
            compared_rxn_counter = 0
            for rxn in model_1.reactions:
                if rxn.id in model2_rxns:
                    if re.sub('(_\w\d)', '', rxn.name)!= re.sub('(_\w\d)', '', model2_rxns[rxn.id].name):
                        print(f'\nmisaligned {rxn.id} names\n', ' model1  ',rxn.name, ' model2  ',model2_rxns[rxn.id].name)
                        misaligned.append({
                                    'model_1': rxn.name,
                                    'model_2': model2_rxns[rxn.id].name,
                                })                            
                    elif rxn.reaction != model2_rxns[rxn.id].reaction:                            
                        print(f'\nmisaligned {rxn.id} reagents\n', 'model1: ',rxn.reaction, 'model2: ',model2_rxns[rxn.id].reaction)
                        misaligned.append({
                                    'model_1': rxn.reaction,
                                    'model_2': model2_rxns[rxn.id].reaction,
                                })

                    compared_rxn_counter += 1
                    
            if self.printing:
                print(f'''\n\n{compared_rxn_counter} of the {len(model_1.reactions)} model_1 reactions and {len(model_2.reactions)} model_2 reactions are shared and were compared.''')
                
        if standardize:
            model_1 = self.standardize_MSD(model_1)
            model_2 = self.standardize_MSD(model_2)
            
        return misaligned, model_1, model_2
        
    def exchanges(self,
                  model # cobrakbase model
                  ):                        
        # homogenize isomeric metabolites in exchange reactions
        unique_base_met_id = OrderedDict()
        unique_base_met_name = []
        for rxn in model.reactions:
            if 'EX_' in rxn.id:
                for met in rxn.metabolites:
                    if '-' in met.name:
                        base_name = ''.join(met.name.split('-')[1:])
                        if base_name not in unique_base_met_name:
                            unique_base_met_name.append(base_name)
                            unique_base_met_id[met.id] = met
                        else:
                            # replace isomers of different IDs with a standard isomer
                            base_name_index = unique_base_met_name.index(base_name)
                            base_name_id = list(unique_base_met_id.keys())[base_name_index]
                            if re.sub('(_\w\d)', '', met.id) != base_name_id:
                                met = self._correct_met(met, unique_base_met_id, base_name_id)
                    
        # standardize model metabolite IDs with the ModelSEED Database
        unknown_met_ids = []
        for met in model.metabolites:
            # correct non-standard IDs
            if 'cpd' not in met.id and len(met.reactions) >= 1:
                for name in self.compound_names:
                    if met.name == name:
                        met = self._correct_met(met, self.compound_names, name)
                        break
           
            if 'cpd' not in met.id and len(met.reactions) >= 1: 
                unknown_met_ids.append(met.id)
                warn(f'The metabolite {met.id} | {met.name} cannot not recognized by the ModelSEED Database')
        
        return model, unknown_met_ids