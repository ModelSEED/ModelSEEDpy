from numpy import negative
from warnings import warn
from pprint import pprint
import pandas
import json, re, os

# define the environment path 
local_cobrakbase_path = os.path.join('C:', 'Users', 'Andrew Freiburger','Documents','Argonne','cobrakbase')
os.environ["HOME"] = local_cobrakbase_path

class align_models:
    def __init__(self,
                 modelseed_db_path: str, # the local path to the ModelSEEDDatabase repository
                 printing = True         # specifies whether results are printed
                 ):       
        self.printing = printing
        with open(os.path.join(modelseed_db_path, 'Biochemistry', 'reactions.json')) as rxns:
            self.reactions = json.load(rxns)
            
    def _parse_ms_rxn_string(self,reaction_string):
        # parse the reaction string
        if '<=>' in reaction_string:
            compounds = reaction_string.split('<=>')
        elif '=>' in reaction_string:
            compounds = reaction_string.split('=>')
        elif '<=' in reaction_string:
            compounds = reaction_string.split('<=')
        else:
            warn(f'The reaction string {reaction_string} has an unexpected reagent delimiter.')
        reactant, product = compounds[0], compounds[1]
        reactants = [x.strip() for x in reactant.split('+')]
        products = [x.strip() for x in product.split('+')]
        reactant_met = [x.split(' ') for x in reactants]
        product_met = [x.split(' ') for x in products]
        
        # assemble a dictionary for the reaction
        reaction_dict = {}
        for met in reactant_met:
            stoich = float(re.search('(\d+)', met[0]).group())
            met[1] = met[1].replace('[0]', '_c0')
            met[1] = met[1].replace('[1]', '_c0')
            reaction_dict[met[1]] = negative(stoich)
        for met in product_met:
            stoich = float(re.search('(\d+)', met[0]).group())
            met[1] = met[1].replace('[0]', '_c0')
            met[1] = met[1].replace('[1]', '_c0')
            reaction_dict[met[1]] = stoich
            
        return reaction_dict
    
    def _parse_modelReactionReagents(self,modelReactionReagents):
        # parse the reaction dictionary
        rxn_dict = {}
        for cpd in modelReactionReagents:
            cpd_id = re.search('(cpd\d{5})', cpd['modelcompound_ref']).group()
            stoich = float(cpd['coefficient'])
            rxn_dict[cpd_id] = stoich
            
        return rxn_dict
        
    def _parse_models(self,
                     model_1, # an arbitrary cobrakbase model
                     model_2, # an arbitrary cobrakbase model
                     metabolites = True, # contrast the metabolic termology between the models
                     ):
        if metabolites:
            misaligned_metabolites = []
            compared_met_counter = 0
            for met_id in model_1.metabolites:
                if met_id in model_2.metabolites:
                    # define model1 specifications
                    model1_met_index = model_1.metabolites.index(met_id)
                    model1_met_dict = model_1.modelcompounds[model1_met_index]
                    model1_met_name = re.sub('(_\w\d)$', '', model1_met_dict['name'])
                    
                    # define model2 specifications
                    model2_met_index = model_2.metabolites.index(str(met_id))
                    model2_met_dict = model_2.modelcompounds[model2_met_index]
                    model2_met_name = re.sub('(_\w\d)$', '', model2_met_dict['name'])    
                    
                    if model1_met_name != model2_met_name:
                        print('\nmisaligned met names', 'model1:',model1_rxn_index, 'model2:',model2_rxn_index)
                        misaligned_metabolites.append({
                                    'model_1': model1_met_name,
                                    'model_2': model2_met_name,
                                })
                        if self.printing:
                            print(model1_met_name, model2_met_name)
                        
                    compared_met_counter += 1
                    
            if self.printing:
                print(f'''\n\n{compared_met_counter}/{len(model_1.metabolites)} model_1 metabolites are shared with model_2 and were compared.''')
            
            return misaligned_metabolites, model_1, model_2

        else:
#            pprint(dir(self.model1))
#            with open('model_reactions.json', 'w') as reactions:
#                json.dump(self.model1.modelreactions, reactions, indent = 3)
#            with open('reactions.json', 'w') as rxns:
#                rxns.write('\n'.join([str(x) for x in self.model1.reactions]))
            model1_rxns = [rxn['id'] for rxn in model_1.modelreactions]
            model2_rxns = misaligned_reactions = model1_exchanges = []
            model2_ex_mets = set()
            model2_ex = {}
            for rxn in model_2.modelreactions:
                 model2_rxns.append(rxn['id'])
                 if re.search('(R-|EX_)', rxn['id']):
                     reaction_dict = self._parse_ms_rxn_string(rxn['equation'])
                     model2_ex_mets.add()
                     for met in model2_ex_metabolites:
                         model2_ex[rxn['id']] = {
                                     met['id']: met['name']
                                 }
            compared_reaction_counter = 0
            for rxn in model_1.modelreactions:
                if re.search('(R-|EX_)', rxn['id']): # exchange reactions
                    model1_exchanges.append(rxn)   # parse the reaction via a function to acquire the set of metabolites
                if rxn['id'] in model2_rxns: # identical cytoplasmic reactions
                    # define model1 specifications
                    model1_rxn_index = model1_rxns.index(rxn['id'])
                    model1_rxn_dict = model_1.modelreactions[rxn]
                    model1_rxn_id = re.sub('(_\w\d)$', '', model1_rxn_dict['id'])
                    model1_rxn_name = re.sub('(_\w\d)$', '', model1_rxn_dict['name'])
                    model1_rxn_reagents = model1_rxn_dict['modelReactionReagents']
                    
                    # define model2 specifications
                    model2_rxn_index = model2_rxns.index(rxn['id'])
                    model2_rxn_dict = model_2.modelreactions[model2_rxn_index]
                    model2_rxn_id = re.sub('(_\w\d)$', '', model2_rxn_dict['id'])  
                    model2_rxn_name = re.sub('(_\w\d)$', '', model2_rxn_dict['name'])
                    model2_rxn_reagents = model2_rxn_dict['modelReactionReagents']
                    
                    # assess the alignment between the reaction descriptions
                    if model1_rxn_name != model2_rxn_name:
                        print(f'\nmisaligned names {model2_rxn_id}', ' model1  ',model1_rxn_name, ' model2  ',model2_rxn_name)
                        misaligned_reactions.append({
                                    'model_1': model1_rxn_name,
                                    'model_2': model2_rxn_name,
                                })
                        # determine the correct reaction name according to the ModelSEEDDatabase
                        for rxn in self.reactions:
                            if model2_rxn_id == rxn['id']:
                                reaction_name = rxn['name']
                                break
                            
                        # correct the erroneous model
                        if model1_rxn_name != reaction_name:
                            model_1.modelreactions[model1_rxn_index]['name'] = reaction_name
                        elif model2_rxn_name != reaction_name:
                            model_2.modelreactions[model2_rxn_index]['name'] = reaction_name
                    elif model1_rxn_reagents != model2_rxn_reagents:                            
                        # determine the correct reaction according to the ModelSEEDDatabase
                        for rxn in self.reactions:
                            if model2_rxn_id == rxn['id']:
                                reaction_dict = self._parse_ms_rxn_string(rxn['equation'])
                                break
                            elif model1_rxn_id == rxn['id']:
                                reaction_dict = self._parse_ms_rxn_string(rxn['equation'])
                                break
                        
                        # standardize the reaction dictionaries of the questioned reactions
                        model1_rxn_dict = self._parse_modelReactionReagents(model1_rxn_reagents)
                        model2_rxn_dict = self._parse_modelReactionReagents(model1_rxn_reagents)
                        if model1_rxn_dict != model2_rxn_dict:
                            print('\nmisaligned reagents', 'model1:',model1_rxn_id, 'model2:',model2_rxn_id)
                            misaligned_reactions.append({
                                        'model_1': model1_rxn_dict,
                                        'model_2': model2_rxn_dict,
                                    })
                            if self.printing:
                                print('\nmodel1:')
                                pprint(model1_rxn_dict)
                                print('\nmodel2:')
                                pprint(model2_rxn_dict)
                            
                            # correct the erroneous model 
                            modelReactionReagents_dict = {}
                            if model2_rxn_dict == reaction_dict:
                                for cpd, stoich in reaction_dict.items():
                                    modelReactionReagents_dict['modelcompound_ref'] = cpd
                                    modelReactionReagents_dict['coefficient'] = stoich
                                model_1.modelreactions[model1_rxn_index] = modelReactionReagents_dict
                            elif model1_rxn_dict == reaction_dict:
                                for cpd, stoich in reaction_dict.items():
                                    modelReactionReagents_dict['modelcompound_ref'] = cpd
                                    modelReactionReagents_dict['coefficient'] = stoich
                                model_2.modelreactions[model1_rxn_index] = modelReactionReagents_dict
                            else:
                                warn(f'The {rxn_id} reactions of model_1 {model1_rxn_dict} and model_2 {model2_rxn_dict} nether align with the ModelSEED Database reaction {reaction_dict} nor align with themselves.')
                        
                    compared_reaction_counter += 1
                    
            if self.printing:
                print(f'''\n\n{compared_reaction_counter}/{len(model_1.reactions)} model_1 reactions are shared with model_2 and were compared.''')
                
            # contrast
            
            return misaligned_reactions, model_1, model_2

        
    def align(self,model1, model2):
        self.model1, self.model2 = model1, model2
        
        # compare metabolites between the models
        misaligned_met_model1_in_model2, self.model1, self.model2 = self._parse_models(model1, model2)
        print(f'Misaligned metabolites: {len(misaligned_met_model1_in_model2)}\n', misaligned_met_model1_in_model2)
        
        # compare reactions between the models
        misaligned_rxn_model1_in_model2, self.model1, self.model2 = self._parse_models(model1, model2, metabolites = False)
        print(f'Misaligned reactions: {len(misaligned_rxn_model1_in_model2)}\n', misaligned_rxn_model1_in_model2)
               

# import the models
import cobrakbase
token = 'WE6CHYRDTJSGOHFIDGPE7WYFT6PRPXJL'
kbase_api = cobrakbase.KBaseAPI(token)
model1 = kbase_api.get_from_ws("iML1515",76994)
model2 = kbase_api.get_from_ws("iSB1139.kb.gf",30650)
#print(model1.modelreactions[8])
#print(model2.modelreactions[591])

#from collections import Counter
#counts = dict(Counter(list(model2.metabolites)))
#duplicates = {key:value for key, value in counts.items() if value > 1}
#print(duplicates)

alignment = align_models(modelseed_db_path = os.path.join('..', 'ModelSEEDDatabase'))
alignment.align(model1, model2)