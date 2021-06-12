import modelseedpy
import cobra
import re
import os

modelseed_path = '..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
bigg_model_path = '.\e_coli_core metabolism from BiGG.json'
'''if os.path.exists(modelseed_path):
    print('yes')'''


model = cobra.io.load_json_model(bigg_model_path)
metabolites = []
for metabolite in model.metabolites:
    metabolites.append(metabolite.id)

modelseed = modelseedpy.biochem.from_local(modelseed_path)

def add_ms_reaction(model,rxn_id,modelseed):
    ''' Add a reaction with ModelSEED parameters to the FBA simulation
    "model" (COBRA object): The metabolic model that is defined by COBRApy
    "rxn_id" (Python object, string): The ModelSEED reaction id that will be added to the model
    "modelseed" (ModelSEED object): The ModelSEED database that describes the reaction and metabolites of the argument    
    '''
    modelseed_reaction = modelseed.get_seed_reaction(rxn_id)
    reaction_stoich = modelseed_reaction.cstoichiometry
    cobra_reaction = cobra.Reaction(rxn_id)
    cobra_reaction.name = modelseed_reaction.data['name']
        
    metabolites_to_add = {}
    for metabolite, stoich in reaction_stoich.items():
        id = metabolite[0]
        compound = modelseed.get_seed_compound(id).data
        abbreviation = compound['abbreviation']
        compartment = metabolite[1]
        
        if abbreviation not in metabolites:
            metabolites_to_add[cobra.Metabolite(abbreviation, name = compound['name'], compartment = compartment)] = stoich
            
    cobra_reaction.add_metabolites(metabolites_to_add)
    cobra_reaction.reaction
    model.add_reactions([cobra_reaction])


add_ms_reaction(model, 'rxn00002', modelseed)
for reaction in model.reactions:
    print(reaction)