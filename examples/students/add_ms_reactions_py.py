import modelseedpy
import cobra

modelseed_path = '..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
bigg_model_path = '.\e_coli_core metabolism from BiGG.json'
model = cobra.io.load_json_model(bigg_model_path)
modelseed = modelseedpy.biochem.from_local(modelseed_path)

compartment_equivalents = {'0':'c'}

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
        compartment_number = metabolite[1]
        compartment_string = compartment_equivalents[compartment_number]        
        compound_abbreviation = '{}_{}'.format(compound['abbreviation'], compartment_string)

        metabolites_to_add[cobra.Metabolite(compound_abbreviation, name = compound['name'], compartment = compartment_string)] = stoich
            
    cobra_reaction.add_metabolites(metabolites_to_add)
    cobra_reaction.reaction
    model.add_reactions([cobra_reaction])


#add_ms_reaction(model, 'rxn00002', modelseed)
for reaction in model.reactions:
    print(reaction)