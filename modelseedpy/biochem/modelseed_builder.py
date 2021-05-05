
class ModelSeedBuilder:
    
    def __init__(self, fbamodel, modelseed):
        self.fbamodel = fbamodel
        self.modelseed = modelseed
        
    def configure_reaction(self, seed_reaction, direction, compartment_config, model_reaction_proteins):
        maxforflux = 1000
        maxrevflux = 1000
        if direction == '>':
            maxrevflux = 0
        elif direction == '<':
            maxforflux = 0

        if not seed_reaction['is_obsolete'] == 0:
            print('warning obsolete reaction', seed_id)

        model_reaction_reagents = configure_stoichiometry(
            seed_reaction, 
            compartment_config)

        compartment = cobrakbase.core.utils.get_reaction_compartment2(compartment_config)

        if len(compartment) > 1:
            compartment = compartment[0]

        modelreaction = {
            'aliases': [], 
            'dblinks': {}, 
            'direction': direction, 
            'edits': {}, 
            'gapfill_data': {}, 
            'id': "{}_{}".format(seed_reaction['id'], compartment),
            'maxforflux': maxforflux, 
            'maxrevflux': maxrevflux, 
            'modelReactionProteins': model_reaction_proteins, 
            'modelReactionReagents': model_reaction_reagents, 
            'modelcompartment_ref': '~/modelcompartments/id/c0', 
            'name': seed_reaction['name'], 
            'numerical_attributes': {}, 
            'probability': 0, 
            'protons': 0, 
            'reaction_ref': '~/template/reactions/id/{}_{}'.format(seed_reaction['id'], compartment), 
            'string_attributes': {}
        }
        
        return modelreaction