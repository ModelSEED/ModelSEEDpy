from collections import OrderedDict
from cobra.core.metabolite import Metabolite
from cobra.io.json import save_json_model
from zipfile import ZipFile, ZIP_LZMA
from warnings import warn
from pprint import pprint
import json, lzma, re, os

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
            self.compounds_cross_references, self.compound_names = OrderedDict(), OrderedDict()
            for cpd in self.compounds:
                self.compounds_cross_references[cpd['id']] = {}
                if cpd['aliases'] is not None:
                    for category in cpd['aliases']:
                        content = category.split(';')
                        if 'Name' in category:
                            content[0] = content[0].split(':')[0].strip()
                            names = [name.strip() for name in content]
                            names.append(cpd['name'])
                            for name in names:
                                if name not in self.compound_names:
                                    self.compound_names[name] = cpd['id']
                        else:
                            first = content[0].split(':')
                            db = first[0].strip()
                            content[0] = first[1]
                            self.compounds_cross_references[cpd['id']][db] = [x.strip() for x in content]
                        
    
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
        
    def standardize(self, models,                      # the collection of cobrakbase models that will be compared
                    metabolites: bool = True,          # specifies whether metabolites or reactions (FALSE) will be standardized
                    exchanges: bool = True,            # specifies whether only the exchange reaction will be standardized
                    conflicts_file_name: str = None,   # the metabolite conflicts are stored and organized, where None does not export
                    model_names: list = None,          # specifies the export names of the models
                    model_format: str = 'json',        # specifies to which format the model will be exported 
                    export_directory: str = None       # specifies the directory to which all of the content will be exported 
                    ):
        self.models = models
        self.unique_mets, self.met_conflicts = OrderedDict(), OrderedDict()
        self.unknown_met_ids, self.changed_metabolites, self.changed_reactions= [], [], []
        self.changed_ids_count = self.changed_rxn_count = 0
        for self.model_index, self.model in enumerate(self.models):
            # standardize metabolites
            if metabolites:
                if exchanges:
                    model_metabolites = [met.id for met in self.model.metabolites]
                    for ex_rxn in self.model.exchanges:
                        for met in ex_rxn.metabolites:
                            met, new_met_id, success = self._fix_met(met)
                            try:
                                ex_rxn.id = 'EX_'+met.id
                            except:
                                ex_rxn.id = 'EX_'+new_met_id
                            if 'cpd' not in met.id and success and new_met_id not in model_metabolites: 
                                self.unknown_met_ids.append(met.id)
                                warn(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                else:
                    for met in self.model.metabolites:
                        met, new_met_id, success = self._fix_met(met)
                        if 'cpd' not in met.id: 
                            self.unknown_met_ids.append(met.id)
                            warn(f'CodeError: The metabolite {met.id} | {met.name} was not corrected to a ModelSEED metabolite.')
                
                if conflicts_file_name is not None:
                    self._export({'metabolite_changes':self.changed_metabolites, 'reaction_changes':self.changed_reactions}, 
                        conflicts_file_name, model_names, model_format, export_directory
                        )
                        
            # standardize reactions
            # else:  #!!! The modelreactions appear to be incorrect
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
            #             reaction_dict = self._parse_modelReactionReagents(
            #                 modelreactions_ids[re.sub('(_\w\d$)', '', rxn.id)]['modelReactionReagents'], model_metabolites
            #                 )
            #         elif rxn.id in modelreactions_ids:
            #             reaction_dict = self._parse_modelReactionReagents(
            #                 modelreactions_ids[rxn.id]['modelReactionReagents'], model_metabolites
            #                 )
            #         else:
            #             warn(f'ModelSEEDError: The reaction ID {rxn.id} is not captured by the modelreactions.')
                    
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
            #         if rxn.id not in self.reaction_ids: 
            #             missed_reactions += 1
            #             # warn(f'ModelSEEDError: The {rxn.id} | {rxn.name} reaction is not recognized by the ModelSEED Database')
            
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
            #             self.changed_reactions.append(change)
                    
            #     if export_directory is not None:
            #         with open(os.path.join(export_directory, 'standardized_reactions.txt'), 'w') as out:
            #             json.dump(self.changed_reactions, out, indent = 3)
                        
            #     total_reactions = 0
            #     for model in models:
            #         total_reactions += len(model.reactions)
                    
            #     warn(f'\nModelSEEDError: {missed_reactions}/{total_reactions} reactions were not captured by the ModelSEED modelreaction IDs.')
        
            self.models[self.model_index] = self.model
        print(f'\n\n{self.changed_rxn_count} reactions were substituted and {self.changed_ids_count} metabolite IDs were redefined.')
        return self.models
       
    def align_exchanges(self, models,                        # the collection of cobrakbase models that will be compared
                          standardize: bool = False,         # standardize the model names and reactions to the ModelSEED Database
                          conflicts_file_name: str = None,   # the metabolite conflicts are stored and organized, where None does not the conflicts
                          model_names: list = None,          # specifies the name of the exported model, where None does not export the models
                          model_format: str = 'json',        # specifies to which format the model will be exported 
                          export_directory: str = None       # specifies the directory to which all of the content will be exported 
                          ): 
        self.models = models
        self.changed_ids_count = self.changed_rxn_count = 0
        if standardize:
            self.standardize_MSD(self.models)
            
        unique_names, established_mets, self.unknown_met_ids, self.changed_metabolites, self.changed_reactions = [], [], [], [], []
        self.unique_mets, self.met_conflicts = OrderedDict(), OrderedDict()
        for self.model_index, self.model in enumerate(self.models):
            model_metabolites = {met.id:met for met in self.model.metabolites}
            for ex_rxn in self.model.exchanges:
                for met in ex_rxn.metabolites:
                    met_name = re.sub('_\w\d$', '', met.name) 
                    if met.id not in self.unique_mets and met.id not in established_mets: 
                        if met_name not in unique_names:
                            # identify the unique metabolite
                            self.unique_mets[met.id] = {
                                f'model{self.model_index}_id': met.id,
                                f'model{self.model_index}_met': met
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
                                        f'model{self.model_index}_id': met.id,
                                        f'model{self.model_index}_met': met
                                    }
                            else:
                                self.met_conflicts[met_name].update({
                                        f'model{self.model_index}_id': met.id,
                                        f'model{self.model_index}_met': met
                                    })
                            met, new_met_id, success = self._fix_met(met)
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
                                        f'model{self.model_index}_name': met.name,
                                        f'model{self.model_index}_met': met
                                    }
                            else:
                                if f'model{self.model_index}_name' not in self.met_conflicts[met.id]:
                                    self.met_conflicts[met.id].update({
                                            f'model{self.model_index}_name': met.name,
                                            f'model{self.model_index}_met': met
                                        })
                                else:
                                    iteration = 0
                                    while f'model{self.model_index}_{iteration}_name' in self.met_conflicts[met.id]:
                                        iteration += 1
                                        
                                    self.met_conflicts[met.id].update({
                                            f'model{self.model_index}_{iteration}_name': met.name,
                                            f'model{self.model_index}_{iteration}_met': met
                                        })
                            met, new_met_id, success = self._fix_met(met)
             
                    self.models[self.model_index] = self.model
                    
                # correct the reaction ID
                if re.sub('(_\w\d$)', '', ex_rxn.id).removeprefix('EX_') in model_metabolites:
                    suffix = re.search('(_\w\d$)', ex_rxn.id).group()
                    rxn_met, new_met_id, success = self._fix_met(re.sub('(_\w\d$)', '', ex_rxn.id).removeprefix('EX_'))
                    ex_rxn.id = 'EX_'+new_met_id+suffix

        if conflicts_file_name:
            export_met_conflicts = {}
            for met_id, content in self.met_conflicts.items():
                export_met_conflicts[met_id] = {}
                for key, val in content.items():
                    if '_met' not in key:
                        export_met_conflicts[met_id][key] = val
                    else:
                        export_met_conflicts[met_id][key.replace('_met','_formula')] = val.formula
                        
            self._export(export_met_conflicts, conflicts_file_name, model_names, model_format, export_directory)

        print(f'\n\n{self.changed_rxn_count} exchange reactions were substituted and {self.changed_ids_count} exchange metabolite IDs were redefined.')
        return self.models
    
    def _fix_met(self,met):
        # correct the conflict
        base_name = ''.join(met.name.split('-')[1:]).capitalize()
        met_name = re.sub('_\w\d$', '', met.name)
        new_met_id = met.id
        success = True
        if met.name in self.compound_names:
            met, new_met_id = self.__correct_met(met, met.name)
        elif met.name.capitalize() in self.compound_names:
            met, new_met_id = self.__correct_met(met, met.name.capitalize())
        elif met_name in self.compound_names:
            met, new_met_id = self.__correct_met(met, met_name)
        elif met_name.capitalize() in self.compound_names:
            met, new_met_id = self.__correct_met(met, met_name.capitalize())
        elif base_name in self.compound_names and base_name != '':
            met, new_met_id = self.__correct_met(met, base_name)
        else:
            self.unknown_met_ids.append(met.id)
            success = False
            warn(f'ModelSEEDError: The metabolite ({" | ".join([x for x in [met.id, met.name, base_name, met_name] if x != ""])}) is not recognized by the ModelSEED Database')    
        return met, new_met_id, success
    
    def _export(self, conflicts,             # the conflicts dictionary that will be exported
                conflicts_file_name,         # the metabolite conflicts are stored and organized, where None does not the conflicts
                model_names,                 # specifies the name of the exported model, where None does not export the models
                model_format,                # specifies to which format the model will be exported 
                export_directory             # specifies the directory to which all of the content will be exported 
                ):
        if export_directory is None:
            export_directory = os.getcwd()
                    
        file_paths = []
        if conflicts_file_name is not None:
            path = os.path.join(export_directory,conflicts_file_name)
            file_paths.append(os.path.relpath(path, export_directory))
            with open(path, 'w') as out:
                json.dump(conflicts, out, indent = 3)
        if model_names is not None:
            for index, model in enumerate(self.models):
                path = os.path.join(export_directory,f'{model_names[index]}.{model_format}')
                file_paths.append(os.path.relpath(path, export_directory))
                save_json_model(model, path)
        with ZipFile('_'.join(model_names[:4])+'.zip', 'w', compression = ZIP_LZMA) as zip:
            for file in file_paths:
                zip.write(file)
                os.remove(file)
        
    def __correct_met(self, met, met_name, standardize = False):
        def check_cross_references(met, general_met):
            for db in self.compounds_cross_references[general_met]:
                for cross_ref in self.compounds_cross_references[general_met][db]:
                    if cross_ref in self.compounds_cross_references[self.compound_names[met_name]][db]:
                        match = True
                        break
                if match:
                    break
            return match, db
        
        original_id = new_met_id = met.id
        compartment = re.search('(_\w\d$)', met.id).group()
        if met.id.removesuffix(compartment) != self.compound_names[met_name]:  # If the ID associated with the name deviates from that in the ModelSEED Database
            new_met_id = self.compound_names[met_name]+compartment
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
                                warn(f'CodeError: The metabolite {new_met_id} replacement for {met.id} already exists in the reaction {rxn.id}, thus the reaction cannot be updated.')
                                break
                            
                            # affirm the match with cross-references, where it is possible for ModelSEED compounds
                            general_met = re.sub("(_\w\d$)", "", met.id)
                            if 'cpd' in met.id and self.compounds_cross_references[general_met] != {}:
                                match, db = check_cross_references(met, general_met)
                                if not match:
                                    warn(f'ModelSEEDError: The old metabolite {met.id} cross-references ({self.compounds_cross_references[general_met]}) do not overlap with those ({self.compounds_cross_references[self.compound_names[met_name]]}) of the new metabolite {new_met_id}.')
                                    
                            # remove duplicate exchange reaction
                            if 'EX_' in rxn.id and 'EX_'+new_met_id in [ex_rxn.id for ex_rxn in self.model.exchanges]:
                                change = {
                                        'original': {
                                                'reaction': original_reaction
                                            },
                                        'new': {
                                                'reaction': None
                                            },
                                        'justification': f'A {new_met_id} exchange reaction already exists in model {self.model_index}, thus this duplicative exchange reaction ({rxn.id}) is deleted.'
                                    }
                                if match:
                                    change['justification'] += f' The ID match was verified with {db} cross-references.'
                                self.model.remove_reactions([rxn.id])
                                self.changed_reactions.append(change)
                                if self.printing:
                                    print('\n')
                                    pprint(change, sort_dicts=False)
                                self.changed_rxn_count += 1
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
                                'justification': f'The {new_met_id} replacement for {met.id} already exists in model {self.model_index}, so each reaction (here {rxn.id}) must be updated.'
                            }
                        if match:
                            change['justification'] += f' The ID match was verified with {db} cross-references.'
                        self.changed_reactions.append(change)
                        if self.printing:
                            print('\n')
                            pprint(change, sort_dicts=False)
                            
                        self.changed_rxn_count += 1
                    else:
                        warn(f'CodeError: The reaction {reaction_dict} | {new_reactants} {new_products} possesses a different number of reagents than the original reaction {original_reaction} | {num_reactants} {num_products}, and is skipped.')
            else:
                # affirm the match with cross-references, where it is possible for ModelSEED compounds
                match = False
                general_met = re.sub("(_\w\d$)", "", met.id)
                if 'cpd' in met.id and self.compounds_cross_references[general_met] != {}:
                    match, db = check_cross_references(met, general_met)
                    if not match:
                        warn(f'ModelSEEDError: The old metabolite {met.id} cross-references ({self.compounds_cross_references[general_met]}) do not overlap with those ({self.compounds_cross_references[self.compound_names[met_name]]}) of the new metabolite {new_met_id}.')
                            
                # rename the undesirable isomer
                met.id = self.compound_names[met_name]+compartment
                change = {
                        'original': {
                                'id': original_id,
                                'name': met.name
                            },
                        'new': {
                                'id': met.id,
                                'name': met_name+compartment
                            },
                        'justification': f'The {original_id} and {met.id} distinction in {self.model_index} is incompatible.'
                    }
                if 'cpd' not in original_id:
                    change['justification'] = f'The {original_id} ID is not a ModelSEED Database ID.'
                if standardize:
                    change['justification'] = f'The {original_id} and {met.id} metabolites were matched via their name.'
                if match:
                    change['justification'] += f' The ID match was verified with {db} cross-references.'
                self.changed_metabolites.append(change)
                if self.printing:
                    print('\n')
                    pprint(change, sort_dicts=False)
                self.changed_ids_count += 1

        return met, new_met_id
            