from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper
from cobra import Model, Reaction, Metabolite
import re


def build_from_species_models(org_models, model_id=None, name=None, names=None, abundances=None, cobra_model=False, standardize=False):
    """Merges the input list of single species metabolic models into a community metabolic model

    Parameters
    ----------
    org_models : list<Cobra.Model> to be merged into a community model
    model_id : string specifying community model ID
    name : string specifying community model name
    names : list<string>  human-readable names for models being merged
    abundances : dict<string,float> relative abundances for input models in community model
    cobra_model : bool for whether the raw COBRA model is returned
    standardize: bool for whether the exchanges of each member model will be standardized (True) or just aligned.

    Returns
    -------
    Cobra.Model for the desired Community

    Raises
    ------
    """
    # construct the new model
    names = names or []
    if standardize:
        models = MSCompatibility.standardize(org_models, conflicts_file_name='exchanges_conflicts.json', model_names=names)
    else:
        models = MSCompatibility.align_exchanges(org_models, 'exchanges_conflicts.json', names)
    newmodel = Model(model_id, name)
    biomass_compounds, biomass_indices = [], []
    biomass_index = minimal_biomass_index = 2
    new_metabolites, new_reactions = set(), set()
    for model_index, org_model in enumerate(models):
        model = org_model.copy()
        model_reaction_ids = [rxn.id for rxn in model.reactions]
        model_index += 1
        # print([rxn.id for rxn in model.reactions if "bio" in rxn.id])
        # print(model_index, model.id)
        # Rename metabolites
        for met in model.metabolites:
            # Renaming compartments
            output = MSModelUtil.parse_id(met)
            if output is None:
                index = 0 if met.compartment[0] == "e" else model_index
                met.compartment = met.compartment[0] + str(index)
                if "_" in met.id:
                    print(met.id)
                    met.id = met.id.split("_")[:-1] + met.compartment
            else:
                name, compartment, index = output
                index = 0 if compartment == "e" else model_index
                if index == "":
                    met.id += str(index)
                    met.compartment += str(index)
                elif compartment == "e":
                    met.compartment = "e0"
                else:
                    met.compartment = compartment + str(index)
                    met.id = name + "_" + met.compartment
            new_metabolites.add(met)
            if "cpd11416_c" in met.id:
                print(met.id, model.id)
                biomass_compounds.append(met)
        # Rename reactions
        for rxn in model.reactions:  # !!! all reactions should have a non-zero compartment index
            if rxn.id[0:3] != "EX_":
                if re.search('^(bio)(\d+)$', rxn.id):
                    index = int(rxn.id.removeprefix('bio'))
                    if index not in biomass_indices and index >= minimal_biomass_index:
                        biomass_indices.append(index)
                        print(rxn.id, '2')
                    else:  # biomass indices can be decoupled from the respective reaction indices of the same model
                        rxn.id = "bio" + str(biomass_index)
                        if rxn.id not in model_reaction_ids:
                            print(rxn.id, '1')
                            biomass_indices.append(biomass_index)
                        else:
                            index = minimal_biomass_index
                            rxn.id = "bio" + str(index)
                            while rxn.id not in model_reaction_ids and index not in biomass_indices:
                                index += 1
                                rxn.id = "bio" + str(index)
                            biomass_indices.append(index)
                            print(rxn.id, '3')
                    biomass_index += 1
                else:
                    output = MSModelUtil.parse_id(rxn)
                    if output is None:
                        if "e" not in rxn.compartment.id and not rxn.compartment.id[-1].isnumeric():
                            rxn.id += str(model_index)
                    else:
                        name, compartment, index = output
                        if compartment != "e":
                            rxn.name = name + "_" + compartment + str(model_index)
                            if index == "":
                                rxn.id += str(model_index)
            new_reactions.add(rxn)
        print(biomass_indices)
    # adds only unique reactions and metabolites to the community model
    newmodel.add_reactions(FBAHelper.filter_cobra_set(new_reactions))
    newmodel.add_metabolites(FBAHelper.filter_cobra_set(new_metabolites))

    # Create community biomass
    comm_biomass = Metabolite("cpd11416_c0", None, "Community biomass", 0, "c0")
    metabolites = {comm_biomass: 1}
    metabolites.update({cpd: -1 / len(biomass_compounds) for cpd in biomass_compounds})
    comm_biorxn = Reaction(id="bio1", name="bio1", lower_bound=0, upper_bound=100)
    comm_biorxn.add_metabolites(metabolites)
    newmodel.add_reactions([comm_biorxn])

    # define the model objective
    FBAHelper.add_objective(newmodel, comm_biorxn.flux_expression)

    # create a biomass sink reaction
    newutl = MSModelUtil(newmodel)
    newutl.add_exchanges_for_metabolites([comm_biomass], 0, 100, 'SK_')
    if cobra_model:
        return newmodel
    return newmodel, names, abundances


class CommHelper:

    @staticmethod
    def placeholder():
        pass
