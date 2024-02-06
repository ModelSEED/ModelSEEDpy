# -*- coding: utf-8 -*-
import logging
import re
import time
import json
import sys
import pandas as pd
import cobra
from cobra import Model, Reaction, Metabolite
from optlang.symbolics import Zero
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem
from modelseedpy.core.fbahelper import FBAHelper
from multiprocessing import Value

# from builtins import None

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO


class MSModelUtil:
    mdlutls = {}

    @staticmethod
    def metabolite_msid(metabolite):
        if re.search("^(cpd\d+)", metabolite.id):
            m = re.search("^(cpd\d+)", metabolite.id)
            return m[1]
        for anno in metabolite.annotation:
            if isinstance(metabolite.annotation[anno], list):
                for item in metabolite.annotation[anno]:
                    if re.search("^(cpd\d+)", item):
                        m = re.search("^(cpd\d+)", item)
                        return m[1]
            elif re.search("^(cpd\d+)", metabolite.annotation[anno]):
                m = re.search("^(cpd\d+)", metabolite.annotation[anno])
                return m[1]
        return None

    @staticmethod
    def reaction_msid(reaction):
        if re.search("^(rxn\d+)", reaction.id):
            m = re.search("^(rxn\d+)", reaction.id)
            return m[1]
        for anno in reaction.annotation:
            if isinstance(reaction.annotation[anno], list):
                for item in reaction.annotation[anno]:
                    if re.search("^(rxn\d+)", item):
                        m = re.search("^(rxn\d+)", item)
                        return m[1]
            elif re.search("^(rxn\d+)", reaction.annotation[anno]):
                m = re.search("^(rxn\d+)", reaction.annotation[anno])
                return m[1]
        return None

    @staticmethod
    def stoichiometry_to_string(stoichiometry):
        reactants = []
        products = []
        for met in stoichiometry:
            coef = stoichiometry[met]
            if not isinstance(met, str):
                if MSModelUtil.metabolite_msid(met) == "cpd00067":
                    met = None
                else:
                    met = met.id
            if met != None:
                if coef < 0:
                    reactants.append(met)
                else:
                    products.append(met)
        reactants.sort()
        products.sort()
        return [
            "+".join(reactants) + "=" + "+".join(products),
            "+".join(products) + "=" + "+".join(reactants),
        ]

    @staticmethod
    def search_name(name):
        name = name.lower()
        name = re.sub(r"_[a-z]\d*$", "", name)
        name = re.sub(r"\W+", "", name)
        return name

    @staticmethod
    def get(model, create_if_missing=True):
        if isinstance(model, MSModelUtil):
            return model
        if model in MSModelUtil.mdlutls:
            return MSModelUtil.mdlutls[model]
        elif create_if_missing:
            MSModelUtil.mdlutls[model] = MSModelUtil(model)
            return MSModelUtil.mdlutls[model]
        else:
            return None

    @staticmethod
    def build_from_kbase_json_file(filename, kbaseapi):
        """
        Builds an MSModelUtil object from a KBase JSON file.
        
        Args:
            filename (str): The path to the KBase JSON file.
            kbaseapi (KBaseAPI): An instance of the KBase API.

        Returns:
            An MSModelUtil object representing the contents of the KBase JSON file.
        """
        factory = kbaseapi.KBaseObjectFactory()
        model = factory.build_object_from_file(filename, "KBaseFBA.FBAModel")
        return MSModelUtil(model)

    def __init__(self, model):
        self.model = model
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.wsid = None
        self.atputl = None
        self.gfutl = None
        self.metabolite_hash = None
        self.search_metabolite_hash = None
        self.test_objective = None
        self.reaction_scores = None
        self.score = None
        self.breaking_reaction = None
        self.integrated_gapfillings = []
        self.attributes = {}
        self.atp_tests = None
        if hasattr(self.model, "computed_attributes"):
            if self.model.computed_attributes:
                self.attributes = self.model.computed_attributes
        if "pathways" not in self.attributes:
            self.attributes["pathways"] = {}
        if "auxotrophy" not in self.attributes:
            self.attributes["auxotrophy"] = {}
        if "fbas" not in self.attributes:
            self.attributes["fbas"] = {}

    ########Functions related to ATP gapfilling method
    def get_atputl(self,atp_media_filename=None,core_template=None,gapfilling_delta=0,max_gapfilling=0,forced_media=[],remake_atputil=False):
        """
        Returns and creates, if needed, an atp correction object for the model

        Parameters
        ----------
        core_template (optional) : MSTemplate object with core reactions
        atp_media_filename (optional) : string to tsv file with ATP media formulations
        gapfilling_delta (optional) : maximum difference in gapfilling to accept ATP condition
        max_gapfilling (optional) : maximum gapfilling allowable to accept an ATP growth condition
        forced_media (optional) : list of media in which model MUST make ATP
        
        Returns
        -------
        MSATPCorrection : Object for ATP correction

        Raises
        ------
        """
        if not self.atputl or remake_atputil:
            from modelseedpy.core.msatpcorrection import MSATPCorrection 
            self.atputl = MSATPCorrection(
                self,core_template,[],
                load_default_medias=True,
                max_gapfilling=max_gapfilling,
                gapfilling_delta=gapfilling_delta,
                forced_media=forced_media,
                default_media_path=atp_media_filename
            )
            self.atputl = MSATPCorrection(self.model)
        return self.atputl
    
    def get_atp_tests(self,core_template=None,atp_media_filename=None,recompute=False,remake_atputil=False):
        """
        Attempts to get ATP tests from attributes and failing that compute denovo using MSATPCorrection

        Parameters
        ----------
        core_template (optional) : MSTemplate object with core reactions
        atp_media_filename (optional) : string to tsv file with ATP media formulations

        Returns
        -------
        list<{"media":obj media,"is_max_threshold":bool,"threshold":float,"objective":string}>
            List of test specifications

        Raises
        ------
        """
        #Creating MSATPCorrection object which we need regardless
        atpcorrection = self.get_atputl(core_template=core_template,atp_media_filename=atp_media_filename,remake_atputil=remake_atputil)
        #Returning cached tests if available
        if self.atp_tests and not recompute:
            return self.atp_tests
        #Attempting to pull ATP tests from attributes
        if not recompute:
            print("Getting tests from attributes")
            atp_analysis = self.get_attributes("ATP_analysis",None)
            if atp_analysis:
                if "tests" in atp_analysis:
                    self.atp_tests = []
                    for item in atp_analysis["tests"]:
                        if item in atpcorrection.media_hash:
                            self.atp_tests.append({
                                "media":atpcorrection.media_hash[item],
                                "is_max_threshold":True,
                                "threshold":atp_analysis["tests"][item]["threshold"],
                                "objective":atp_analysis["tests"][item]["objective"]
                            })
                    return self.atp_tests
                else:
                    logger.warning("tests attribute missing in ATP analysis. Must recalculate ATP tests!")
            else:
                logger.warning("ATP analysis attributes missing. Must recalculate ATP tests!")
        #If recompute called for or if attributes are missing, recompute tests
        if not core_template:
            logger.warning("Cannot recompute ATP tests without a core template!")
            return None
        self.atp_tests = atpcorrection.build_tests()
        return self.atp_tests   

    def compute_automated_reaction_scores(self):
        """
        Computes reaction scores automatically from model data
        :return:
        """
        self.reaction_scores = {}

    def printlp(self, lpfilename="debug.lp"):
        with open(lpfilename, "w") as out:
            out.write(str(self.model.solver))

    def build_metabolite_hash(self):
        self.metabolite_hash = {}
        self.search_metabolite_hash = {}
        for met in self.model.metabolites:
            self.add_name_to_metabolite_hash(met.id, met)
            self.add_name_to_metabolite_hash(met.name, met)
            for anno in met.annotation:
                if isinstance(met.annotation[anno], list):
                    for item in met.annotation[anno]:
                        self.add_name_to_metabolite_hash(item, met)
                else:
                    self.add_name_to_metabolite_hash(met.annotation[anno], met)

    def add_name_to_metabolite_hash(self, name, met):
        if name not in self.metabolite_hash:
            self.metabolite_hash[name] = []
        self.metabolite_hash[name].append(met)
        sname = MSModelUtil.search_name(name)
        if sname not in self.search_metabolite_hash:
            self.search_metabolite_hash[sname] = []
        self.search_metabolite_hash[sname].append(met)

    def find_met(self, name, compartment=None):
        if self.metabolite_hash == None:
            self.build_metabolite_hash()
        if name in self.metabolite_hash:
            if not compartment:
                return self.metabolite_hash[name]
            for met in self.metabolite_hash[name]:
                array = met.id.split("_")
                if array[1] == compartment or met.compartment == compartment:
                    return [met]
            return None
        sname = MSModelUtil.search_name(name)
        if sname in self.search_metabolite_hash:
            if not compartment:
                return self.search_metabolite_hash[sname]
            for met in self.search_metabolite_hash[sname]:
                array = met.id.split("_")
                if array[1] == compartment or met.compartment == compartment:
                    return [met]
            return None
        logger.info(name + " not found in model!")
        return []

    def rxn_hash(self):
        output = {}
        for rxn in self.model.reactions:
            strings = MSModelUtil.stoichiometry_to_string(rxn.metabolites)
            output[strings[0]] = [rxn, 1]
            output[strings[1]] = [rxn, -1]
        return output

    def find_reaction(self, stoichiometry):
        output = MSModelUtil.stoichiometry_to_string(stoichiometry)
        atpstring = output[0]
        rxn_hash = self.rxn_hash()
        if atpstring in rxn_hash:
            return rxn_hash[atpstring]
        return None

    def msid_hash(self):
        output = {}
        for cpd in self.model.metabolites:
            msid = MSModelUtil.metabolite_msid(cpd)
            if msid != None:
                if msid not in output:
                    output[msid] = []
                output[msid].append(cpd)
        return output

    def exchange_list(self):
        exchange_reactions = []
        for reaction in self.model.reactions:
            if reaction.id[:3] == "EX_":
                exchange_reactions.append(reaction)
        return exchange_reactions

    def nonexchange_reaction_count(self):
        count = 0
        for reaction in self.model.reactions:
            if (
                reaction.id[:3] != "EX_"
                and reaction.id[:3] != "SK_"
                and reaction.id[:3] != "DM_"
                and reaction.id[:3] != "bio"
            ):
                if reaction.upper_bound > 0 or reaction.lower_bound < 0:
                    count += 1
        return count

    def exchange_hash(self):
        exchange_reactions = {}
        exlist = self.exchange_list()
        for reaction in exlist:
            for met in reaction.metabolites:
                if reaction.metabolites[met] == -1:
                    exchange_reactions[met] = reaction
                else:
                    logger.warn("Nonstandard exchange reaction ignored:" + reaction.id)
        return exchange_reactions

    def add_missing_exchanges(self, media):
        output = []
        exchange_hash = self.exchange_hash()
        exchange_list = []
        self.build_metabolite_hash()
        for mediacpd in media.mediacompounds:
            mets = self.find_met(mediacpd.id)
            if len(mets) > 0:
                found = 0
                cpd = None
                for met in mets:
                    if met in exchange_hash:
                        found = 1
                    elif met.compartment[0:1] == "c":
                        # We prefer to add a transport for the cytosol compound
                        cpd = met
                if cpd == None:
                    # No cytosol compound exists so choosing the first version we found that does exist
                    cpd = mets[0]
                if found == 0:
                    # No transporter currently exists - adding exchange reaction for the compound that does exist
                    output.append(cpd.id)
                    exchange_list.append(cpd)
        if len(exchange_list) > 0:
            self.add_exchanges_for_metabolites(exchange_list)
        return output

    def add_exchanges_for_metabolites(
        self, cpds, uptake=0, excretion=0, prefix="EX_", prefix_name="Exchange for "
    ):
        drains = []
        for cpd in cpds:
            drain_reaction = Reaction(
                id=f"{prefix}{cpd.id}",
                name=prefix_name + cpd.name,
                lower_bound=-1 * uptake,
                upper_bound=excretion,
            )
            drain_reaction.add_metabolites({cpd: -1})
            drain_reaction.annotation["sbo"] = "SBO:0000627"
            if drain_reaction.id not in self.model.reactions:
                drains.append(drain_reaction)
        self.model.add_reactions(drains)
        return drains

    def reaction_scores(self):
        return {}

    #################################################################################
    # Functions related to editing the model
    #################################################################################
    def get_attributes(self, key=None, default=None):
        if not key:
            return self.attributes
        if key not in self.attributes:
            self.attributes[key] = default
        return self.attributes[key]

    def save_attributes(self, value=None, key=None):
        if value:
            if key:
                self.attributes[key] = value
            else:
                self.attributes = value
        if hasattr(self.model, "computed_attributes"):
            logger.info("Setting FBAModel computed_attributes to mdlutl attributes")
            self.attributes["gene_count"] = len(self.model.genes)
            self.model.computed_attributes = self.attributes

    def add_ms_reaction(self, rxn_dict, compartment_trans=["c0", "e0"]):
        modelseed = ModelSEEDBiochem.get()
        output = []
        for rxnid, compartment in rxn_dict.items():
            fullid = rxnid + "_" + compartment
            modelseed_reaction = modelseed.get_seed_reaction(rxnid)
            reaction_stoich = modelseed_reaction.cstoichiometry
            cobra_reaction = Reaction(fullid)
            output.append(cobra_reaction)
            cobra_reaction.name = modelseed_reaction.data["name"] + "_" + compartment
            metabolites_to_add = {}
            for metabolite, stoich in reaction_stoich.items():
                id = metabolite[0]
                compound = modelseed.get_seed_compound(id).data
                compartment_number = int(metabolite[1])
                if compartment_number > len(compartment_trans):
                    logger.critical(
                        "Compartment index " + str(compartment_number) + " out of range"
                    )
                compartment_string = compartment_trans[compartment_number]
                met_output = self.find_met(id, compartment_string)
                cobramet = None
                if met_output:
                    cobramet = met_output[0]
                else:
                    cobramet = Metabolite(
                        id + "_" + compartment_string,
                        name=compound["name"] + "_" + compartment_string,
                        compartment=compartment_string,
                    )
                metabolites_to_add[cobramet] = stoich
            cobra_reaction.add_metabolites(metabolites_to_add)
            cobra_reaction.reaction
        self.model.add_reactions(output)
        return output

    #################################################################################
    # Functions related to utility functions
    #################################################################################
    def build_model_data_hash(self):
        data = {
            "Model": self.id,
            "Genome": self.genome.info.metadata["Name"],
            "Genes": self.genome.info.metadata["Number of Protein Encoding Genes"],
        }
        return data

    def compare_reactions(self, reaction_list, filename):
        data = {}
        for rxn in reaction_list:
            for met in rxn.metabolites:
                if met.id not in data:
                    data[met.id] = {}
                    for other_rxn in reaction_list:
                        data[met.id][other_rxn.id] = 0
                data[met.id][rxn.id] = rxn.metabolites[met]
        df = pd.DataFrame(data)
        df = df.transpose()
        df.to_csv(filename)

    #################################################################################
    # Functions related to managing biomass reactions
    #################################################################################
    def evaluate_biomass_reaction_mass(self, biomass_rxn_id, normalize=False):
        biorxn = self.model.reactions.get_by_id(biomass_rxn_id)
        # First computing energy biosynthesis coefficients
        atp = None
        atp_compounds = {
            "cpd00002": -1,
            "cpd00001": -1,
            "cpd00008": 1,
            "cpd00009": 1,
            "cpd00067": 1,
        }
        mass_compounds = {"cpd11463": 1, "cpd11461": 1, "cpd11462": 1}
        process_compounds = {"cpd17041": 1, "cpd17042": 1, "cpd17043": 1}
        for met in biorxn.metabolites:
            msid = self.metabolite_msid(met)
            if msid == "cpd00008":
                atp = abs(biorxn.metabolites[met])
        # Computing non ATP total mass
        total = 0
        for met in biorxn.metabolites:
            msid = self.metabolite_msid(met)
            if msid == "cpd11416":
                continue
            coef = biorxn.metabolites[met]
            if msid in mass_compounds:
                total += coef
            elif msid in process_compounds:
                total += 0
            else:
                mw = FBAHelper.metabolite_mw(met)
                if msid in atp_compounds:
                    if coef < 0:
                        coef += atp
                    else:
                        coef += -1 * atp
                total += mw * coef / 1000
        return {"ATP": atp, "Total": total}

    # Required this function to add gapfilled compounds to a KBase model for saving gapfilled model
    def convert_cobra_compound_to_kbcompound(self, cpd, kbmodel, add_to_model=1):
        refid = "cpd00000"
        if re.search("cpd\d+_[a-z]+", cpd.id):
            refid = cpd.id
            refid = re.sub("_[a-z]\d+$", "", refid)
        cpd_data = {
            "aliases": [],
            "charge": cpd.charge,
            "compound_ref": "~/template/compounds/id/" + refid,
            "dblinks": {},
            "formula": cpd.formula,
            "id": cpd.id,
            "modelcompartment_ref": "~/modelcompartments/id/" + cpd.id.split("_").pop(),
            "name": cpd.name,
            "numerical_attributes": {},
            "string_attributes": {},
        }
        if add_to_model == 1:
            kbmodel["modelcompounds"].append(cpd_data)
        return cpd_data

    def compute_flux_values_from_variables(self):
        """Returns a hash of reaction fluxes from model object

        Parameters
        ----------

        Returns
        -------
        dict<string reaction ID,{'reverse':float flux,'forward':float flux}>
            Hash of reactions and their associated flux values

        Raises
        ------
        """
        flux_values = {}
        for rxn in self.model.reactions:
            flux_values[rxn.id] = {
                "reverse": rxn.reverse_variable.primal,
                "forward": rxn.forward_variable.primal,
            }
        return flux_values

    # Required this function to add gapfilled reactions to a KBase model for saving gapfilled model
    def convert_cobra_reaction_to_kbreaction(
        self, rxn, kbmodel, cpd_hash, direction="=", add_to_model=1, reaction_genes=None
    ):
        rxnref = "~/template/reactions/id/rxn00000_c"
        if re.search("rxn\d+_[a-z]+", rxn.id):
            rxnref = "~/template/reactions/id/" + rxn.id
            rxnref = re.sub("\d+$", "", rxnref)
        rxn_data = {
            "id": rxn.id,
            "aliases": [],
            "dblinks": {},
            "direction": direction,
            "edits": {},
            "gapfill_data": {},
            "maxforflux": 1000000,
            "maxrevflux": 1000000,
            "modelReactionProteins": [],
            "modelReactionReagents": [],
            "modelcompartment_ref": "~/modelcompartments/id/" + rxn.id.split("_").pop(),
            "name": rxn.name,
            "numerical_attributes": {},
            "probability": 0,
            "protons": 0,
            "reaction_ref": rxnref,
            "string_attributes": {},
        }
        for cpd in rxn.metabolites:
            if cpd.id not in kbmodel["modelcompounds"]:
                cpd_hash[cpd.id] = self.convert_cobra_compound_to_kbcompound(
                    cpd, kbmodel, 1
                )
            rxn_data["modelReactionReagents"].append(
                {
                    "coefficient": rxn.metabolites[cpd],
                    "modelcompound_ref": "~/modelcompounds/id/" + cpd.id,
                }
            )
        if reaction_genes != None and rxn.id in reaction_genes:
            best_gene = None
            for gene in reaction_genes[rxn.id]:
                if (
                    best_gene == None
                    or reaction_genes[rxn.id][gene] > reaction_genes[rxn.id][best_gene]
                ):
                    best_gene = gene
            rxn_data["modelReactionProteins"] = [
                {
                    "note": "Added from gapfilling",
                    "modelReactionProteinSubunits": [],
                    "source": "Unknown",
                }
            ]
            rxn_data["modelReactionProteins"][0]["modelReactionProteinSubunits"] = [
                {
                    "note": "Added from gapfilling",
                    "optionalSubunit": 0,
                    "triggering": 1,
                    "feature_refs": ["~/genome/features/id/" + best_gene],
                    "role": "Unknown",
                }
            ]
        if add_to_model == 1:
            kbmodel["modelreactions"].append(rxn_data)
        return rxn_data

    #################################################################################
    # Functions related to gapfilling of models
    #################################################################################
    def convert_solution_to_list(self,solution):
        """Converting solution to list format, which is easier to work with
        Parameters
        ----------
        solution : dict
            Specifies the reactions to be added to the model to implement the gapfilling solution
        """
        output = []
        for label in ["new","reversed"]:
            for rxn_id in solution[label]:
                output.append([rxn_id, solution[label][rxn_id],label])
        return output

    def find_item_in_solution(self,input_list,input,ignore_dir=False):
        for item in input_list:
            if input[0] == item[0] and input[1] == item[1]:
                return True
            elif ignore_dir and input[0] == item[0]:
                return True
        return False
    
    def test_solution(self,solution,targets,medias,thresholds=[0.1],remove_unneeded_reactions=False,do_not_remove_list=[]):
        """Tests if every reaction in a given gapfilling solution is actually needed for growth. Note, this code assumes the gapfilling solution is already integrated.

        Parameters
        ----------
        solution : {"new":{string reaction_id: string direction},"reversed":{string reaction_id: string direction}}
                    or 
                    list<tuple<string - reaction id, string direction>>
            Data for gapfilling solution to be tested
        target : string,
        media : MSMedia,
        threshold : float, default 0.1

        Returns
        -------
        list<tuple<string - reaction id, string direction>>
            List of unneeded reactions

        Raises
        ------
        """
        #Saving the current objective
        current_objective = self.model.objective
        #Saving the current media
        current_media = self.pkgmgr.getpkg("KBaseMediaPkg").current_media
        #Computing the initial objective values
        initial_objectives = []
        for (i,target) in enumerate(targets):
            #Setting the media
            self.pkgmgr.getpkg("KBaseMediaPkg").build_package(medias[i])
            #Setting the objective
            self.model.objective = target
            #Computing the objective value
            objective = self.model.slim_optimize()
            initial_objectives.append(objective)
            logger.debug("Starting objective for " + medias[i].id + "/"+target+" = " + str(objective))
        #Iterating through solution reactions and flagging them if they are unneeded to achieve the specified minimum objective
        unneeded = []
        #If object is a dictionary, convert to a list
        if isinstance(solution,dict):
            solution = self.convert_solution_to_list(solution)
        #Processing solution in standardized format
        for item in solution:
            rxn_id = item[0]    
            rxnobj = self.model.reactions.get_by_id(rxn_id)
            #Testing all media and target and threshold combinations to see if the reaction is needed
            needed = False
            for (i,target) in enumerate(targets):
                if len(targets) > 1:#If there's only one target, then these steps were done earlier already
                    #Setting the media
                    self.pkgmgr.getpkg("KBaseMediaPkg").build_package(medias[i])
                    #Setting the objective
                    self.model.objective = target
                #Knocking out the reaction to test for the impact on the objective
                #This has to happen after media is applied in case the reaction is an exchange
                if item[1] == ">":
                    original_bound = rxnobj.upper_bound
                    rxnobj.upper_bound = 0
                else:
                    original_bound = rxnobj.lower_bound
                    rxnobj.lower_bound = 0
                #Computing the objective value
                objective = self.model.slim_optimize()
                if objective < thresholds[i]:
                    needed = True
                    logger.info(
                        medias[i].id + "/" + target + ":" +rxn_id
                        + item[1]
                        + " needed:"
                        + str(objective)
                        + " with min obj:"
                        + str(thresholds[i])
                    )
            #If the reaction isn't needed for any media and target combinations, add it to the unneeded list
            if not needed:
                unneeded.append([rxn_id, item[1], item[2],original_bound])
                logger.info(
                    rxn_id
                    + item[1]
                    + " not needed:"
                    + str(objective)
                )
                #VERY IMPORTANT: Leave the reaction knocked out for now so we screen for combinatorial effects
            else:
                #Restore the reaction if it is needed
                if item[1] == ">":
                    rxnobj.upper_bound = original_bound
                else:
                    rxnobj.lower_bound = original_bound
        if not remove_unneeded_reactions:
            #Restoring the bounds on the unneeded reactions
            for item in unneeded:
                rxnobj = self.model.reactions.get_by_id(item[0])
                if item[1] == ">":
                    rxnobj.upper_bound = item[3]
                else:
                    rxnobj.lower_bound = item[3]
        else:
            #Do not restore bounds on unneeded reactions and remove reactions from model if their bounds are zero
            removed_rxns = []
            for item in unneeded:
                rxnobj = self.model.reactions.get_by_id(item[0])
                if self.find_item_in_solution(do_not_remove_list,item):
                    if item[1] == ">":
                        rxnobj.upper_bound = item[3]
                    else:
                        rxnobj.lower_bound = item[3]
                elif rxnobj.lower_bound == 0 and rxnobj.upper_bound == 0 and not self.find_item_in_solution(do_not_remove_list,item,ignore_dir=True):
                    removed_rxns.append(rxnobj)
            if len(removed_rxns) > 0:
                self.model.remove_reactions(removed_rxns)
        #Restoring the original objective
        self.model.objective = current_objective
        #Restoring the original media
        if current_media:
            self.pkgmgr.getpkg("KBaseMediaPkg").build_package(current_media)
        #Returning the unneeded list
        return unneeded

    def add_gapfilling(self, solution):
        print("Adding gapfilling",str(solution))
        self.integrated_gapfillings.append(solution)

    def create_kb_gapfilling_data(self, kbmodel, atpmedia_ws="94026"):
        gapfilling_hash = {}
        if "gapfillings" not in kbmodel:
            kbmodel["gapfillings"] = []
        for gapfilling in kbmodel["gapfillings"]:
            gapfilling_hash[gapfilling["id"]] = gapfilling
        rxn_hash = {}
        for rxn in kbmodel["modelreactions"]:
            rxn_hash[rxn["id"]] = rxn
        for gf in self.integrated_gapfillings:
            media_ref = "KBaseMedia/Empty"
            gf["media"].id.replace("/", ".")
            gfid = gf["media"].id
            if self.atputl:
                for item in self.atputl.atp_medias:
                    if item[0] == gf["media"]:
                        gfid = "ATP-" + gfid
                        media_ref = atpmedia_ws + "/" + gf["media"].id + ".atp"
                        break
            if hasattr(gf["media"], "info"):
                media_ref = gf["media"].info.workspace_id + "/" + gf["media"].info.id
            suffix = 0
            while gfid in gapfilling_hash:
                suffix += 1
                gfid += "." + str(suffix)
            gapfilling_hash[gfid] = 1
            gapfilling_obj = {
                "gapfill_id": gfid,
                "id": gfid,
                "integrated": 1,
                "integrated_solution": "0",
                "target": gf["target"],
                "minobjective": gf["minobjective"],
                "binary_check": gf["binary_check"],
                "media_ref": media_ref,
            }
            kbmodel["gapfillings"].append(gapfilling_obj)
            for rxn in gf["new"]:
                if rxn in rxn_hash:
                    rxnobj = rxn_hash[rxn]
                    if "gapfill_data" not in rxnobj:
                        rxnobj["gapfill_data"] = {}
                    if gfid not in rxnobj["gapfill_data"]:
                        rxnobj["gapfill_data"][gfid] = {"0": [gf["new"][rxn], 1, []]}
            for rxn in gf["reversed"]:
                if rxn in rxn_hash:
                    rxnobj = rxn_hash[rxn]
                    if "gapfill_data" not in rxnobj:
                        rxnobj["gapfill_data"] = {}
                    if gfid not in rxnobj["gapfill_data"]:
                        rxnobj["gapfill_data"][gfid] = {
                            "0": [gf["reversed"][rxn], 1, []]
                        }

    #################################################################################
    # Functions related to applying, running, and expanding with test conditions
    #################################################################################
    def apply_test_condition(self, condition, model=None):
        """Applies constraints and objective of specified condition to model

        Parameters
        ----------
        condition : dict
            Specifies condition to be tested with media, objective, is_max_threshold, threshold.
        model : cobra.Model, optional
            Specific instance of model to apply conditions to (useful if using "with model")

        Returns
        -------
        boolean
            True if threshold is NOT exceeded, False if threshold is exceeded

        Raises
        ------
        """
        if model is None:
            model = self.model
            pkgmgr = self.pkgmgr
        else:
            pkgmgr = MSPackageManager.get_pkg_mgr(model)
        model.objective = condition["objective"]
        #if condition["is_max_threshold"]:
        model.objective.direction = "max"
        #else: TODO - need to revisit this
        #    model.objective.direction = "min"
        pkgmgr.getpkg("KBaseMediaPkg").build_package(condition["media"])

    def test_single_condition(self, condition, apply_condition=True, model=None):
        """Runs a single test condition to determine if objective value on set media exceeds threshold

        Parameters
        ----------
        condition : dict
            Specifies condition to be tested with media, objective, is_max_threshold, threshold.
        apply_condition : bool,optional
            Indicates if condition constraints and objective should be applied.
        model : cobra.Model, optional
            Specific instance of model to apply tests to (useful if using "with model")

        Returns
        -------
        boolean
            True if threshold is NOT exceeded, False if threshold is exceeded

        Raises
        ------
        """
        if model is None:
            model = self.model
        if apply_condition:
            self.apply_test_condition(condition, model)
        new_objective = model.slim_optimize()
        value = new_objective
        if "change" in condition and condition["change"]:
            if self.test_objective:
                value = new_objective - self.test_objective
                logger.debug(
                    condition["media"].id
                    + " testing for change:"
                    + str(value)
                    + "="
                    + str(new_objective)
                    + "-"
                    + str(self.test_objective)
                )
        self.score = value
        if model.solver.status != "optimal":
            self.printlp(condition["media"].id + "-Testing-Infeasible.lp")
            logger.critical(
                condition["media"].id
                + "testing leads to infeasible problem. LP file printed to debug!"
            )
            return False
        if value >= condition["threshold"] and condition["is_max_threshold"]:
            logger.debug("Failed high:"+condition["media"].id+":"+str(new_objective)+";"+str(condition["threshold"]))
            return False
        elif value <= condition["threshold"] and not condition["is_max_threshold"]:
            logger.info("Failed low:"+condition["media"].id+":"+str(new_objective)+";"+str(condition["threshold"]))
            return False
        self.test_objective = new_objective
        logger.debug("Passed:"+condition["media"].id+":"+str(new_objective)+";"+str(condition["threshold"]))
        return True

    def test_condition_list(self, condition_list, model=None,positive_growth=[]):
        """Runs a set of test conditions to determine if objective values on set medias exceed thresholds

        Parameters
        ----------
        condition_list : list<dict>
            Specifies set of conditions to be tested with media, objective, is_max_threshold, threshold.
        model : cobra.Model, optional
            Specific instance of model to apply tests to (useful if using "with model")

        Returns
        -------
        boolean
            True if ALL tests pass, False if any test returns false

        Raises
        ------
        """
        if model == None:
            model = self.model
        for condition in condition_list:
            if not self.test_single_condition(condition,apply_condition=True,model=model):
                return False
        return True

    def linear_expansion_test(self, reaction_list, condition, currmodel,positive_growth=[]):
        """Tests addition of reactions one at a time

        Parameters
        ----------
        reaction_list : list<[obj reaction,{>|>}]>
            List of reactions and directions to test for addition in the model (should already be in model)

        Returns
        -------
        list<[obj reaction,{>|>}]>
            List of reactions and directions filtered because they fail tests when in the model

        Raises
        ------
        """
        # First run the full test
        if self.test_single_condition(condition, apply_condition=False, model=currmodel,positive_growth=positive_growth):
            return []
        # First knockout all reactions in the input list and save original bounds
        filtered_list = []
        original_bound = []
        for item in reaction_list:
            if item[1] == ">":
                original_bound.append(item[0].upper_bound)
                item[0].upper_bound = 0
            else:
                original_bound.append(item[0].lower_bound)
                item[0].lower_bound = 0
        # Now restore reactions one at a time
        count = 0
        for item in reaction_list:
            if item[1] == ">":
                item[0].upper_bound = original_bound[count]
                if not self.test_single_condition(condition, apply_condition=False, model=currmodel):
                    # logger.debug(item[0].id+":"+item[1])
                    item[0].upper_bound = 0
                    if item not in filtered_list:
                        item.append(original_bound[count])
                        item.append(self.score)
                        filtered_list.append(item)
            else:
                item[0].lower_bound = original_bound[count]
                if not self.test_single_condition(condition, apply_condition=False, model=currmodel):
                    # logger.debug(item[0].id+":"+item[1])
                    item[0].lower_bound = 0
                    if item not in filtered_list:
                        item.append(original_bound[count])
                        item.append(self.score)
                        filtered_list.append(item)
            count += 1
        return filtered_list

    def binary_expansion_test(self, reaction_list, condition, currmodel, depth=0,positive_growth=[]):
        """Conducts a binary search for bad reaction combinations
        Parameters
        ----------
        reaction_list : list<[obj reaction,{>|>}]>
            List of reactions and directions to test for addition in the model (should already be in model)
        condition_list : list<dict>
            Specifies set of conditions to be tested with media, objective, is_max_threshold, threshold.

        Returns
        -------
        list<[obj reaction,{>|>}]>
            List of reactions and directions filtered because they fail tests when in the model

        Raises
        ------
        """
        newdepth = depth + 1
        filtered_list = []
        # First run the full test
        if self.test_single_condition(condition,apply_condition=False,model=currmodel):
            return []
        # Check if input list contains only one reaction:
        if len(reaction_list) == 1:
            print("Failed:"+reaction_list[0][1]+reaction_list[0][0].id)
            if reaction_list[0][1] == ">":
                reaction_list[0].append(reaction_list[0][0].upper_bound)
                reaction_list[0][0].upper_bound = 0
            else:
                reaction_list[0].append(reaction_list[0][0].lower_bound)
                reaction_list[0][0].lower_bound = 0
            #Check if the reaction passes the positive growth test
            success = True
            if len(positive_growth) > 0:
                #Testing positive growth conditions
                for pos_condition in positive_growth:
                    if not self.test_single_condition(pos_condition,apply_condition=True,model=currmodel):
                        print("Does not pass positive growth tests:"+reaction_list[0][1]+reaction_list[0][0].id)
                        success = False
                        break
                #Restoring current test condition
                self.apply_test_condition(condition)
            if success:
                reaction_list[0].append(self.score)
                filtered_list.append(reaction_list[0])
            else:
                #Restoring reaction
                if reaction_list[0][1] == ">":
                    reaction_list[0][0].upper_bound = reaction_list[0][2]
                else:
                    reaction_list[0][0].lower_bound = reaction_list[0][2]
                self.breaking_reaction = reaction_list[0][0]
            return filtered_list
        # Break reaction list into two
        original_bound = []
        sub_lists = [[], []]
        midway_point = int(len(reaction_list) / 2)
        for i, item in enumerate(reaction_list):
            if item[1] == ">":
                original_bound.append(item[0].upper_bound)
            else:
                original_bound.append(item[0].lower_bound)
            if i < midway_point:
                sub_lists[0].append(item)
            else:
                sub_lists[1].append(item)
                if item[1] == ">":
                    item[0].upper_bound = 0
                else:
                    item[0].lower_bound = 0
        # Submitting first half of reactions for testing
        new_filter = self.binary_expansion_test(
            sub_lists[0], condition, currmodel,depth=newdepth,positive_growth=positive_growth
        )
        for item in new_filter:
            filtered_list.append(item)
        if self.breaking_reaction != None:
            print("Ending early due to breaking reaction:"+self.breaking_reaction.id)
            return filtered_list
        # Submitting second half of reactions for testing - now only breaking reactions are removed from the first list
        for i, item in enumerate(reaction_list):
            if i >= midway_point:
                if item[1] == ">":
                    item[0].upper_bound = original_bound[i]
                else:
                    item[0].lower_bound = original_bound[i]
        new_filter = self.binary_expansion_test(
            sub_lists[1], condition, currmodel,depth=newdepth,positive_growth=positive_growth
        )
        for item in new_filter:
            filtered_list.append(item)
        return filtered_list

    def check_if_solution_exists(self, reaction_list, condition, model):
        original_bound = []
        for i, item in enumerate(reaction_list):
            if item[1] == ">":
                original_bound.append(item[0].upper_bound)
                item[0].upper_bound = 0
            else:
                original_bound.append(item[0].lower_bound)
                item[0].lower_bound = 0
        result = self.test_single_condition(condition,model=model)
        for i, item in enumerate(reaction_list):
            if item[1] == ">":
                item[0].upper_bound = original_bound[i]
            else:
                item[0].lower_bound = original_bound[i]
        return result          
    
    def reaction_expansion_test(
        self,
        reaction_list,
        condition_list,
        binary_search=True,
        attribute_label="gf_filter",
        positive_growth=[]
    ):
        """Adds reactions in reaction list one by one and appplies tests, filtering reactions that fail

        Parameters
        ----------
        reaction_list : list<[obj reaction,{>|>}]>
            List of reactions and directions to test for addition in the model (should already be in model)
        condition_list : list<dict>
            Specifies set of conditions to be tested with media, objective, is_max_threshold, threshold.

        Returns
        -------
        list<[obj reaction,{>|>}]>
            List of reactions and directions filtered because they fail tests when in the model

        Raises
        ------
        """
        logger.debug(f"Expansion started! Binary = {binary_search}")
        self.breaking_reaction = None
        filtered_list = []
        for condition in condition_list:
            logger.debug(f"testing condition {condition}")
            currmodel = self.model
            tic = time.perf_counter()
            new_filtered = []
            if not self.check_if_solution_exists(reaction_list, condition, currmodel):
                print("No solution exists that passes tests for condition "+condition["media"].id)
                return None
            with currmodel:
                self.apply_test_condition(condition)
                if binary_search:
                    done = False
                    while not done:
                        new_filtered = self.binary_expansion_test(
                            reaction_list, condition, currmodel,positive_growth=positive_growth
                        )
                        for item in new_filtered:
                            if item not in filtered_list:
                                filtered_list.append(item)
                        if self.breaking_reaction == None:
                            done = True
                        else:
                            #Remove breaking reaction from reaction_list
                            print("Keeping breaking reaction:"+self.breaking_reaction.id)
                            for i in range(len(reaction_list)):
                                if reaction_list[i][0] == self.breaking_reaction:
                                    del reaction_list[i]
                                    break
                            self.breaking_reaction = None
                            if not self.check_if_solution_exists(reaction_list, condition, currmodel):
                                print("No solution exists after retaining breaking reaction:"+self.breaking_reaction.id)
                                return None
                else:
                    new_filtered = self.linear_expansion_test(
                        reaction_list, condition, currmodel,positive_growth=positive_growth
                    )
                    for item in new_filtered:
                        if item not in filtered_list:
                            filtered_list.append(item)
            # Restoring knockout of newly filtered reactions, which expire after exiting the "with" block above
            for item in new_filtered:
                if item[1] == ">":
                    item[0].upper_bound = 0
                else:
                    item[0].lower_bound = 0
            toc = time.perf_counter()
            logger.info(
                "Expansion time:" + condition["media"].id + ":" + str((toc - tic))
            )
            logger.info(
                "Filtered count:"
                + str(len(filtered_list))
                + " out of "
                + str(len(reaction_list))
            )
            # Adding filter results to attributes
            gf_filter_att = self.get_attributes(attribute_label, {})
            if condition["media"].id not in gf_filter_att:
                gf_filter_att[condition["media"].id] = {}
            if condition["objective"] not in gf_filter_att[condition["media"].id]:
                gf_filter_att[condition["media"].id][condition["objective"]] = {}
            if (
                condition["threshold"]
                not in gf_filter_att[condition["media"].id][condition["objective"]]
            ):
                gf_filter_att[condition["media"].id][condition["objective"]][
                    condition["threshold"]
                ] = {}
            for item in new_filtered:
                if (
                    item[0].id
                    not in gf_filter_att[condition["media"].id][condition["objective"]][
                        condition["threshold"]
                    ]
                ):
                    gf_filter_att[condition["media"].id][condition["objective"]][
                        condition["threshold"]
                    ][item[0].id] = {}
                if (
                    item[1]
                    not in gf_filter_att[condition["media"].id][condition["objective"]][
                        condition["threshold"]
                    ][item[0].id]
                ):
                    if len(item) < 3:
                        gf_filter_att[condition["media"].id][condition["objective"]][
                            condition["threshold"]
                        ][item[0].id][item[1]] = None
                    else:
                        gf_filter_att[condition["media"].id][condition["objective"]][
                            condition["threshold"]
                        ][item[0].id][item[1]] = item[2]
        return filtered_list

    #################################################################################
    # Functions related to biomass sensitivity analysis
    #################################################################################
    def find_unproducible_biomass_compounds(self, target_rxn="bio1", ko_list=None):
        # Cloning the model because we don't want to modify the original model with this analysis
        tempmodel = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
        # Getting target reaction and making sure it exists
        if target_rxn not in tempmodel.reactions:
            logger.critical(target_rxn + " not in model!")
            return None
        target_rxn_obj = tempmodel.reactions.get_by_id(target_rxn)
        tempmodel.objective = target_rxn
        original_objective = tempmodel.objective
        pkgmgr = MSPackageManager.get_pkg_mgr(tempmodel)
        rxn_list = [target_rxn, "rxn05294_c0", "rxn05295_c0", "rxn05296_c0"]
        for rxn in rxn_list:
            if rxn in tempmodel.reactions:
                pkgmgr.getpkg("FlexibleBiomassPkg").build_package(
                    {
                        "bio_rxn_id": rxn,
                        "flex_coefficient": [0, 1],
                        "use_rna_class": None,
                        "use_dna_class": None,
                        "use_protein_class": None,
                        "use_energy_class": [0, 1],
                        "add_total_biomass_constraint": False,
                    }
                )

        # Creating min flex objective
        min_flex_obj = tempmodel.problem.Objective(Zero, direction="min")
        obj_coef = dict()
        for reaction in tempmodel.reactions:
            if reaction.id[0:5] == "FLEX_" or reaction.id[0:6] == "energy":
                obj_coef[reaction.forward_variable] = 1
                obj_coef[reaction.reverse_variable] = 1
        # Temporarily setting flex objective so I can set coefficients
        tempmodel.objective = min_flex_obj
        min_flex_obj.set_linear_coefficients(obj_coef)
        if not ko_list:
            return self.run_biomass_dependency_test(
                target_rxn_obj, tempmodel, original_objective, min_flex_obj, rxn_list
            )
        else:
            output = {}
            for item in ko_list:
                logger.debug("KO:" + item[0] + item[1])
                if item[0] not in output:
                    output[item[0]] = {}
                if item[0] in tempmodel.reactions:
                    rxnobj = tempmodel.reactions.get_by_id(item[0])
                    if item[1] == ">":
                        original_bound = rxnobj.upper_bound
                        rxnobj.upper_bound = 0
                        output[item[0]][item[1]] = self.run_biomass_dependency_test(
                            target_rxn_obj,
                            tempmodel,
                            original_objective,
                            min_flex_obj,
                            rxn_list,
                        )
                        rxnobj.upper_bound = original_bound
                    else:
                        original_bound = rxnobj.lower_bound
                        rxnobj.lower_bound = 0
                        output[item[0]][item[1]] = self.run_biomass_dependency_test(
                            target_rxn_obj,
                            tempmodel,
                            original_objective,
                            min_flex_obj,
                            rxn_list,
                        )
                        rxnobj.lower_bound = original_bound
                else:
                    logger.info("Reaction "+item[0]+" not in model during sensitivity analysis!")
                    output[item[0]][item[1]] = []
            return output

    def run_biomass_dependency_test(
        self, target_rxn, tempmodel, original_objective, min_flex_obj, rxn_list
    ):
        tempmodel.objective = original_objective
        objective = tempmodel.slim_optimize()
        if objective > 0:
            target_rxn.lower_bound = 0.1
            tempmodel.objective = min_flex_obj
            solution = tempmodel.optimize()
            biocpds = []
            for reaction in tempmodel.reactions:
                if reaction.id[0:5] == "FLEX_" and (
                    reaction.forward_variable.primal > Zero
                    or reaction.reverse_variable.primal > Zero
                ):
                    logger.debug("Depends on:" + reaction.id)
                    label = reaction.id[5:]
                    for item in rxn_list:
                        if label[0 : len(item)] == item:
                            biocpds.append(label[len(item) + 1 :])
            target_rxn.lower_bound = 0
            return biocpds
        else:
            logger.debug("Cannot grow")
            return None

    def add_atp_hydrolysis(self, compartment):
        # Searching for ATP hydrolysis compounds
        coefs = {
            "cpd00002": [-1, compartment],
            "cpd00001": [-1, compartment],
            "cpd00008": [1, compartment],
            "cpd00009": [1, compartment],
            "cpd00067": [1, compartment],
        }
        msids = ["cpd00002", "cpd00001", "cpd00008", "cpd00009", "cpd00067"]
        stoichiometry = {}
        id_hash = self.msid_hash()
        for msid in msids:
            if msid not in id_hash:
                logger.warning("Compound " + msid + " not found in model!")
                return None
            else:
                for cpd in id_hash[msid]:
                    if cpd.compartment == coefs[msid][1]:
                        stoichiometry[cpd] = coefs[msid][0]
        output = self.find_reaction(stoichiometry)
        if output != None and output[1] == ">":
            return {"reaction": output[0], "direction": ">", "new": False}
        cobra_reaction = Reaction(
            "rxn00062_" + compartment,
            name="ATP hydrolysis",
            lower_bound=0,
            upper_bound=1000,
        )
        cobra_reaction.annotation["sbo"] = "SBO:0000176"  # biochemical reaction
        cobra_reaction.annotation["seed.reaction"] = "rxn00062"
        cobra_reaction.add_metabolites(stoichiometry)
        self.model.add_reactions([cobra_reaction])
        return {"reaction": cobra_reaction, "direction": ">", "new": True}

    @staticmethod
    def parse_id(object):
        if re.search("(.+)_([a-z]+)(\d*)$", object.id) != None:
            m = re.search("(.+)_([a-z]+)(\d*)$", object.id)
            return (m[1], m[2], m[3])
        return None
