# -*- coding: utf-8 -*-
import logging
import re
import time
import json
import sys
from cobra import Model, Reaction, Metabolite
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem
from modelseedpy.core.exceptions import *
from modelseedpy.core.fbahelper import FBAHelper
from optlang.symbolics import Zero
from optlang import Constraint

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# handler = logging.StreamHandler(sys.stdout)
# handler.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)
# logger.addHandler(handler)


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
        if model in MSModelUtil.mdlutls:
            return MSModelUtil.mdlutls[model]
        elif create_if_missing:
            MSModelUtil.mdlutls[model] = MSModelUtil(model)
            return MSModelUtil.mdlutls[model]
        else:
            return None

    def __init__(self, model):
        self.model = model
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.atputl = None
        self.gfutl = None
        self.metabolite_hash = None
        self.search_metabolite_hash = None
        self.test_objective = None
        self.reaction_scores = None
        self.score = None
        self.integrated_gapfillings = []

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
        return [rxn for rxn in self.model.reactions if 'EX_' in rxn.id]

    def carbon_mets(self):
        return [met for met in self.model.metabolites if 'C' in met.elements]

    def carbon_exchange_list(self, include_unknown=True):
        if not include_unknown:
            return [ex for ex in self.exchange_list() if "C" in ex.reactants[0].elements]
        return [ex for ex in self.exchange_list() if not ex.reactants[0].elements or "C" in ex.reactants[0].elements]

    def exchange_mets_list(self):
        return [met for ex_rxn in self.exchange_list() for met in ex_rxn.metabolites]

    def media_exchanges_list(self):
        return [exRXN for exRXN in self.exchange_list() if exRXN.id in self.model.medium]

    def bio_rxns_list(self):
        return [rxn for rxn in self.model.reactions if re.search(r"(^bio\d+)", rxn.id)]

    def transport_list(self):
        return [rxn for rxn in self.model.reactions if any(["_e0" in met.id for met in rxn.metabolites])]

    def compatibilize(self, conflicts_file_name="orig_conflicts.json", printing=False):
        from modelseedpy.community.mscompatibility import MSCompatibility
        self.model = MSCompatibility.standardize(
            [self.model], conflicts_file_name=conflicts_file_name, printing=printing)[0]
        return self.model

    def standard_exchanges(self):
        for ex in self.exchange_list():
            if len(ex.reactants) != 1 and len(ex.products) != 0:
                raise ModelError(f"The ex {ex.id} possesses {len(ex.reactants)} reactants and "
                                 f"{len(ex.products)} products, which are non-standard and are incompatible"
                                 f" with various ModelSEED operations.")

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
        for ex_rxn in self.exchange_list():
            for met in ex_rxn.metabolites:
                if ex_rxn.metabolites[met] == -1:
                    exchange_reactions[met] = ex_rxn
                else:
                    logger.warn("Nonstandard exchange reaction ignored:" + reaction.id)
        return exchange_reactions
    
    def var_names_list(self):
        return [var.name for var in self.model.variables]

    def add_missing_exchanges(self, media):
        output = []
        exchange_hash = self.exchange_hash()
        exchange_list = []
        self.build_metabolite_hash()
        for mediacpd in media.mediacompounds:
            mets = self.find_met(mediacpd.id)
            if len(mets) > 0:
                found = False
                cpd = None
                for met in mets:
                    if met in self.exchange_hash():
                        found = True
                    elif met.compartment[0:1] == "c":
                        # We prefer to add a transport for the cytosol compound
                        cpd = met
                if cpd == None:
                    # No cytosol compound exists so choosing the first version we found that does exist
                    cpd = mets[0]
                if found:
                    #No transporter currently exists - adding exchange reaction for the compound that does exist
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
        
    def reaction_scores(self):  #!!! Can this be deleted?
        return {}
    
    #adding gapfilling compounds to a KBase model saves gapfilled models 
    def convert_cobra_compound_to_kbcompound(self, cpd, kbmodel, add_to_model=1):
        refid = "cpd00000"
        if re.search('cpd\d+_[a-z]+',cpd.id):
            refid = re.sub("_[a-z]\d+$","",cpd.id)

    #################################################################################
    # Functions related to editing the model
    #################################################################################
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
            print(cobra_reaction.id)
        print(len(output))
        self.model.add_reactions(output)
        return output

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
            rxnref = re.sub("\d+$","",f"~/template/reactions/id/{rxn.id}")
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
    """Tests if every reaction in a given gapfilling solution is actually needed for growth
        Optionally can remove unneeded reactions from the model AND the solution object.
        Note, this code assumes the gapfilling solution is already integrated.

        Parameters
        ----------
        {"new":{string reaction_id: string direction},"reversed":{string reaction_id: string direction}} - solution
            Data for gapfilling solution to be tested
        bool - keep_changes
            Set this bool to True to remove the unneeded reactions from the solution and model
        Returns
        -------
        list<tuple<string - reaction id, string direction>>
            List of unneeded reactions

        Raises
        ------
        """

    def test_solution(self, solution, keep_changes=False):
        unneeded = []
        tempmodel = self.model
        if not keep_changes:
            tempmodel = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
        tempmodel.objective = solution["target"]
        pkgmgr = MSPackageManager.get_pkg_mgr(tempmodel)
        pkgmgr.getpkg("KBaseMediaPkg").build_package(solution["media"])
        objective = tempmodel.slim_optimize()
        logger.debug("Starting objective:" + str(objective))
        types = ["new", "reversed"]
        for key in types:
            for rxn_id in solution[key]:
                rxnobj = tempmodel.reactions.get_by_id(rxn_id)
                if solution[key][rxn_id] == ">":
                    original_bound = rxnobj.upper_bound
                    rxnobj.upper_bound = 0
                    objective = tempmodel.slim_optimize()
                    if objective < solution["minobjective"]:
                        logger.debug(
                            rxn_id
                            + solution[key][rxn_id]
                            + " needed:"
                            + str(objective)
                            + " with min obj:"
                            + str(solution["minobjective"])
                        )
                        rxnobj.upper_bound = original_bound
                    else:
                        unneeded.append([rxn_id, solution[key][rxn_id], key])
                        logger.debug(
                            rxn_id
                            + solution[key][rxn_id]
                            + " not needed:"
                            + str(objective)
                        )
                else:
                    original_bound = rxnobj.lower_bound
                    rxnobj.lower_bound = 0
                    objective = tempmodel.slim_optimize()
                    if objective < solution["minobjective"]:
                        logger.debug(
                            rxn_id
                            + solution[key][rxn_id]
                            + " needed:"
                            + str(objective)
                            + " with min obj:"
                            + str(solution["minobjective"])
                        )
                        rxnobj.lower_bound = original_bound
                    else:
                        unneeded.append([rxn_id, solution[key][rxn_id], key])
                        logger.debug(
                            rxn_id
                            + solution[key][rxn_id]
                            + " not needed:"
                            + str(objective)
                        )
        if keep_changes:
            for items in unneeded:
                del solution[items[2]][items[0]]
        return unneeded

    def add_gapfilling(self, solution):
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
        if condition["is_max_threshold"]:
            model.objective.direction = "max"
        else:
            model.objective.direction = "min"
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
                ondition["media"].id
                + "testing leads to infeasible problem. LP file printed to debug!"
            )
            return False
        if value >= condition["threshold"] and condition["is_max_threshold"]:
            # logger.debug("Failed high:"+condition["media"].id+":"+str(new_objective)+";"+str(condition["threshold"]))
            return False
        elif value <= condition["threshold"] and not condition["is_max_threshold"]:
            # logger.debug("Failed low:"+condition["media"].id+":"+str(new_objective)+";"+str(condition["threshold"]))
            return False
        self.test_objective = new_objective
        return True

    def test_condition_list(self, condition_list, model=None):
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
            if not self.test_single_condition(condition, True, model):
                return False
        return True

    def linear_expansion_test(self, reaction_list, condition, currmodel):
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
        if self.test_single_condition(condition, False, currmodel):
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
                if not self.test_single_condition(condition, False, currmodel):
                    # logger.debug(item[0].id+":"+item[1])
                    item[0].upper_bound = 0
                    if item not in filtered_list:
                        item.append(original_bound[count])
                        item.append(self.score)
                        filtered_list.append(item)
            else:
                item[0].lower_bound = original_bound[count]
                if not self.test_single_condition(condition, False, currmodel):
                    # logger.debug(item[0].id+":"+item[1])
                    item[0].lower_bound = 0
                    if item not in filtered_list:
                        item.append(original_bound[count])
                        item.append(self.score)
                        filtered_list.append(item)
            count += 1
        return filtered_list

    def binary_expansion_test(self, reaction_list, condition, currmodel, depth=0):
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
        if self.test_single_condition(condition, False, currmodel):
            return []
        # Check if input list contains only one reaction:
        if len(reaction_list) == 1:
            if reaction_list[0][1] == ">":
                reaction_list[0].append(reaction_list[0][0].upper_bound)
                reaction_list[0][0].upper_bound = 0
            else:
                reaction_list[0].append(reaction_list[0][0].lower_bound)
                reaction_list[0][0].lower_bound = 0
            reaction_list[0].append(self.score)
            filtered_list.append(reaction_list[0])
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
            sub_lists[0], condition, currmodel, newdepth
        )
        for item in new_filter:
            filtered_list.append(item)
        # Submitting second half of reactions for testing - now only breaking reactions are removed from the first list
        for i, item in enumerate(reaction_list):
            if i >= midway_point:
                if item[1] == ">":
                    item[0].upper_bound = original_bound[i]
                else:
                    item[0].lower_bound = original_bound[i]
        new_filter = self.binary_expansion_test(
            sub_lists[1], condition, currmodel, newdepth
        )
        for item in new_filter:
            filtered_list.append(item)
        return filtered_list

    def reaction_expansion_test(
        self, reaction_list, condition_list, binary_search=True
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
        logger.debug("Expansion started!")
        filtered_list = []
        for condition in condition_list:

            logger.debug(f"testing condition {condition}")

            currmodel = self.model
            tic = time.perf_counter()
            new_filtered = []
            with currmodel:
                self.apply_test_condition(condition)
                if binary_search:
                    new_filtered = self.binary_expansion_test(
                        reaction_list, condition, currmodel
                    )
                    for item in new_filtered:
                        filtered_list.append(item)
                else:
                    new_filtered = self.linear_expansion_test(
                        reaction_list, condition, currmodel
                    )
                    for item in new_filtered:
                        filtered_list.append(item)
            # Restoring knockout of newly filtered reactions, which expire after exiting the "with" block above
            for item in new_filtered:
                if item[1] == ">":
                    item[0].upper_bound = 0
                else:
                    item[0].lower_bound = 0
            toc = time.perf_counter()
            logger.debug(
                "Expansion time:" + condition["media"].id + ":" + str((toc - tic))
            )
            logger.debug(
                "Filtered count:"
                + str(len(filtered_list))
                + " out of "
                + str(len(reaction_list))
            )
        return filtered_list

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
        if re.search('(.+)_([a-z]+)(\d*)$', object.id) != None:
            m = re.search('(.+)_([a-z]+)(\d*)$', object.id)
            return (m[1],m[2],m[3])
        return None

    def add_kbase_media(self, kbase_media):
        exIDs = [exRXN.id for exRXN in self.exchange_list()]
        self.model.medium = {"EX_"+exID: -bound[0] for exID, bound in kbase_media.get_media_constraints().items()
                             if "EX_"+exID in exIDs}
        return self.model.medium

    def add_medium(self, media):
        # add the new media and its flux constraints
        exIDs = [exRXN.id for exRXN in self.exchange_list()]
        self.model.medium = {ex: uptake for ex, uptake in media.items() if ex in exIDs}
        return self.model.medium
