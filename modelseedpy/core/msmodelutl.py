# -*- coding: utf-8 -*-
import logging
import re
import time
import json
import sys
import pandas as pd
import cobra
from cobra import Model, Reaction, Metabolite
from cobra.io.json import from_json, to_json
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.exceptions import *
from modelseedpy.core.fbahelper import FBAHelper
from itertools import chain
from optlang.symbolics import Zero
from optlang import Constraint, Objective
from math import isclose
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
        reactants, products = [], []
        for met in stoichiometry:
            stoich = stoichiometry[met]
            if not isinstance(met, str):
                met = (
                    None
                    if FBAHelper.modelseed_id_from_cobra_metabolite(met) == "cpd00067"
                    else met.id
                )
            if met:
                if stoich < 0:
                    reactants.append(met)
                else:
                    products.append(met)
        return [
            "+".join(sorted(reactants)) + "=" + "+".join(sorted(products)),
            "+".join(sorted(products)) + "=" + "+".join(sorted(reactants)),
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


    ########################### CLASS METHODS ###########################

    def __init__(self, model, copy=False, environment=None):
        self.model = model
        if environment is not None:  self.add_medium(environment)
        self.id = model.id
        if copy:
            org_obj_val = model.slim_optimize()
            self.model = model.copy()  ;  self.model.objective = model.objective
            new_obj_val = self.model.slim_optimize()
            if not isclose(org_obj_val, new_obj_val, rel_tol=1e-2) and org_obj_val > 1e-2:
                raise ModelError(f"The {model.id} objective value is corrupted by being copied,"
                                 f" where the original objective value is {org_obj_val}"
                                 f" and the new objective value is {new_obj_val}.")
        self.pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
        self.wsid = None
        self.atputl = None
        self.gfutl = None
        self.metabolite_hash = None
        self.search_metabolite_hash = None
        self.test_objective = None
        self.reaction_scores = None
        self.score = None
        self.integrated_gapfillings = []
        self.attributes = {}
        if hasattr(self.model, "computed_attributes"):
            if self.model.computed_attributes:
                self.attributes = self.model.computed_attributes
        if "pathways" not in self.attributes:
            self.attributes["pathways"] = {}
        if "auxotrophy" not in self.attributes:
            self.attributes["auxotrophy"] = {}
        if "fbas" not in self.attributes:
            self.attributes["fbas"] = {}

    def compute_automated_reaction_scores(self):
        """
        Computes reaction scores automatically from model data
        :return:
        """
        self.reaction_scores = {}

    def add_timeout(self, timeout_s=10):
        from optlang.interface import Configuration
        Configuration(self.model.problem, timeout=timeout_s)

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
            if not compartment:  return self.metabolite_hash[name]
            for met in self.metabolite_hash[name]:
                array = met.id.split("_")
                if array[1] == compartment or met.compartment == compartment:  return [met]
            return None
        sname = MSModelUtil.search_name(name)
        if sname in self.search_metabolite_hash:
            if not compartment:  return self.search_metabolite_hash[sname]
            for met in self.search_metabolite_hash[sname]:
                array = met.id.split("_")
                if array[1] == compartment or met.compartment == compartment:  return [met]
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
                if msid not in output:  output[msid] = []
                output[msid].append(cpd)
        return output

    def exchange_list(self):
        return [rxn for rxn in self.model.reactions if 'EX_' in rxn.id]
    
    def internal_list(self):
        exchanges, transports = self.exchange_list(), self.transport_list()
        return [rxn for rxn in self.model.reactions if rxn not in exchanges and rxn not in transports]

    def transport_list(self):
        all_transports = [rxn for rxn in self.model.reactions if len(set([
            met.id.split("_")[0] for met in rxn.reactants]).intersection(set([
            met.id.split("_")[0] for met in rxn.products]))) > 0]
        # TODO look for compounds that have compounds in different compartments
        # TODO PTS transporters would fail this logic
        # remove biomass reactions
        for rxn in all_transports:
            if "cpd11416" in [met.id.split("_")[0] for met in rxn.metabolites]:  all_transports.remove(rxn)
        return all_transports

    def carbon_mets(self):
        return [met for met in self.model.metabolites if 'C' in met.elements]

    def carbon_exchange_list(self, include_unknown=True):
        if not include_unknown:
            return [ex for ex in self.exchange_list() if "C" in ex.reactants[0].elements]
        return [ex for ex in self.exchange_list() if not ex.reactants[0].elements or "C" in ex.reactants[0].elements]

    def carbon_exchange_mets_list(self, include_unknown=True):
        return self.metabolites_set(self.carbon_exchange_list(include_unknown))

    def exchange_mets_list(self):
        return self.metabolites_set(self.exchange_list())

    def media_exchanges_list(self):
        return [exRXN for exRXN in self.exchange_list() if exRXN.id in self.model.medium]

    def metabolites_set(self, reactions_set=None, ids=False):
        rxns = reactions_set or self.model.reactions
        if ids:  return {met.id for rxn in rxns for met in rxn.metabolites}
        return {met for rxn in rxns for met in rxn.metabolites}

    def bio_rxns_list(self):
        return [rxn for rxn in self.model.reactions if re.search(r"(^bio\d+)", rxn.id)]

    def compatibilize(self, conflicts_file_name="orig_conflicts.json", printing=False):
        from commscores import GEMCompatibility
        self.model = GEMCompatibility.standardize(
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
                    logger.warning("Nonstandard exchange reaction ignored:" + ex_rxn.id)
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
            logger.info(
                "Setting FBAModel computed_attributes to mdlutl attributes"
            )
            self.attributes["gene_count"] = len(self.model.genes)
            self.model.computed_attributes = self.attributes

    def add_ms_reaction(self, rxn_dict, msdb_path=None, msdb_object=None, comp_trans=["c0", "e0"]):
        if msdb_object:  modelseed = msdb_object
        else:
            # from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem
            from modelseedpy.biochem import from_local
            # modelseed = ModelSEEDBiochem.get()
            modelseed = from_local(msdb_path)
        output = []
        model_mets = self.metabolites_set(ids=True)
        for rxnid, comp in rxn_dict.items():
            fullid = f"{rxnid}_{comp}"
            rxn = modelseed.reactions.get_by_id(rxnid)
            new_reaction = Reaction(id=fullid, name=f"{rxn.name}_{comp}")
            metabolites_to_add = {}
            for met, stoich in rxn.metabolites.items():
                comp_num = FBAHelper.compartment_index(met.id)
                if comp_num > len(comp_trans):
                    logger.critical(f"The compartment index {comp_num} is out of range")
                comp_str = comp_trans[comp_num]
                met_output = self.find_met(met.id, comp_str)
                new_met = Metabolite(f"{met.id}_{comp_str}", name=f"{met.name}_{comp_str}",
                                     compartment=comp_str) if not met_output else met_output[0]
                metabolites_to_add[new_met] = stoich
                if new_met.id not in model_mets:  self.model.add_metabolites([new_met])
            new_reaction.add_metabolites(metabolites_to_add)
            output.append(new_reaction)
            print(f"The {new_reaction.id} reaction is defined.")
        self.model.add_reactions(output)
        print(f"{len(output)} reactions and {len(self.model.metabolites)-len(model_mets)} metabolites"
              f" were added to the model.")
        return output

    def create_constraint(self, constraint, coef=None, sloppy=False, printing=False):
        if printing:   print(coef)
        self.model.add_cons_vars(constraint, sloppy=sloppy)
        self.model.solver.update()
        if coef:   constraint.set_linear_coefficients(coef)
        self.model.solver.update()

            # self.model.solver.update()
            # for cons in self.model.constraints:
            #     if cons.name == constraint.name:
            #         cons.set_linear_coefficients(coef)
            #         self.model.solver.update()
        
        # self.model.add_cons_vars(constraint, sloppy=sloppy)
        # self.model.solver.update()

    def add_cons_vars(self, vars_cons, sloppy=False):
        self.model.add_cons_vars(vars_cons, sloppy=sloppy)
        self.model.solver.update()

    def remove_cons_vars(self, vars_cons):
        self.model.remove_cons_vars(vars_cons)
        self.model.solver.update()

    def add_objective(self, objective, direction="max", coef=None):
        self.model.objective = Objective(objective, direction=direction)
        self.model.solver.update()
        if coef:
            self.model.objective.set_linear_coefficients(coef)
            self.model.solver.update()

    def set_objective_from_target_reaction(self, target_reaction, minimize=False):
        target_reaction = self.model.reactions.get_by_id(target_reaction)
        sense = "max" if not minimize else "min"
        self.model.objective = self.model.problem.Objective(
            target_reaction.flux_expression, direction=sense
        )
        return target_reaction

    def biomass_expression(self):
        for met in self.model.metabolites:
            if "cpd11416" in met.id:
                # returns the biomass expression of the lowest cytoplasmic compartment
                return met.constraint.expression

    def add_minimal_objective_cons(self, min_value=0.1, objective_expr=None):
        if "min_value" not in self.model.constraints:
            objective_expr = objective_expr or self.model.objective.expression
            self.create_constraint(Constraint(objective_expr, lb=min_value, ub=None, name="min_value"))
            # print(self.model.constraints["min_value"])
        else:
            print(f"The min_value constraint already exists in {self.model.id}, "
                  f"hence the lb is simply updated from"
                  f" {self.model.constraints['min_value'].lb} to {min_value}.\n")
            self.model.constraints["min_value"].lb = min_value

    def add_exchange_to_model(self, cpd, rxnID):
        self.model.add_boundary(metabolite=Metabolite(id=cpd.id, name=cpd.name, compartment="e0"),
                                reaction_id=rxnID, type="exchange", lb=cpd.minFlux, ub=cpd.maxFlux)

    def update_model_media(self, media):
        medium = self.model.medium
        model_reactions = [rxn.id for rxn in self.model.reactions]
        for cpd in media.data["mediacompounds"]:
            ex_rxn = f"EX_{cpd.id}_e0"
            if ex_rxn not in model_reactions:
                self.add_exchange_to_model(cpd, ex_rxn)
            medium[ex_rxn] = cpd.maxFlux
        self.model.medium = medium
        return self.model

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
    # adding gapfilling compounds to a KBase model saves gapfilled models
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
        removed_rxns = []
        tempmodel = self.model
        if not keep_changes:
            tempmodel = from_json(to_json(self.model))
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
                        logger.info(
                            rxn_id
                            + solution[key][rxn_id]
                            + " needed:"
                            + str(objective)
                            + " with min obj:"
                            + str(solution["minobjective"])
                        )
                        rxnobj.upper_bound = original_bound
                    else:
                        removed_rxns.append(rxnobj)
                        unneeded.append([rxn_id, solution[key][rxn_id], key])
                        logger.info(
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
                        logger.info(
                            rxn_id
                            + solution[key][rxn_id]
                            + " needed:"
                            + str(objective)
                            + " with min obj:"
                            + str(solution["minobjective"])
                        )
                        rxnobj.lower_bound = original_bound
                    else:
                        removed_rxns.append(rxnobj)
                        unneeded.append([rxn_id, solution[key][rxn_id], key])
                        logger.info(
                            rxn_id
                            + solution[key][rxn_id]
                            + " not needed:"
                            + str(objective)
                        )
        if keep_changes:
            tempmodel.remove_reactions(removed_rxns)
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
    def run_fba(self, media=None, pfba=False, fva_reactions=None):
        from cobra import flux_analysis
        if media:
            self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        if pfba:
            return flux_analysis.pfba(self.model)
        if fva_reactions is not None:
            return flux_analysis.variability.flux_variability_analysis(self.model, fva_reactions)
        return self.model.optimize()

    def resource_balance_constraint(self, flux_limit=140):
        vars_coef = {}
        for rxn in self.model.reactions:
            if "EX_" not in rxn.id:
                vars_coef[rxn.forward_variable] = vars_coef[rxn.reverse_variable] = 1
        self.create_constraint(Constraint(Zero, lb=0, ub=flux_limit, name="resource_balance_limit"), coef=vars_coef)

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
            print("applying - bad")
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
        logger.debug(f"Expansion started! Binary = {binary_search}")
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
                        if item not in filtered_list:
                            filtered_list.append(item)
                else:
                    new_filtered = self.linear_expansion_test(
                        reaction_list, condition, currmodel
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
            gf_filter_att = self.get_attributes("gf_filter", {})
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
            gf_filter_att = self.save_attributes(gf_filter_att, "gf_filter")
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
                rxnobj = tempmodel.reactions.get_by_id(item[0])
                if item[1] == ">":
                    original_bound = rxnobj.upper_bound
                    rxnobj.upper_bound = 0
                    if item[0] not in output:
                        output[item[0]] = {}
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
                    if item[0] not in output:
                        output[item[0]] = {}
                    output[item[0]][item[1]] = self.run_biomass_dependency_test(
                        target_rxn_obj,
                        tempmodel,
                        original_objective,
                        min_flex_obj,
                        rxn_list,
                    )
                    rxnobj.lower_bound = original_bound
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
        stoichiometry = {}
        id_hash = self.msid_hash()
        for msid, content in coefs.items():
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

    def costless_excreta(self, pfba=False):
        # the double optimization is intended to truly find the maximal biomass growth
        original_objective = self.model.objective
        minObj_cons = Constraint(self.model.objective.expression, lb=self.model.slim_optimize(), name="minObj")
        self.add_cons_vars([minObj_cons])
        if pfba:
            self.model.problem.constraints.minObj_cons.lb = self.model.slim_optimize()
            reaction_variables = ((rxn.forward_variable, rxn.reverse_variable) for rxn in self.model.reactions)
            self.model.problem.Objective(sum(chain(*reaction_variables)), direction="min")
        sol = self.model.optimize()
        # revert conditions to before the simulation
        self.model.objective = original_objective
        self.remove_cons_vars([minObj_cons])
        return [rxnID.replace("EX_", "").replace("_e0", "") for rxnID, flux in sol.fluxes.items()
                if "EX_" in rxnID and flux > 0]

    @staticmethod
    def parse_id(cobra_obj):
        MSID = re.search("(.+)_([a-z])(\d+)$", cobra_obj.id)
        if MSID is not None:  return (MSID[1], MSID[2], int(MSID[3]))
        nonMSID = re.search("(.+)\[([a-z])\]$", cobra_obj.id)
        if nonMSID is not None:  return (nonMSID[1], nonMSID[2])
        return (cobra_obj.id.replace("EX_", ""), "c" if "EX_" not in cobra_obj.id else "e")

    def add_kbase_media(self, kbase_media):
        exIDs = [exRXN.id for exRXN in self.exchange_list()]
        self.model.medium = {"EX_"+exID: -bound[0] for exID, bound in kbase_media.get_media_constraints().items()
                             if "EX_"+exID in exIDs}
        return self.model.medium

    def add_medium(self, media, uniform_uptake=None):
        # add the new media and its flux constraints
        if media is None:  return self.model.medium
        exIDs = [exRXN.id for exRXN in self.exchange_list()]
        if not hasattr(media, "items"):  media = FBAHelper.convert_kbase_media(media)
        self.model.medium = {ex: uptake for ex, uptake in media.items() if ex in exIDs}
        if uniform_uptake is not None:  self.model.medium = dict(zip(
            list(self.model.medium.keys()), [uniform_uptake]*len(self.model.medium)))
        return self.model.medium