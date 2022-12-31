import logging
import re
from cobra import Model, Reaction, Metabolite  # !!! Model and Metabolite are not used
from modelseedpy.core.exceptions import ModelError
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.fbahelper import FBAHelper
from optlang.symbolics import Zero
from optlang import Constraint
import time

logger = logging.getLogger(__name__)

def metabolite_msid(metabolite):
    if re.search('^(cpd\d+)', metabolite.id):
        m = re.search('^(cpd\d+)', metabolite.id)
        return m[1]
    for anno in metabolite.annotation:
        if isinstance(metabolite.annotation[anno], list):
            for item in metabolite.annotation[anno]:
                if re.search('^(cpd\d+)', item):
                    m = re.search('^(cpd\d+)', item)
                    return m[1]
        elif re.search('^(cpd\d+)', metabolite.annotation[anno]):
            m = re.search('^(cpd\d+)', metabolite.annotation[anno])
            return m[1]
    return None
    
def reaction_msid(reaction):
    if re.search('^(rxn\d+)', reaction.id):
        m = re.search('^(rxn\d+)', reaction.id)
        return m[1]
    for anno in reaction.annotation:
        if isinstance(reaction.annotation[anno], list):
            for item in reaction.annotation[anno]:
                if re.search('^(rxn\d+)', item):
                    m = re.search('^(rxn\d+)', item)
                    return m[1]
        elif re.search('^(rxn\d+)', reaction.annotation[anno]):
            m = re.search('^(rxn\d+)', reaction.annotation[anno])
            return m[1]
    return None

def stoichiometry_to_string(stoichiometry):
    reactants = []
    products = []
    for met in stoichiometry:
        coef = stoichiometry[met]
        if not isinstance(met, str):
            if metabolite_msid(met) == "cpd00067":
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
    return ["+".join(reactants)+"="+"+".join(products),"+".join(products)+"="+"+".join(reactants)]

def search_name(name):
    name = name.lower()
    name = re.sub(r'_[a-z]\d*$', '', name)
    name = re.sub(r'\W+', '', name)
    return name


class MSModelUtil:

    def __init__(self,model):
        self.model = model.copy()
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.atputl = self.gfutl = self.metabolite_hash = self.search_metabolite_hash = None
        self.test_objective = self.score = None
        if self.model.slim_optimize() == 0:
            raise ModelError(f"The {self.model.id} model is broken and does not grow in complete media.")
    
    def printlp(self,lpfilename="debug.lp"):
        with open(lpfilename, 'w') as out:
            out.write(str(self.model.solver))

    def build_metabolite_hash(self):
        self.metabolite_hash = {}
        self.search_metabolite_hash = {}
        for met in self.model.metabolites:
            self.add_name_to_metabolite_hash(met.id,met)
            self.add_name_to_metabolite_hash(met.name,met)
            for item in met.annotation.values():
                if isinstance(item, list):
                    for entry in item:
                        self.add_name_to_metabolite_hash(entry,met)
                else:
                    self.add_name_to_metabolite_hash(item,met)

    def add_name_to_metabolite_hash(self,name,met):
        if name not in self.metabolite_hash:
            self.metabolite_hash[name] = []
        self.metabolite_hash[name].append(met)
        sname = search_name(name)
        if sname not in self.search_metabolite_hash:
            self.search_metabolite_hash[sname] = []
        self.search_metabolite_hash[sname].append(met)
        
    def find_met(self,name):
        if self.metabolite_hash == None:
            self.build_metabolite_hash()
        if name in self.metabolite_hash:
            return self.metabolite_hash[name]
        sname = search_name(name)
        if sname in self.search_metabolite_hash:
            return self.search_metabolite_hash[sname]
        logger.info(name," not found in model!")
        return []
    
    def rxn_hash(self):
        output = {}
        for rxn in self.model.reactions:
            strings = stoichiometry_to_string(rxn.metabolites)
            output[strings[0]] = [rxn,1]
            output[strings[1]] = [rxn,-1]
        return output
    
    def find_reaction(self,stoichiometry):
        output = stoichiometry_to_string(stoichiometry)
        atpstring = output[0]
        rxn_hash = self.rxn_hash()
        if atpstring in rxn_hash:
            return rxn_hash[atpstring]
        return None
    
    def msid_hash(self): 
        output = {}
        for cpd in self.model.metabolites:
            msid = metabolite_msid(cpd)
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

    def exchange_hash(self):
        exchange_reactions = {}
        for ex_rxn in self.exchange_list():
            for met in ex_rxn.metabolites:
                if ex_rxn.metabolites[met] == -1:
                    exchange_reactions[met] = ex_rxn
                else:
                    logger.warn("Nonstandard exchange reaction ignored:"+ex_rxn.id)
        return exchange_reactions

    def var_names_list(self):
        return [var.name for var in self.model.variables]
    
    def add_missing_exchanges(self,media):
        output, exchange_list = [], []
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
                        #We prefer to add a transport for the cytosol compound
                        cpd = met
                if cpd == None:
                    #No cytosol compound exists so choosing the first version we found that does exist
                    cpd = mets[0]
                if found:
                    #No transporter currently exists - adding exchange reaction for the compound that does exist
                    output.append(cpd.id)
                    exchange_list.append(cpd)
        if len(exchange_list) > 0:
            self.add_exchanges_for_metabolites(exchange_list)
        return output
        
    def add_exchanges_for_metabolites(self,cpds,uptake=0,excretion=0,prefix='EX_', prefix_name='Exchange for '):  # !!! why not use the existing COBRApy infrastructure of add_boundary() ?
        drains = []
        for cpd in cpds:
            drain_reaction = Reaction(
                id=f'{prefix}{cpd.id}', name=prefix_name + cpd.name, lower_bound=-uptake, upper_bound=excretion)
            drain_reaction.add_metabolites({cpd : -1})
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
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
        cpd_data = {
            "aliases": [],
            "charge": cpd.charge,
            "compound_ref": "~/template/compounds/id/"+refid,
            "dblinks": {},
            "formula": cpd.formula,
            "id": cpd.id,
            "modelcompartment_ref": "~/modelcompartments/id/"+cpd.id.split("_").pop(),
            "name": cpd.name,
            "numerical_attributes": {},
            "string_attributes": {}
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
                'reverse': rxn.reverse_variable.primal,
                'forward': rxn.forward_variable.primal
            }
        return flux_values

    #Required this function to add gapfilled reactions to a KBase model for saving gapfilled model    
    def convert_cobra_reaction_to_kbreaction(self, rxn, kbmodel=None, cpd_hash={}, direction="=", add_to_model=1, reaction_genes=None):
        rxnref = "~/template/reactions/id/rxn00000_c"
        if re.search('rxn\d+_[a-z]+',rxn.id):
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
            "modelcompartment_ref": "~/modelcompartments/id/"+rxn.id.split("_").pop(),
            "name": rxn.name,
            "numerical_attributes": {},
            "probability": 0,
            "protons": 0,
            "reaction_ref": rxnref,
            "string_attributes": {}
        }
        for met in rxn.metabolites:
            if met.id not in kbmodel["modelcompounds"]:
                cpd_hash[met.id] = self.convert_cobra_compound_to_kbcompound(met,kbmodel)
            rxn_data["modelReactionReagents"].append({
                "coefficient" : rxn.metabolites[met],
                "modelcompound_ref" : "~/modelcompounds/id/"+met.id
            })
        if reaction_genes and rxn.id in reaction_genes:
            best_gene = None
            for gene, val in reaction_genes[rxn.id].items():
                if not best_gene or val > reaction_genes[rxn.id][best_gene]:
                    best_gene = gene
            rxn_data["modelReactionProteins"] = [{"note":"Added from gapfilling","modelReactionProteinSubunits":[],"source":"Unknown"}]
            rxn_data["modelReactionProteins"][0]["modelReactionProteinSubunits"] = [
                {"note":"Added from gapfilling","optionalSubunit":0,"triggering":1,"feature_refs":["~/genome/features/id/"+best_gene],"role":"Unknown"}]
        if add_to_model:
            kbmodel["modelreactions"].append(rxn_data)
        return rxn_data
    
    def add_gapfilling_solution_to_kbase_model(self,newmodel,gapfilled_reactions,gfid=None,media_ref = None,reaction_genes = None):
        """
        NOTE: to be moved to cobrakbase
        """
        rxn_table = []
        gapfilling_obj = None
        if gfid is None:
            largest_index = max([int(gapfilling["id"].split(".").pop()) for gapfilling in newmodel["gapfillings"]]) + 1
            gfid = "gf."+str(largest_index)
        else:
            for gapfilling in newmodel["gapfillings"]:
                if gapfilling["id"] == gfid:
                    gapfilling_obj = gapfilling
        if not gapfilling_obj:
            gapfilling_obj = {
                "gapfill_id": newmodel["id"]+"."+gfid,
                "id": gfid,
                "integrated": 1,
                "integrated_solution": "0",
                "media_ref": media_ref
            }
            newmodel["gapfillings"].append(gapfilling_obj)
        cpd_hash = {}
        for cpd in newmodel["modelcompounds"]:
            cpd_hash[cpd["id"]] = cpd
        for rxn in gapfilled_reactions["new"]:
            reaction = self.model.reactions.get_by_id(rxn)
            kbrxn = self.convert_cobra_reaction_to_kbreaction(reaction,newmodel,cpd_hash,gapfilled_reactions["new"][rxn],reaction_genes)
            kbrxn["gapfill_data"][gfid] = dict()
            kbrxn["gapfill_data"][gfid]["0"] = [gapfilled_reactions["new"][rxn],1,[]]
            #rxn_table.append({
            #    'id':kbrxn["id"],
            #    'name':kbrxn["name"],
            #    'direction':format_direction(kbrxn["direction"]),
            #    'gene':format_gpr(kbrxn),
            #    'equation':format_equation(kbrxn,cpd_hash),
            #    'newrxn':1
            #})
        for rxn in gapfilled_reactions["reversed"]:
            for kbrxn in newmodel["modelreactions"]:
                if kbrxn["id"] == rxn:
                    kbrxn["direction"] = "="
                    #rxn_table.append({
                    #    'id':kbrxn["id"],
                    #    'name':kbrxn["name"],
                    #    'direction':format_direction(kbrxn["direction"]),
                    #    'gene':format_gpr(kbrxn),
                    #    'equation':format_equation(kbrxn,cpd_hash),
                    #    'newrxn':0
                    #})
                    kbrxn["gapfill_data"][gfid] = dict()
                    kbrxn["gapfill_data"][gfid]["0"] = [gapfilled_reactions["reversed"][rxn],1,[]]
        return rxn_table

    def apply_test_condition(self,condition,model = None):
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
        if model == None:
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
    
    def test_single_condition(self,condition,apply_condition=True,model=None):
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
        if model == None:
            model = self.model
        if apply_condition:
            self.apply_test_condition(condition,model)
        new_objective = model.slim_optimize()
        value = new_objective
        if "change" in condition and condition["change"]:
            if self.test_objective != None:
                value = new_objective - self.test_objective
        self.score = value
        if model.solver.status != 'optimal':
            self.printlp("Infeasible.lp")
            logger.critical("Infeasible problem - LP file printed to debug!")
            return False
        if value >= condition["threshold"] and condition["is_max_threshold"]:
            logger.debug("Failed high:"+str(self.test_objective)+";"+str(condition["threshold"]))
            return False
        elif value <= condition["threshold"] and not condition["is_max_threshold"]:
            logger.debug("Failed low:"+str(self.test_objective)+";"+str(condition["threshold"]))
            return False
        self.test_objective = new_objective
        return True
    
    def test_condition_list(self,condition_list,model=None):
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
            if not self.test_single_condition(condition,True,model):
                return False
        return True
    
    def reaction_expansion_test(self,reaction_list,condition_list):
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
        print("Expansion started!")
        tic = time.perf_counter()
        filtered_list = []
        for condition in condition_list:
            currmodel = self.model
            with currmodel:
                self.apply_test_condition(condition)
                # First knockout all reactions in the input list and save original bounds
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
                        if not self.test_single_condition(condition,False,currmodel):
                            item[0].upper_bound = 0
                            if item not in filtered_list:
                                item.append(original_bound[count])
                                item.append(self.score)
                                filtered_list.append(item)
                    else:
                        item[0].lower_bound = original_bound[count]
                        if not self.test_single_condition(condition,False,currmodel):
                            item[0].lower_bound = 0
                            if item not in filtered_list:
                                item.append(original_bound[count])
                                item.append(self.score)
                                filtered_list.append(item)
                    count += 1
        toc = time.perf_counter()
        print("Expansion time:",(toc-tic))
        print("Filtered count:",len(filtered_list)," out of ",len(reaction_list))
        return filtered_list

    def add_atp_hydrolysis(self,compartment):
        #Searching for ATP hydrolysis compounds
        coefs = {"cpd00002":[-1,compartment],"cpd00001":[-1,compartment],"cpd00008":[1,compartment],"cpd00009":[1,compartment],"cpd00067":[1,compartment]}
        msids = ["cpd00002","cpd00001","cpd00008","cpd00009","cpd00067"]
        stoichiometry = {}
        id_hash = self.msid_hash()
        for msid in msids:
            if msid not in id_hash:
                logger.warning("Compound "+msid+" not found in model!")
                return None
            else:
                for cpd in id_hash[msid]:
                    if cpd.compartment == coefs[msid][1]:
                        stoichiometry[cpd] = coefs[msid][0]
        output = self.find_reaction(stoichiometry)
        if output != None and output[1] == ">":
            return {"reaction":output[0],"direction":">","new":False}
        cobra_reaction = Reaction("rxn00062_"+compartment, 
                                  name="ATP hydrolysis", 
                                  lower_bound=0, 
                                  upper_bound=1000)
        cobra_reaction.annotation["sbo"] = "SBO:0000176" #biochemical reaction
        cobra_reaction.annotation["seed.reaction"] = "rxn00062"
        cobra_reaction.add_metabolites(stoichiometry)
        self.model.add_reactions([cobra_reaction])
        return {"reaction":cobra_reaction,"direction":">","new":True}
    
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
