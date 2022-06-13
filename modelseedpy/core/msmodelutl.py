import logging
logger = logging.getLogger(__name__)

from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from cobra import Reaction
import re

def search_name(name):
    name = name.lower()
    name = re.sub(r'_[a-z]\d*$', '', name)
    name = re.sub(r'\W+', '', name)
    return name

class MSModelUtil:

    def __init__(self,model):
        self.model = model
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.metabolite_hash = self.search_metabolite_hash = None
    
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
    
    def exchange_list(self): 
        return [rxn for rxn in self.model.reactions if 'EX_' in rxn.id]
    
    def exchange_hash(self):
        exchange_reactions = {}
        for ex_rxn in self.exchange_list():
            for met in ex_rxn.metabolites:
                if ex_rxn.metabolites[met] == -1:
                    exchange_reactions[met] = ex_rxn
                else:
                    logger.warn("Nonstandard exchange reaction ignored:"+ex_rxn.id)
        return exchange_reactions
    
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
        
    def add_exchanges_for_metabolites(self,cpds,uptake=0,excretion=0,prefix='EX_', prefix_name='Exchange for '):
        drains = []
        for cpd in cpds:
            drain_reaction = Reaction(
                id=f'{prefix}{cpd.id}', name=prefix_name + cpd.name, lower_bound=-uptake, upper_bound=excretion)
            drain_reaction.add_metabolites({cpd : -1})
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
            drains.append(drain_reaction)
        self.model.add_reactions(drains)
        return drains
        
    # def reaction_scores(self):  #!!! Can this be deleted?
    #     return {}
    
    #adding gapfilling compounds to a KBase model saves gapfilled models 
    def convert_cobra_compound_to_kbcompound(self, cpd, kbmodel, add_to_model = True):
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
        if add_to_model:
            kbmodel["modelcompounds"].append(cpd_data)
        return cpd_data

    #adding gapfilling reactions to a KBase model saves gapfilled models   
    def convert_cobra_reaction_to_kbreaction(self,rxn,kbmodel = None,cpd_hash = {},direction = "=",reaction_genes = {}, add_to_model = True):
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
        if rxn.id in reaction_genes:
            best_gene = None
            for gene, val in reaction_genes[rxn.id].items():
                if best_gene == None or val > reaction_genes[rxn.id][best_gene]:
                    best_gene = gene
            rxn_data["modelReactionProteins"] = [{"note":"Added from gapfilling","modelReactionProteinSubunits":[],"source":"Unknown"}]
            rxn_data["modelReactionProteins"][0]["modelReactionProteinSubunits"] = [
                {"note":"Added from gapfilling","optionalSubunit":0,"triggering":1,"feature_refs":["~/genome/features/id/"+best_gene],"role":"Unknown"}]
        if add_to_model:
            kbmodel["modelreactions"].append(rxn_data)
        return rxn_data
    
    def add_gapfilling_solution_to_kbase_model(self,newmodel,gapfilled_reactions,gfid=None,media_ref = None,reaction_genes = None):
        rxn_table = []
        gapfilling_obj = None
        if gfid is None:
            largest_index = max([int(gapfilling["id"].split(".").pop()) for gapfilling in newmodel["gapfillings"]]) + 1
            gfid = "gf."+str(largest_index)
        else:
            for gapfilling in newmodel["gapfillings"]:
                if gapfilling["id"] == gfid:
                    gapfilling_obj = gapfilling
        if gapfilling_obj == None:    
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
            rxn_table.append({
                'id':kbrxn["id"],
                'name':kbrxn["name"],
                'direction':format_direction(kbrxn["direction"]),
                'gene':format_gpr(kbrxn),
                'equation':format_equation(kbrxn,cpd_hash),
                'newrxn':1
            })
        for rxn in gapfilled_reactions["reversed"]:
            for kbrxn in newmodel["modelreactions"]:
                if kbrxn["id"] == rxn:
                    kbrxn["direction"] = "="
                    rxn_table.append({
                        'id':kbrxn["id"],
                        'name':kbrxn["name"],
                        'direction':format_direction(kbrxn["direction"]),
                        'gene':format_gpr(kbrxn),
                        'equation':format_equation(kbrxn,cpd_hash),
                        'newrxn':0
                    })
                    kbrxn["gapfill_data"][gfid] = dict()
                    kbrxn["gapfill_data"][gfid]["0"] = [gapfilled_reactions["reversed"][rxn],1,[]]
        return rxn_table
