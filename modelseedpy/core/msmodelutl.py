import logging
import re
from cobra import Model, Reaction, Metabolite
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

logger = logging.getLogger(__name__)

def search_name(name):
    name = name.lower()
    name = re.sub(r'_[a-z]\d*$', '', name)
    name = re.sub(r'\W+', '', name)
    return name

class MSModelUtil:

    def __init__(self,model):
        self.model = model
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.metabolite_hash = None
        self.search_metabolite_hash = None
    
    def build_metabolite_hash(self):
        self.metabolite_hash = {}
        self.search_metabolite_hash = {}
        for met in self.model.metabolites:
            self.add_name_to_metabolite_hash(met.id,met)
            self.add_name_to_metabolite_hash(met.name,met)
            for anno in met.annotation:
                if isinstance(met.annotation[anno], list):
                    for item in met.annotation[anno]:
                        self.add_name_to_metabolite_hash(item,met)
                else:
                    self.add_name_to_metabolite_hash(met.annotation[anno],met)
    
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
        exchange_reactions = []
        for reaction in self.model.reactions:
            if reaction.id[:3] == 'EX_':
                exchange_reactions.append(reaction)
        return exchange_reactions
    
    def exchange_hash(self):
        exchange_reactions = {}
        exlist = self.exchange_list()
        for reaction in exlist:
            for met in reaction.metabolites:
                if reaction.metabolites[met] == -1:
                    exchange_reactions[met] = reaction
                else:
                    logger.warn("Nonstandard exchange reaction ignored:"+reaction.id)
        return exchange_reactions
    
    def add_missing_exchanges(self,media):
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
                        #We prefer to add a transport for the cytosol compound
                        cpd = met
                if cpd == None:
                    #No cytosol compound exists so choosing the first version we found that does exist
                    cpd = mets[0]
                if found == 0:
                    #No transporter currently exists - adding exchange reaction for the compound that does exist
                    output.append(cpd.id)
                    exchange_list.append(cpd)
        if len(exchange_list) > 0:
            self.add_exchanges_for_metabolites(exchange_list)
        return output
        
    def add_exchanges_for_metabolites(self,cpds,uptake=0,excretion=0,prefix='EX_', prefix_name='Exchange for '):
        drains = []
        for cpd in cpds:
            drain_reaction = Reaction(id=f'{prefix}{cpd.id}',
                                      name=prefix_name + cpd.name,
                                      lower_bound=-1*uptake, 
                                      upper_bound=excretion)
            drain_reaction.add_metabolites({cpd : -1})
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
            drains.append(drain_reaction)
        self.model.add_reactions(drains)
        return drains
        
    def reaction_scores(self):
        return {}
    
    #Required this function to add gapfilled compounds to a KBase model for saving gapfilled model    
    def convert_cobra_compound_to_kbcompound(self,cpd,kbmodel,add_to_model = 1):
        refid = "cpd00000"
        if re.search('cpd\d+_[a-z]+',cpd.id):
            refid = cpd.id
            refid = re.sub("_[a-z]\d+$","",refid)
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

    #Required this function to add gapfilled reactions to a KBase model for saving gapfilled model    
    def convert_cobra_reaction_to_kbreaction(self,rxn,kbmodel,cpd_hash,direction = "=",add_to_model = 1,reaction_genes = None):
        rxnref = "~/template/reactions/id/rxn00000_c"
        if re.search('rxn\d+_[a-z]+',rxn.id):
            rxnref = "~/template/reactions/id/"+rxn.id
            rxnref = re.sub("\d+$","",rxnref)
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
        for cpd in rxn.metabolites:
            if cpd.id not in kbmodel["modelcompounds"]:
                cpd_hash[cpd.id] = self.convert_cobra_compound_to_kbcompound(cpd,kbmodel,1)
            rxn_data["modelReactionReagents"].append({
                "coefficient" : rxn.metabolites[cpd],
                "modelcompound_ref" : "~/modelcompounds/id/"+cpd.id
            })
        if reaction_genes != None and rxn.id in reaction_genes:
            best_gene = None
            for gene in reaction_genes[rxn.id]:
                if best_gene == None or reaction_genes[rxn.id][gene] > reaction_genes[rxn.id][best_gene]:
                    best_gene = gene
            rxn_data["modelReactionProteins"] = [{"note":"Added from gapfilling","modelReactionProteinSubunits":[],"source":"Unknown"}]
            rxn_data["modelReactionProteins"][0]["modelReactionProteinSubunits"] = [{"note":"Added from gapfilling","optionalSubunit":0,"triggering":1,"feature_refs":["~/genome/features/id/"+best_gene],"role":"Unknown"}]
        if add_to_model == 1:
            kbmodel["modelreactions"].append(rxn_data)
        return rxn_data
    
    def add_gapfilling_solution_to_kbase_model(self,newmodel,gapfilled_reactions,gfid=None,media_ref = None,reaction_genes = None):
        rxn_table = []
        gapfilling_obj = None
        if gfid == None:
            largest_index = 0
            for gapfilling in newmodel["gapfillings"]:
                current_index = int(gapfilling["id"].split(".").pop())
                if largest_index == 0 or largest_index < current_index:
                    largest_index = current_index
            largest_index += 1
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
            kbrxn = self.convert_cobra_reaction_to_kbreaction(reaction,newmodel,cpd_hash,gapfilled_reactions["new"][rxn],1,reaction_genes)
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