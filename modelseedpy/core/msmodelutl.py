import logging
import re
from cobra import Model, Reaction, Metabolite

logger = logging.getLogger(__name__)

def search_name(name):
    name = name.lower()
    name = re.sub(r'_[a-z]\d*$', '', name)
    name = re.sub(r'\W+', '', name)
    return name

class MSModelUtil:

    def __init__(self,model):
        self.model = model
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