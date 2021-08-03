import logging
import itertools
import cobra
from modelseedpy.core import FBAHelper
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

logger = logging.getLogger(__name__)

#Adding a few exception classes to handle different types of errors
class GapfillingError(Exception):
    """Error in model gapfilling"""
    pass

default_blacklist = ["rxn12985", "rxn00238", "rxn07058", "rxn05305", "rxn00154", "rxn09037", "rxn10643",
                     "rxn11317", "rxn05254", "rxn05257", "rxn05258", "rxn05259", "rxn05264", "rxn05268",
                     "rxn05269", "rxn05270", "rxn05271", "rxn05272", "rxn05273", "rxn05274", "rxn05275",
                     "rxn05276", "rxn05277", "rxn05278", "rxn05279", "rxn05280", "rxn05281", "rxn05282",
                     "rxn05283", "rxn05284", "rxn05285", "rxn05286", "rxn05963", "rxn05964", "rxn05971",
                     "rxn05989", "rxn05990", "rxn06041", "rxn06042", "rxn06043", "rxn06044", "rxn06045",
                     "rxn06046", "rxn06079", "rxn06080", "rxn06081", "rxn06086", "rxn06087", "rxn06088",
                     "rxn06089", "rxn06090", "rxn06091", "rxn06092", "rxn06138", "rxn06139", "rxn06140",
                     "rxn06141", "rxn06145", "rxn06217", "rxn06218", "rxn06219", "rxn06220", "rxn06221",
                     "rxn06222", "rxn06223", "rxn06235", "rxn06362", "rxn06368", "rxn06378", "rxn06474",
                     "rxn06475", "rxn06502", "rxn06562", "rxn06569", "rxn06604", "rxn06702", "rxn06706",
                     "rxn06715", "rxn06803", "rxn06811", "rxn06812", "rxn06850", "rxn06901", "rxn06971",
                     "rxn06999", "rxn07123", "rxn07172", "rxn07254", "rxn07255", "rxn07269", "rxn07451",
                     "rxn09037", "rxn10018", "rxn10077", "rxn10096", "rxn10097", "rxn10098", "rxn10099",
                     "rxn10101", "rxn10102", "rxn10103", "rxn10104", "rxn10105", "rxn10106", "rxn10107",
                     "rxn10109", "rxn10111", "rxn10403", "rxn10410", "rxn10416", "rxn11313", "rxn11316",
                     "rxn11318", "rxn11353", "rxn05224", "rxn05795", "rxn05796", "rxn05797", "rxn05798",
                     "rxn05799", "rxn05801", "rxn05802", "rxn05803", "rxn05804", "rxn05805", "rxn05806",
                     "rxn05808", "rxn05812", "rxn05815", "rxn05832", "rxn05836", "rxn05851", "rxn05857",
                     "rxn05869", "rxn05870", "rxn05884", "rxn05888", "rxn05896", "rxn05898", "rxn05900",
                     "rxn05903", "rxn05904", "rxn05905", "rxn05911", "rxn05921", "rxn05925", "rxn05936",
                     "rxn05947", "rxn05956", "rxn05959", "rxn05960", "rxn05980", "rxn05991", "rxn05992",
                     "rxn05999", "rxn06001", "rxn06014", "rxn06017", "rxn06021", "rxn06026", "rxn06027",
                     "rxn06034", "rxn06048", "rxn06052", "rxn06053", "rxn06054", "rxn06057", "rxn06059",
                     "rxn06061", "rxn06102", "rxn06103", "rxn06127", "rxn06128", "rxn06129", "rxn06130",
                     "rxn06131", "rxn06132", "rxn06137", "rxn06146", "rxn06161", "rxn06167", "rxn06172",
                     "rxn06174", "rxn06175", "rxn06187", "rxn06189", "rxn06203", "rxn06204", "rxn06246",
                     "rxn06261", "rxn06265", "rxn06266", "rxn06286", "rxn06291", "rxn06294", "rxn06310",
                     "rxn06320", "rxn06327", "rxn06334", "rxn06337", "rxn06339", "rxn06342", "rxn06343",
                     "rxn06350", "rxn06352", "rxn06358", "rxn06361", "rxn06369", "rxn06380", "rxn06395",
                     "rxn06415", "rxn06419", "rxn06420", "rxn06421", "rxn06423", "rxn06450", "rxn06457",
                     "rxn06463", "rxn06464", "rxn06466", "rxn06471", "rxn06482", "rxn06483", "rxn06486",
                     "rxn06492", "rxn06497", "rxn06498", "rxn06501", "rxn06505", "rxn06506", "rxn06521",
                     "rxn06534", "rxn06580", "rxn06585", "rxn06593", "rxn06609", "rxn06613", "rxn06654",
                     "rxn06667", "rxn06676", "rxn06693", "rxn06730", "rxn06746", "rxn06762", "rxn06779",
                     "rxn06790", "rxn06791", "rxn06792", "rxn06793", "rxn06794", "rxn06795", "rxn06796",
                     "rxn06797", "rxn06821", "rxn06826", "rxn06827", "rxn06829", "rxn06839", "rxn06841",
                     "rxn06842", "rxn06851", "rxn06866", "rxn06867", "rxn06873", "rxn06885", "rxn06891",
                     "rxn06892", "rxn06896", "rxn06938", "rxn06939", "rxn06944", "rxn06951", "rxn06952",
                     "rxn06955", "rxn06957", "rxn06960", "rxn06964", "rxn06965", "rxn07086", "rxn07097",
                     "rxn07103", "rxn07104", "rxn07105", "rxn07106", "rxn07107", "rxn07109", "rxn07119",
                     "rxn07179", "rxn07186", "rxn07187", "rxn07188", "rxn07195", "rxn07196", "rxn07197",
                     "rxn07198", "rxn07201", "rxn07205", "rxn07206", "rxn07210", "rxn07244", "rxn07245",
                     "rxn07253", "rxn07275", "rxn07299", "rxn07302", "rxn07651", "rxn07723", "rxn07736",
                     "rxn07878", "rxn11417", "rxn11582", "rxn11593", "rxn11597", "rxn11615", "rxn11617",
                     "rxn11619", "rxn11620", "rxn11624", "rxn11626", "rxn11638", "rxn11648", "rxn11651",
                     "rxn11665", "rxn11666", "rxn11667", "rxn11698", "rxn11983", "rxn11986", "rxn11994",
                     "rxn12006", "rxn12007", "rxn12014", "rxn12017", "rxn12022", "rxn12160", "rxn12161",
                     "rxn01267", "rxn05294", "rxn04656"]

class MSGapfill:

    def __init__(self,model,default_gapfill_templates = [],default_gapfill_models = [],test_conditions = [],reaction_scores = {},blacklist = []):
        self.auto_sink = ["cpd02701", "cpd11416", "cpd15302"]
        self.model = model
        self.gfmodel = None
        self.model_penalty = 1
        self.default_gapfill_models = default_gapfill_models
        self.default_gapfill_templates = default_gapfill_templates
        self.gapfill_templates_by_index = {}
        self.gapfill_models_by_index = {}
        self.gapfill_all_indecies_with_default_templates = True
        self.gapfill_all_indecies_with_default_models = True
        self.blacklist = default_blacklist
        for rxnid in blacklist:
            if rxnid not in self.blacklist:
                self.blacklist.append(rxnid)
        self.lp_filename = None
        self.test_condition_iteration_limit = 10
        self.test_conditions = test_conditions
        self.reaction_scores = reaction_scores
        self.solutions = {}
        
    def run_gapfilling(self,media = None,target_reaction = "bio1",minimum_obj = 0.01,binary_check = False):
        self.model.objective = target_reaction
        self.gfmodel = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
        pkgmgr = MSPackageManager.get_pkg_mgr(self.gfmodel)
        pkgmgr.getpkg("GapfillingPkg").build_package({
            "auto_sink":self.auto_sink,
            "model_penalty":self.model_penalty,
            "default_gapfill_models":self.default_gapfill_models,
            "default_gapfill_templates":self.default_gapfill_templates,
            "gapfill_templates_by_index":self.gapfill_templates_by_index,
            "gapfill_models_by_index":self.gapfill_models_by_index,
            "gapfill_all_indecies_with_default_templates":self.gapfill_all_indecies_with_default_templates,
            "gapfill_all_indecies_with_default_models":self.gapfill_all_indecies_with_default_models,
            "default_excretion":100,
            "default_uptake":100,
            "minimum_obj":minimum_obj,
            "blacklist":self.blacklist,
            "reaction_scores":self.reaction_scores,
            "set_objective": 1
        })
        pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        if self.lp_filename != None:
            with open(self.lp_filename, 'w') as out:
                out.write(str(self.gfmodel.solver))
        sol = self.gfmodel.optimize()
        if media not in self.solutions:
            self.solutions[media] = {}
        self.solutions[media][target_reaction] = pkgmgr.getpkg("GapfillingPkg").compute_gapfilled_solution()
        if self.test_conditions != None:
            self.solutions[media][target_reaction] = pkgmgr.getpkg("GapfillingPkg").run_test_conditions(self.test_conditions,self.solutions[media][target_reaction],self.test_condition_iteration_limit)
            if self.solutions[media][target_reaction] == None:
                print("Warning - no solution could be found that satisfied all specified test conditions in specified iterations!")
                return None
        if binary_check:
            return pkgmgr.getpkg("GapfillingPkg").binary_check_gapfilling_solution()
        return self.solutions[media][target_reaction] 
    
    def integrate_gapfill_solution(self,solution):
        for rxnid in solution["reversed"]:
            rxn = self.model.reactions.get_by_id(rxnid)
            if gfresults["reversed"][rxnid] == ">":
                rxn.upper_bound = 100
            else:
                rxn.lower_bound = -100
        for rxnid in solution["new"]:
            rxn = self.gfmodel.reactions.get_by_id(rxnid)
            rxn = rxn.copy()
            self.model.add_reactions([rxn])
            if solution["new"][rxnid] == ">":
                rxn.upper_bound = 100
                rxn.lower_bound = 0 
            else:
                rxn.upper_bound = 0
                rxn.lower_bound = -100
        return self.model
    
    @staticmethod
    def gapfill(model,media = None,target_reaction = "bio1",default_gapfill_templates = [],default_gapfill_models = [],test_conditions = [],reaction_scores = {},blacklist = []):
        gapfiller = MSGapfill(model,default_gapfill_templates,default_gapfill_models,test_conditions,reaction_scores,blacklist)
        gfresults = gapfiller.run_gapfilling(media,target_reaction)
        return gapfiller.integrate_gapfill_solution(gfresults)
    
    
    