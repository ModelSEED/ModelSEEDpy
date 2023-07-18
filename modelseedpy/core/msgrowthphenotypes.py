# -*- coding: utf-8 -*-
import pandas as pd
import logging
import cobra
from cobra.core.dictlist import DictList
from modelseedpy.core.msmedia import MSMedia
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgapfill import MSGapfill

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO

class MSGrowthPhenotype:
    def __init__(
        self,
        id,
        media=None,
        growth=None,
        gene_ko=[],
        additional_compounds=[],
        parent=None,
        name=None,
    ):
        self.id = id
        self.name = name
        if name == None:
            self.name = self.id
        self.growth = growth
        self.media = media
        self.gene_ko = gene_ko
        self.gapfilling = None
        self.additional_compounds = additional_compounds
        self.parent = parent

    def build_media(self,include_base_media=True):
        """Builds media object to use when simulating the phenotype
        Parameters
        ----------
        include_base_media : bool
            Indicates whether to include the base media for the phenotype set in the formulation
        """
        cpd_hash = {}
        for cpd in self.additional_compounds:
            cpd_hash[cpd] = 100
        full_media = MSMedia.from_dict(cpd_hash)
        if self.media:
            full_media.merge(self.media, overwrite_overlap=False)
        if include_base_media:
            if self.parent and self.parent.base_media:
                full_media.merge(self.parent.base_media, overwrite_overlap=False)
        return full_media

    def simulate(
        self,
        model_or_mdlutl,
        objective,
        growth_multiplier=3,
        add_missing_exchanges=False,
        save_fluxes=False,
        pfba=False,
    ):
        """Simulates a single phenotype
        Parameters
        ----------
        model_or_modelutl : Model | MSModelUtl
            Model to use to run the simulations
        add_missing_exchanges : bool
            Boolean indicating if exchanges for compounds mentioned explicitly in phenotype media should be added to the model automatically
        growth_multiplier : double
            Indicates a multiplier to use for positive growth above the growth on baseline media
        save_fluxes : bool
            Indicates if the fluxes should be saved and returned with the results
        pfba : bool
            Runs pFBA to compute fluxes after initially solving for growth
        """
        modelutl = model_or_mdlutl
        if not isinstance(model_or_mdlutl, MSModelUtil):
            modelutl = MSModelUtil.get(model_or_mdlutl)
        
        #Setting objective
        if objective:
            modelutl.model.objective = objective
        
        #Building full media and adding missing exchanges
        output = {"growth": None, "class": None, "missing_transports": [], "baseline_growth": None}
        full_media = self.build_media()
        if add_missing_exchanges:
            output["missing_transports"] = modelutl.add_missing_exchanges(full_media)
        
        #Getting basline growth
        output["baseline_growth"] = 0.01
        if self.parent:
            output["baseline_growth"] = self.parent.baseline_growth(modelutl,objective)
        
        #Building specific media and setting compound exception list
        if self.parent and self.parent.atom_limits and len(self.parent.atom_limits) > 0:
            reaction_exceptions = []
            specific_media = self.build_media(False)
            for mediacpd in specific_media.mediacompounds:
                ex_hash = mediacpd.get_mdl_exchange_hash(modelutl)
                for mdlcpd in ex_hash:
                    reaction_exceptions.append(ex_hash[mdlcpd])
            modelutl.pkgmgr.getpkg("ElementUptakePkg").build_package(self.parent.atom_limits,exception_reactions=reaction_exceptions)
        
        #Applying media
        if self.parent:
            modelutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(
                full_media, self.parent.base_uptake, self.parent.base_excretion
            )
        else:
            modelutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(
                full_media,0,1000
            )
        
        with modelutl.model:
            #Applying gene knockouts
            for gene in self.gene_ko:
                if gene in modelutl.model.genes:
                    geneobj = modelutl.model.genes.get_by_id(gene)
                    geneobj.knock_out()
            
            #Optimizing model
            solution = modelutl.model.optimize()
            output["growth"] = solution.objective_value
            if solution.objective_value > 0 and pfba:
                solution = cobra.flux_analysis.pfba(modelutl.model)
            if save_fluxes:
                output["fluxes"] = solution.fluxes
        
        #Determining phenotype class
        if output["growth"] >= output["baseline_growth"]*growth_multiplier:
            if not self.growth:
                output["class"] = "GROWTH"
            elif self.growth > 0:
                output["class"] = "CP"
            elif self.growth == 0:
                output["class"] = "FP"
        else:
            if not self.growth:
                output["class"] = "NOGROWTH"
            elif self.growth > 0:
                output["class"] = "FN"
            elif self.growth == 0:
                output["class"] = "CN"
        return output

    def gapfill_model_for_phenotype(
        self,
        msgapfill,
        objective,
        test_conditions,
        growth_multiplier=10,
        add_missing_exchanges=False,
    ):
        """Gapfills the model to permit this single phenotype to be positive
        Parameters
        ----------
        msgapfill : MSGapfill
            Fully configured gapfilling object
        add_missing_exchanges : bool
            Boolean indicating if exchanges for compounds mentioned explicitly in phenotype media should be added to the model automatically
        growth_multiplier : double
            Indicates a multiplier to use for positive growth above the growth on baseline media
        objective : string
            Expression for objective to be activated by gapfilling
        """
        #First simulate model without gapfilling to assess ungapfilled growth
        output = self.simulate(msgapfill.mdlutl,objective,growth_multiplier,add_missing_exchanges)
        if output["growth"] >= output["baseline_growth"]*growth_multiplier:
            #No gapfilling needed - original model grows without gapfilling
            return {"reversed": {}, "new": {},"media": self.build_media(), "target":objective, "minobjective": output["baseline_growth"]*growth_multiplier, "binary_check":False}
        
        #Now pulling the gapfilling configured model from MSGapfill
        gfmodelutl = MSModelUtil.get(msgapfill.gfmodel)
        #Saving the gapfill objective because this will be replaced when the simulation runs
        gfobj = gfmodelutl.model.objective
        #Running simulate on gapfill model to add missing exchanges and set proper media and uptake limit constraints
        output = self.simulate(modelutl,objective,growth_multiplier,add_missing_exchanges)
        #If the gapfilling model fails to achieve the minimum growth, then no solution exists
        if output["growth"] < output["baseline_growth"]*growth_multiplier:
            logger.warning(
                "Gapfilling failed with the specified model, media, and target reaction."
            )
            return None
        
        #Running the gapfilling itself
        full_media = self.build_media()
        with modelutl.model:
            #Applying gene knockouts
            for gene in self.gene_ko:
                if gene in modelutl.model.genes:
                    geneobj = modelutl.model.genes.get_by_id(gene)
                    geneobj.knock_out()
    
            gfresults = self.gapfilling.run_gapfilling(media,None,minimum_obj=output["baseline_growth"]*growth_multiplier)
            if gfresults is None:
                logger.warning(
                    "Gapfilling failed with the specified model, media, and target reaction."
                )
                    
        return gfresults

class MSGrowthPhenotypes:
    def __init__(self, base_media=None, base_uptake=0, base_excretion=1000,global_atom_limits={}):
        self.base_media = base_media
        self.phenotypes = DictList()
        self.base_uptake = base_uptake
        self.base_excretion = base_excretion
        self.atom_limits = global_atom_limits
        self.baseline_growth_data = {}
        self.cached_based_growth = {}

    @staticmethod
    def from_compound_hash(compounds,base_media=None, base_uptake=0, base_excretion=1000,global_atom_limits={}):
        growthpheno = MSGrowthPhenotypes(base_media, base_uptake, base_excretion,global_atom_limits)
        new_phenos = []
        for cpd in compounds:
            newpheno = MSGrowthPhenotype(cpd, None, compounds[cpd], [], [cpd])
            new_phenos.append(newpheno)
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    @staticmethod
    def from_kbase_object(data, kbase_api,base_media=None, base_uptake=0, base_excretion=1000,global_atom_limits={}):
        growthpheno = MSGrowthPhenotypes(base_media,base_uptake, base_excretion,global_atom_limits)
        new_phenos = []
        for pheno in data["phenotypes"]:
            media = kbase_api.get_from_ws(pheno["media_ref"], None)
            geneko = []
            for gene in pheno["geneko_refs"]:
                geneko.append(added_cpd.split("/").pop())
            added_compounds = []
            for added_cpd in pheno["additionalcompound_refs"]:
                added_compounds.append(added_cpd.split("/").pop())
            newpheno = MSGrowthPhenotype(
                pheno["id"], media, pheno["normalizedGrowth"], geneko, added_compounds
            )
            new_phenos.append(newpheno)
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    @staticmethod
    def from_kbase_file(filename, kbase_api,base_media=None, base_uptake=0, base_excretion=1000,global_atom_limits={}):
        # TSV file with the following headers:media    mediaws    growth    geneko    addtlCpd
        growthpheno = MSGrowthPhenotypes(base_media,base_uptake, base_excretion,global_atom_limits)
        headings = []
        new_phenos = []
        with open(filename) as f:
            lines = f.readlines()
            for line in lines:
                items = line.split("\t")
                if headings == None:
                    headings = items
                else:
                    data = {}
                    for i in range(0, len(items)):
                        data[headings[i]] = items[i]
                    data = FBAHelper.validate_dictionary(
                        headings,
                        ["media", "growth"],
                        {"mediaws": None, "geneko": [], "addtlCpd": []},
                    )
                    media = kbase_api.get_from_ws(data["media"], data["mediaws"])
                    id = data["media"]
                    if len(data["geneko"]) > 0:
                        id += "-" + ",".join(data["geneko"])
                    if len(data["addtlCpd"]) > 0:
                        id += "-" + ",".join(data["addtlCpd"])
                    newpheno = MSGrowthPhenotype(
                        id, media, data["growth"], data["geneko"], data["addtlCpd"]
                    )
                    new_phenos.append(newpheno)
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    @staticmethod
    def from_ms_file(filename,base_media=None, base_uptake=0, base_excretion=100,global_atom_limits={}):
        growthpheno = MSGrowthPhenotypes(base_media,base_uptake, base_excretion,global_atom_limits)
        df = pd.read_csv(filename)
        required_headers = ["Compounds", "Growth"]
        for item in required_headers:
            if item not in df:
                raise ValueError("Required header " + item + " is missing!")
        new_phenos = []
        for row in df.rows:
            cpds = row["Compounds"].split(";")
            id = row["Compounds"]
            if "ID" in row:
                id = row["ID"]
            geneko = []
            if "GeneKO" in row:
                geneko = row["GeneKO"].split(";")
            newpheno = MSGrowthPhenotype(id, None, row["Growth"], geneko, cpds)
            new_phenos.append(newpheno)
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    def build_super_media(self):
        super_media = None
        for pheno in self.phenotypes:
            if not super_media:
                super_media = pheno.build_media()
            else:
                super_media.merge(pheno.build_media(), overwrite_overlap=False)
        return super_media
    
    def add_phenotypes(self, new_phenotypes):
        keep_phenos = []
        for pheno in new_phenotypes:
            if pheno.id not in self.phenotypes:
                pheno.parent = self
                keep_phenos.append(pheno)
        additions = DictList(keep_phenos)
        self.phenotypes += additions

    def baseline_growth(
        self,
        model_or_mdlutl,
        objective
    ):
        """Simulates all the specified phenotype conditions and saves results
        Parameters
        ----------
        model_or_modelutl : Model | MSModelUtl
            Model to use to run the simulations
        """
        # Discerning input is model or mdlutl and setting internal links
        modelutl = model_or_mdlutl
        if not isinstance(model_or_mdlutl, MSModelUtil):
            modelutl = MSModelUtil.get(model_or_mdlutl)
        #Checking if base growth already computed
        if modelutl in self.cached_based_growth:
            if objective in self.cached_based_growth[modelutl]:
                return self.cached_based_growth[modelutl][objective]
        else:
            self.cached_based_growth[modelutl] = {}
        #Setting objective
        modelutl.objective = objective
        #Setting media
        modelutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(
            self.base_media, self.base_uptake, self.base_excretion
        )
        #Adding uptake limits
        if len(self.atom_limits) > 0:
            modelutl.pkgmgr.getpkg("ElementUptakePkg").build_package(self.atom_limits)
        #Simulating
        self.cached_based_growth[modelutl][objective] = modelutl.model.slim_optimize()
        return self.cached_based_growth[modelutl][objective]

    def simulate_phenotypes(
        self,
        model_or_modelutl,
        objective,
        add_missing_exchanges=False,
        correct_false_negatives=False,
        template=None,
        growth_threshold=0.01,
        save_fluxes=False
    ):
        """Simulates all the specified phenotype conditions and saves results
        Parameters
        ----------
        model_or_modelutl : Model | MSModelUtl
            Model to use to run the simulations
        objective : string
            Expression for objective to maximize in simulations
        add_missing_exchanges : bool
            Boolean indicating if exchanges for compounds mentioned explicitly in phenotype media should be added to the model automatically
        growth_multiplier : double
            Indicates a multiplier to use for positive growth above the growth on baseline media
        save_fluxes : bool
            Indicates if the fluxes should be saved and returned with the results
        """
        # Discerning input is model or mdlutl and setting internal links
        modelutl = model_or_mdlutl
        if not isinstance(model_or_mdlutl, MSModelUtil):
            modelutl = MSModelUtil.get(model_or_mdlutl)
        #Setting objective
        modelutl.objective = objective
        #Getting basline growth
        
        summary = {
            "Label": ["Accuracy", "CP", "CN", "FP", "FN","Growth","No growth"],
            "Count": [0, 0, 0, 0, 0,0,0],
        }
        data = {
            "Phenotype": [],
            "Observed growth": [],
            "Simulated growth": [],
            "Class": [],
            "Transports missing": [],
            "Gapfilled reactions": [],
        }
        for pheno in self.phenotypes:
            with model:
                result = pheno.simulate(
                    modelutl, growth_threshold, add_missing_exchanges, save_fluxes
                )  # Result should have "growth" and "class"
                if result["class"] == "FN" and correct_false_negatives:
                    pheno.gapfill_model_for_phenotype(modelutl, [template], None)
                    if pheno.gapfilling.last_solution != None:
                        list = []
                        for rxn_id in pheno.gapfilling.last_solution["reversed"]:
                            list.append(
                                pheno.gapfilling.last_solution["reversed"][rxn_id]
                                + rxn_id
                            )
                        for rxn_id in pheno.gapfilling.last_solution["new"]:
                            list.append(
                                pheno.gapfilling.last_solution["new"][rxn_id] + rxn_id
                            )
                        data["Gapfilled reactions"].append(";".join(list))
                    else:
                        data["Gapfilled reactions"].append(None)
                else:
                    data["Gapfilled reactions"].append(None)
                result = pheno.simulate(
                    modelutl, growth_threshold, add_missing_exchanges, save_fluxes
                )  # Result should have "growth" and "class"
                data["Class"].append(result["class"])
                data["Phenotype"].append(pheno.id)
                data["Observed growth"].append(pheno.growth)
                data["Simulated growth"].append(result["growth"])
                data["Transports missing"].append(
                    ";".join(result["missing_transports"])
                )
                if result["class"] == "CP":
                    summary["Count"][1] += 1
                    summary["Count"][0] += 1
                if result["class"] == "CN":
                    summary["Count"][2] += 1
                    summary["Count"][0] += 1
                if result["class"] == "FP":
                    summary["Count"][3] += 1
                if result["class"] == "FN":
                    summary["Count"][4] += 1

        summary["Count"][0] = summary["Count"][0] / len(self.phenotypes)
        sdf = pd.DataFrame(summary)
        df = pd.DataFrame(data)
        logger.info(df)
        return {"details": df, "summary": sdf}

    def fit_model_to_phenotypes(
        self,
        msgapfill,
        objective,
        grow_multiplier,
        correct_false_positives=False,
        minimize_new_false_positives=True,
        atp_safe=True,
        integrate_results=True,
        global_gapfilling=True
    ):
        
        """Simulates all the specified phenotype conditions and saves results
        Parameters
        ----------
        msgapfill : MSGapfill
            Gapfilling object used for the gapfilling process
        correct_false_positives : bool
            Indicates if false positives should be corrected
        minimize_new_false_positives : bool
            Indicates if new false positivies should be avoided
        integrate_results : bool
            Indicates if the resulting modifications to the model should be integrated
        """
        #Create super media for all 
        super_media = self.build_super_media()
        #Adding missing exchanges
        msgapfill.gfmodel.add_missing_exchanges(super_media)
        #Adding elemental constraints
        self.add_elemental_constraints()
        #Getting ATP tests
        
        #Filtering database for ATP tests
        
        #Penalizing database to avoid creating false positives
        
        #Building additional tests from current correct negatives
        
        #Computing base-line growth
        
        #Computing growth threshold
        
        #Running global gapfill
        
        #Integrating solution
        

    def gapfill_all_phenotypes(
        self,
        model_or_mdlutl,
        msgapfill=None,  # Needed if the gapfilling object in model utl is not initialized
        growth_threshold=None,
        add_missing_exchanges=False,
    ):
        mdlutl = MSModelUtil.get(model_or_mdlutl)
        # if msgapfill:
        #    mdlutl.gfutl = msgapfill
        # if not mdlutl.gfutl:
        #    logger.critical(
        #        "Must either provide a gapfilling object or provide a model utl with an existing gapfilling object"
        #    )
        # media_list = []
        # for pheno in self.phenotypes:
        #
        #
        # output = mdlutl.gfutl.run_multi_gapfill(
        #    media_list,
        #    default_minimum_objective=growth_threshold
        #    target=mdlutl.primary_biomass(),
        #
        #    binary_check=False,
        #    prefilter=True,
        #    check_for_growth=True,
        # )
