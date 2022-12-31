# -*- coding: utf-8 -*-
import pandas as pd
import logging
import cobra
from cobra.core.dictlist import DictList
from modelseedpy.core.msmedia import MSMedia
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper

logger = logging.getLogger(__name__)


class MSGrowthPhenotype:
    def __init__(
        self,
        obj_id,
        media=None,
        growth=None,
        gene_ko=[],
        additional_compounds=[],
        parent=None,
        name=None,
    ):
        self.id = obj_id
        self.name = name or obj_id
        self.growth = growth
        self.media = media
        self.gene_ko = gene_ko
        self.gapfilling = None
        self.additional_compounds = additional_compounds
        self.parent = parent

    def build_media(self):
        cpd_hash = {}
        for cpd in self.additional_compounds:
            cpd_hash[cpd] = 100
        full_media = MSMedia.from_dict(cpd_hash)
        if self.media is not None:
            full_media.merge(self.media, overwrite_overlap=False)
        if self.parent is not None and self.parent.base_media is not None:
                full_media.merge(self.parent.base_media, overwrite_overlap=False)
        return full_media

    def simulate(
        self,
        modelutl,
        growth_threshold=0.001,
        add_missing_exchanges=False,
        save_fluxes=False,
        pfba=False,
    ):
        if not isinstance(modelutl, MSModelUtil):
            modelutl = MSModelUtil(modelutl)
        media = self.build_media()
        output = {"growth": None, "class": None, "missing_transports": []}
        if add_missing_exchanges:
            output["missing_transports"] = modelutl.add_missing_exchanges(media)
        pkgmgr = MSPackageManager.get_pkg_mgr(modelutl.model)
        pkgmgr.getpkg("KBaseMediaPkg").build_package(
            media, self.parent.base_uptake, self.parent.base_excretion
        )
        for gene in self.gene_ko:
            if gene in modelutl.model.genes:
                geneobj = modelutl.model.genes.get_by_id(gene)
                geneobj.knock_out()
        solution = modelutl.model.optimize()
        output["growth"] = solution.objective_value
        if solution.objective_value > 0 and pfba:
            solution = cobra.flux_analysis.pfba(modelutl.model)
        if save_fluxes:
            output["fluxes"] = solution.fluxes
        if output["growth"] >= growth_threshold:
            if self.growth > 0:
                results["class"] = "CP"
        else:    
            results["class"] = "CN"
            if self.growth > 0:
                results["class"] = "FN"
        return results  

    def gapfill_model_for_phenotype(
        self,
        modelutl,
        default_gapfill_templates,
        test_conditions,
        default_gapfill_models=[],
        blacklist=[],
        growth_threshold=0.001,
        add_missing_exchanges=False,
    ):
        if not isinstance(modelutl, MSModelUtil):
            modelutl = MSModelUtil(modelutl)
        self.gapfilling = MSGapfill(
            modelutl.model,
            default_gapfill_templates,
            default_gapfill_models,
            test_conditions,
            modelutl.reaction_scores(),
            blacklist,
        )
        media = self.build_media()
        if add_missing_exchanges:
            modelutl.add_missing_exchanges(media)
        for gene in self.gene_ko:
            if gene in modelutl.model.genes:
                geneobj = modelutl.model.genes.get_by_id(gene)
                geneobj.knock_out()
        gfresults = self.gapfilling.run_gapfilling(media, None)
        if gfresults is None:
            logger.warning(
                "Gapfilling failed with the specified model, media, and target reaction."
            )
        return self.gapfilling.integrate_gapfill_solution(gfresults)


class MSGrowthPhenotypes:
    def __init__(self, base_media=None, base_uptake=0, base_excretion=1000):
        self.base_media = base_media
        self.phenotypes = DictList()
        self.base_uptake = base_uptake
        self.base_excretion = base_excretion

    @staticmethod
    def from_compound_hash(compounds, base_media, base_uptake=0, base_excretion=1000):
        growthpheno = MSGrowthPhenotypes(base_media, base_uptake, base_excretion)
        new_phenos = []
        for cpd in compounds:
            newpheno = MSGrowthPhenotype(cpd, None, compounds[cpd], [], [cpd])
            new_phenos.append(newpheno)
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    @staticmethod
    def from_kbase_object(data, kbase_api):
        growthpheno = MSGrowthPhenotypes(None, 0, 1000)
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
    def from_kbase_file(filename, kbase_api):
        # TSV file with the following headers:media    mediaws    growth    geneko    addtlCpd
        from pandas import read_table
        table_df = read_table(filename)
        
        growthpheno = MSGrowthPhenotypes(base_media, 0, 1000)
        new_phenos = []
        for index, row in table_df.iterrows():
            data = FBAHelper.validate_dictionary(
                row.to_dict(), ["media","growth"], {"mediaws":None, "geneko":[], "addtlCpd":[]})
            media = kbase_api.get_from_ws(data["media"], data["mediaws"])
            media_id = data["media"]
            if len(data["geneko"]) > 0:
                media_id += "-"+",".join(data["geneko"])
            if len(data["addtlCpd"]) > 0:
                media_id += "-"+",".join(data["addtlCpd"])
            new_phenos.append(MSGrowthPhenotype(
                media_id, media, data["growth"], data["geneko"], data["addtlCpd"]))
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    @staticmethod
    def from_ms_file(filename, base_media, base_uptake=0, base_excretion=100):
        from pandas import read_csv
        df = read_csv(filename)
        
        growthpheno = MSGrowthPhenotypes(base_media,base_uptake,base_excretion)
        required_headers = ["Compounds","Growth"]
        for item in required_headers:
            if item not in df:
                raise ValueError(f'The required header {item} is missing from the CSV file.')
        new_phenos = []
        for index, row in df.iterrows():
            cpds = row["Compounds"].split(";")
            cpd_id = row["Compounds"]
            if "ID" in row:
                cpd_id = row["ID"]
            geneko = []
            if "GeneKO" in row:
                geneko = row["GeneKO"].split(";")
            newpheno = MSGrowthPhenotype(id, None, row["Growth"], geneko, cpds)
            new_phenos.append(newpheno)
        growthpheno.add_phenotypes(new_phenos)
        return growthpheno

    def add_phenotypes(self, new_phenotypes):
        keep_phenos = []
        for pheno in new_phenotypes:
            if pheno.id not in self.phenotypes:
                pheno.parent = self
                keep_phenos.append(pheno)
        additions = DictList(keep_phenos)
        self.phenotypes += additions

    def simulate_phenotypes(
        self,
        model,
        biomass,
        add_missing_exchanges=False,
        correct_false_negatives=False,
        template=None,
        growth_threshold=0.001,
        save_fluxes=False,
    ):
        from pandas import DataFrame
        data_df = DataFrame(columns = ["Phenotype","Observed growth","Simulated growth","Class","Transports missing","Gapfilled reactions"])
        
        model.objective = biomass
        modelutl = MSModelUtil(model)
        summary = {
            "Label": ["Accuracy", "CP", "CN", "FP", "FN"],
            "Count": [0, 0, 0, 0, 0],
        }
        for index, pheno in enumerate(self.phenotypes):
            with model:
                result = pheno.simulate(modelutl,growth_threshold,add_missing_exchanges)
                gfl = None
                if result["class"] == "FN" and correct_false_negatives:
                    pheno.gapfill_model_for_phenotype(modelutl,[template],None)
                    if pheno.gapfilling.last_solution:
                        gpfl_rxns = []
                        for rxn_id in pheno.gapfilling.last_solution["reversed"]:
                            gpfl_rxns.append(pheno.gapfilling.last_solution["reversed"][rxn_id]+rxn_id)
                        for rxn_id in pheno.gapfilling.last_solution["new"]:
                            gpfl_rxns.append(pheno.gapfilling.last_solution["new"][rxn_id]+rxn_id)
                        gfl = ";".join(gpfl_rxns)
                result = pheno.simulate(modelutl,growth_threshold,add_missing_exchanges)
                data_df.loc[index] = {
                    "Phenotype":pheno.id, "Observed growth": pheno.growth,
                    "Simulated growth":result["growth"], "class":result["class"], 
                    "Transports missing": ";".join(result["missing_transports"]), 
                    "Gapfilled reactions":gfl}
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
                
        summary["Count"][0] = summary["Count"][0]/len(self.phenotypes) 
        sdf = DataFrame(summary)
        logger.info(data_df)
        return {"details":data_df,"summary":sdf}
