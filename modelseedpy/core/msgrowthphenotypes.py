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

    def build_media(self):
        cpd_hash = {}
        for cpd in self.additional_compounds:
            cpd_hash[cpd] = 100
        full_media = MSMedia.from_dict(cpd_hash)
        if self.media != None:
            full_media.merge(self.media, overwrite_overlap=False)
        if self.parent != None and self.parent.base_media != None:
            full_media.merge(parent.base_media, overwrite_overlap=False)
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
                output["class"] = "CP"
            else:
                output["class"] = "FP"
        else:
            if self.growth > 0:
                output["class"] = "FN"
            else:
                output["class"] = "CN"
        return output

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
        growthpheno = MSGrowthPhenotypes(base_media, 0, 1000)
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
    def from_ms_file(filename, basemedia, base_uptake=0, base_excretion=100):
        growthpheno = MSGrowthPhenotypes(base_media, base_uptake, base_excretion)
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
        model.objective = biomass
        modelutl = MSModelUtil(model)
        summary = {
            "Label": ["Accuracy", "CP", "CN", "FP", "FN"],
            "Count": [0, 0, 0, 0, 0],
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
