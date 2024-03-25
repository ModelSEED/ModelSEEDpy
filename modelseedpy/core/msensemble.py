# -*- coding: utf-8 -*-
import logging
import re
import time
import json
import sys
import pandas as pd
import cobra
import random
from cobra.core.dictlist import DictList
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msfba import MSFBA
from modelseedpy.core.msatpcorrection import MSATPCorrection

# from builtins import None

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO

class MSEnsemble:
    @staticmethod
    def from_models(models):
        #Converting models to MSModelUtil
        if not isinstance(model_or_mdlutl, MSModelUtil):
            for (i,mdl) in enumerate(models):
                models[i] = MSModelUtil.get(mdl)
        #Cloning the first model as a starting point
        clone_model = cobra.io.json.from_json(cobra.io.json.to_json(models[0].model))
        clone_mdlutl = MSModelUtil.get(clone_model)
        ensemble = MSEnsemble(clone_mdlutl)
        ensemble.rebuild_from_models(models)

    def from_annotation(model_or_mdlutl,reaction_probability_hash,sample_count=100):
        #Create genome from probabilities
        mdl = MSBuilder(genome,template).build(base_model, '0', False, False)
        mdl.template = self.gs_template
        mdlutl = MSModelUtil.get(mdl)
        ensemble = MSEnsemble(mdlutl)
        ensemble.build_ensemble(reaction_probability_hash, gpr_level_sampling, sample_count)
    
    def __init__(self,model_or_mdlutl):
        # Discerning input is model or mdlutl and setting internal links
        if isinstance(model_or_mdlutl, MSModelUtil):
            self.model = model_or_mdlutl.model
            self.mdlutl = model_or_mdlutl
        else:
            self.model = model_or_mdlutl
            self.mdlutl = MSModelUtil.get(model_or_mdlutl)
        self.data = {
            "size": self.size,
            "reactions": {}
        }
        attributes = self.mdlutl.get_attributes()
        if "ensemble" not in attributes:
            for rxn in self.model.reactions:
                self.data["reactions"][rxn.id] = {
                    "presence": "",
                    "gapfilling":"",
                    "genes": {}
                }
                for gene in rxn.genes:
                    self.data["reactions"][rxn.id]["genes"][gene.id] = {
                        "presence": ""
                    }
            logger.warning("Input model is not an ensemble model. You will need to run build_ensemble() to create an ensemble model.")
        else:
            self.data = attributes["ensemble"]

    def rebuild_from_models(self,models):#DONE
        #Clearing existing data
        self.data["ATP_analysis"] = {"core_atp_gapfilling":{},"selected_media":{},"tests":{}}
        for rxnid in self.data["reactions"]:
            self.data["reactions"][rxnid]["presence"] = ""
            self.data["reactions"][rxnid]["gapfilling"] = ""
            if "genes" in self.data["reactions"][rxnid]:
                for geneid in self.data["reactions"][rxnid]["genes"]:
                    self.data["reactions"][rxnid]["genes"][geneid]["presence"] = ""
            else:
                self.data["reactions"][rxnid]["genes"] = {}
        #Building presence strings from models
        self.data["size"] = len(models)
        for (i,mdlutl) in enumerate(models):
            attributes = mdlutl.get_attributes()
            if "ATP_analysis" in attributes:
                if "core_atp_gapfilling" in attributes["ATP_analysis"]:
                    for media in attributes["ATP_analysis"]["core_atp_gapfilling"]:
                        if media not in self.data["ATP_analysis"]["core_atp_gapfilling"]:
                            self.data["ATP_analysis"]["core_atp_gapfilling"][media] = []
                            for j in range(i):
                                self.data["ATP_analysis"]["core_atp_gapfilling"][media].append(None)
                        self.data["ATP_analysis"]["core_atp_gapfilling"][media].append(attributes["ATP_analysis"]["core_atp_gapfilling"][media])

                if "selected_media" in attributes["ATP_analysis"]:
                    for media in attributes["ATP_analysis"]["selected_media"]:
                        if media not in self.data["ATP_analysis"]["selected_media"]:
                            self.data["ATP_analysis"]["selected_media"][media] = []
                            for j in range(i):
                                self.data["ATP_analysis"]["selected_media"][media].append(None)
                        self.data["ATP_analysis"]["selected_media"][media].append(attributes["ATP_analysis"]["selected_media"][media])
                if "tests" in attributes["ATP_analysis"]:
                    for media in attributes["ATP_analysis"]["tests"]:
                        if media not in self.data["ATP_analysis"]["tests"]:
                            self.data["ATP_analysis"]["tests"][media] = {"objective":attributes["ATP_analysis"]["tests"][media]["objective"],"threshold":[]}
                            for j in range(i):
                                self.data["ATP_analysis"]["tests"][media]["threshold"].append(None)
                        self.data["ATP_analysis"]["tests"][media]["threshold"].append(attributes["ATP_analysis"]["tests"][media]["threshold"])
            add_reactions = []
            for rxn in mdlutl.model.reactions:
                if rxn.id not in self.mdlutl.model.reactions:
                    add_reactions.append(rxn)
                if rxn.id not in self.data["reactions"]:
                    self.data["reactions"][rxn.id] = {
                        "presence":'0' * i,
                        "genes":{}
                    }
                self.data["reactions"][rxn.id]["presence"] += "1"
                for gene in rxn.genes:
                    if gene.id not in self.data["reactions"][rxn.id]["genes"]:
                        self.data["reactions"][rxn.id]["genes"][gene.id] = '0' * i
                    self.data["reactions"][rxn.id]["genes"][gene.id] += "1"
            self.mdlutl.model.add_reactions(add_reactions)
        #Updating GPR of base model
        for rxnid in self.data["reactions"]:
            rxn = self.mdlutl.model.reactions.get_by_id(rxnid)
            rxn.gene_reaction_rule = " or ".join(self.data["reactions"][rxnid]["genes"].keys())
        #Computing probabilities from presence if missing
        for rxnid in self.ensemble_data["reactions"]:
            if "probabilty" not in self.ensemble_data["reactions"][rxnid]:
                self.ensemble_data["reactions"][rxnid]["probabilty"] = self.ensemble_data["reactions"][rxnid]["presence"].count('1')/len(self.ensemble_data["reactions"][rxnid]["presence"])
            for geneid in self.ensemble_data["reactions"][rxnid]["genes"]:
                if "probabilty" not in self.ensemble_data["reactions"][rxnid]["genes"][geneid]:
                    self.ensemble_data["reactions"][rxnid]["genes"][geneid]["probabilty"] = self.ensemble_data["reactions"][rxnid]["genes"][geneid]["presence"].count('1')/len(self.ensemble_data["reactions"][rxnid]["genes"][geneid]["presence"])

    def sample_from_probabilities(self,from_reaction_probabilities=False,sample_count=1000):
        #Scrolling through ensemble data with probabilities
        for rxnid in self.data["reactions"]:
            if "probability" not in self.data["reactions"][rxnid]:
                logger.critical("Reaction probability missing for "+rxnid+"!")
                return None
            if rxnid not in self.mdlutl.model.reactions:
                logger.critical("Reaction probability for "+rxnid+" but reaction not in base model!")
                return None
            rxn = self.mdlutl.model.reactions.get_by_id(rxnid)
            #Clearing existing data
            self.data["reactions"][rxnid]["presence"] = ""
            self.data["reactions"][rxnid]["gapfilling"] = ""
            #Loading gene-level data
            if "genes" not in self.data["reactions"][rxnid]:
                self.data["reactions"][rxnid]["genes"] = {}
            for gene in rxn.genes:
                if gene.id not in self.data["reactions"][rxnid]["genes"] or "probability" not in self.data["reactions"][rxnid]["genes"][gene.id]:
                    logger.warning("Reaction "+rxnid+" has gene "+gene.id+" but no associated probability data!")
                    self.data["reactions"][rxnid]["genes"][gene.id] = {"presence":"","probablity":1}
                self.data["reactions"][rxnid]["genes"][gene.id]["presence"] = ""
        #Sampling from probabilities
        for i in range(sample_count):
            if from_reaction_probabilities:
                for rxnid in self.data["reactions"]:
                    if random.uniform(0,1) < self.data["reactions"][rxnid]["probability"]:
                        self.data["reactions"][rxnid]["presence"] += "1"
                    else:
                        self.data["reactions"][rxnid]["presence"] += "0"
                for geneid in self.data["reactions"][rxnid]["genes"]:
                    self.data["reactions"][rxnid]["genes"][geneid]["presence"] += "1"
            else:
                present = False
                for geneid in self.data["reactions"][rxnid]["genes"]:
                    if random.uniform(0,1) < self.data["reactions"][rxnid]["genes"][geneid]["probability"]:
                        present = True
                        self.data["reactions"][rxnid]["genes"][geneid]["presence"] += "1"
                    else:
                        self.data["reactions"][rxnid]["genes"][geneid]["presence"] += "0"            
                if present:
                    self.data["reactions"][rxnid]["presence"] += "1"
                else:
                    self.data["reactions"][rxnid]["presence"] += "0"

    def unpack_models(self,model_list=None):
        output_models = [None]*self.size
        for i in range(self.size):
            if not model_list or i in model_list:
                clone_mdl = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
                clone_mdl_utl = MSModelUtil.get(clone_mdl)
                remove_reactions = []
                for rxn in clone_mdl_utl.model.reactions:
                    if rxn.id in self.data["reactions"]:
                        if self.data["reactions"][rxn.id][i] == "0":
                            remove_reactions.append(rxn)
                        else:
                            new_genes = []
                            for gene in rxn.genes:
                                if gene.id in self.data["reactions"][rxn.id]["genes"]:
                                    if self.data["reactions"][rxn.id]["genes"][gene.id][i] == "1":
                                        new_genes.append(gene)
                            rxn.gene_reaction_rule = " or ".join([gene.id for gene in new_genes])
                    else:
                        logger.warning("Ensemble model contains reaction not included in ensemble data. Removing reaction "+rxn.id+" from ensemble model.")
                        remove_reactions.append(rxn)
                clone_mdl.remove_reactions(remove_reactions)
                if "ATP_analysis" in self.data:
                    attributes = clone_mdl_utl.get_attributes()
                    attributes["ATP_analysis"] = {"core_atp_gapfilling":{},"selected_media":{},"tests":{}}
                    for media in self.data["ATP_analysis"]["core_atp_gapfilling"]:
                        if self.data["ATP_analysis"]["core_atp_gapfilling"][media][i] != None:
                            attributes["ATP_analysis"]["core_atp_gapfilling"][media] = self.data["ATP_analysis"]["core_atp_gapfilling"][media][i]
                    for media in self.data["ATP_analysis"]["selected_media"]:
                        if self.data["ATP_analysis"]["selected_media"][media][i] != None:
                            attributes["ATP_analysis"]["selected_media"][media] = self.data["ATP_analysis"]["selected_media"][media][i]
                    for media in self.data["ATP_analysis"]["tests"]:
                        if self.data["ATP_analysis"]["tests"][media]["threshold"][i] != None:
                            attributes["ATP_analysis"]["tests"][media] = {
                                "objective":self.data["ATP_analysis"]["tests"][media]["objective"],
                                "threshold":self.data["ATP_analysis"]["tests"]["threshold"][media][i]
                            }
                    clone_mdl_utl.save_attributes(attributes)
                output_models[i] = clone_mdl_utl
        return output_models

    def save_ensemble_model(self):
        self.mdlutl.save_attributes(self.data,"ensemble")
        return self.mdlutl

    def run_fba(self,media,objective,maximize,gene_ko=[],reaction_ko=[],pfba=True,fva=True):
        msfba = MSFBA(self.model,media,objective,maximize,gene_ko,reaction_ko,pfba,fva,clone=True)
        msfba.run()
        models = self.unpack_models()
        #Iterating over each model to run FBA on each
        for mdlutl in models:
            subfba = MSFBA(mdlutl,media,objective,maximize,gene_ko,reaction_ko,pfba,fva,clone=False)
            subfba.run()
            msfba.add_secondary_solution(subfba.primary_solution,subfba.fva_results)
        return msfba

    def run_atp_method(
        self,
        core_template=None,
        atp_medias=[],
        compartment="c0",
        max_gapfilling=10,
        gapfilling_delta=0,
        atp_hydrolysis_id=None,
        load_default_medias=True,
        forced_media=[],
        default_media_path=None,
    ):
        models = self.unpack_models()
        for mdlutl in models:
            atpcorrection = MSATPCorrection(
                core_template,
                atp_medias,
                compartment,
                max_gapfilling,
                gapfilling_delta,
                atp_hydrolysis_id,
                load_default_medias,
                forced_media,
                default_media_path
            )
            tests = atpcorrection.run_atp_correction()
        self.rebuild_from_models(models)

    def run_gapfilling(self):
        pass