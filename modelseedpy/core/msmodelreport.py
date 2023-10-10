# -*- coding: utf-8 -*-
import pandas as pd
import logging
import os
import re
import jinja2
from os.path import dirname
from pandas.io.formats.style import Styler
from modelseedpy.core.msmodelutl import MSModelUtil

module_path = dirname(os.path.abspath(__file__))

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO


class MSModelReport:
    def __init__(self, model_or_mdlutl):
        if isinstance(model_or_mdlutl, MSModelUtil):
            self.model = model_or_mdlutl.model
            self.modelutl = model_or_mdlutl
        else:
            self.model = model_or_mdlutl
            self.modelutl = MSModelUtil.get(model_or_mdlutl)

    def generate_reports(self, report_path, multi_tab_report_path):
        self.build_report(report_path)
        self.build_multitab_report(multi_tab_report_path)

    # Helper function to build overview data
    def build_overview_data(self):
        # Get the number of compartments
        number_compartments = len(
            set([metabolite.compartment for metabolite in self.model.metabolites])
        )

        # Extract gapfilling information
        core_gapfilling_media = []
        gapfilling_media = []
        gf_sensitivity = self.modelutl.attributes.get("gf_sensitivity", None)
        if gf_sensitivity:
            for media in gf_sensitivity:
                if (
                    "bio1" in self.modelutl.attributes["gf_sensitivity"][media]
                    and "success"
                    in self.modelutl.attributes["gf_sensitivity"][media]["bio1"]
                ):
                    gapfilling_media.append(media)
                if (
                    "rxn00062_c0" in self.modelutl.attributes["gf_sensitivity"][media]
                    and "success"
                    in self.modelutl.attributes["gf_sensitivity"][media]["rxn00062_c0"]
                ):
                    core_gapfilling_media.append(media)

        # Count the number of gapfills
        number_gapfills = len(gapfilling_media)

        # Convert the lists to strings
        core_gapfilling_str = (
            "; ".join(core_gapfilling_media)
            if core_gapfilling_media
            else "No core gapfilling needed."
        )
        gapfilling_media_str = (
            "; ".join(gapfilling_media)
            if gapfilling_media
            else "No genome-scale gapfilling."
        )

        overview = {
            "Model ID": self.model.id,
            "Full Gapfilling and ATP Analysis Report": "TBD",  # You may replace 'TBD' with actual data when available
            "Genome Scale Template": self.model.notes.get(
                "kbase_template_refs", "Data Not Available"
            ),
            "Core Gapfilling Media": core_gapfilling_str,
            "Gapfilling Media": gapfilling_media_str,
            "Source Genome": self.model.notes.get(
                "kbase_genome_ref", "Data Not Available"
            ),
            "Total Number of reactions": self.modelutl.nonexchange_reaction_count(),
            "Number compounds": len(self.model.metabolites),
            "Number compartments": number_compartments,
            "Number biomass": len(
                [
                    rxn
                    for rxn in self.model.reactions
                    if rxn.annotation.get("sbo") == "SBO:0000629"
                ]
            ),
            "Number gapfills": number_gapfills,
        }
        return overview

    # Helper function for extracting gapfilling data
    def extract_gapfilling_data(self, gf_sensitivity):
        if gf_sensitivity is None:
            return [], {}

        gapfilling_dict = {}
        gapfilling_summary = {}

        for media, media_data in gf_sensitivity.items():
            for target, target_data in media_data.items():
                gf_data = target_data.get("success", {})
                if isinstance(gf_data, dict):
                    for reaction_id, reaction_data in gf_data.items():
                        for direction, metabolites in reaction_data.items():
                            # If metabolites is None, set to empty string
                            if metabolites is None:
                                metabolites = ""

                            # Extract both IDs and Names for Gapfilling Sensitivity
                            sensitivity_ids = []
                            sensitivity_names = []
                            if isinstance(metabolites, (list, tuple)):
                                for met_id in metabolites:
                                    sensitivity_ids.append(met_id)
                                    met_name = (
                                        self.model.metabolites.get_by_id(met_id).name
                                        if met_id in self.model.metabolites
                                        else met_id
                                    )
                                    sensitivity_names.append(met_name)
                            else:
                                metabolites = str(metabolites)
                            entry = {
                                "reaction_id": reaction_id,
                                "reaction_name": self.model.reactions.get_by_id(
                                    reaction_id
                                ).name
                                if reaction_id in self.model.reactions
                                else reaction_id,
                                "media": media,
                                "direction": direction,
                                "target": target,
                                "gapfilling_sensitivity_id": "; ".join(sensitivity_ids)
                                if sensitivity_ids
                                else metabolites,
                                "gapfilling_sensitivity_name": "; ".join(
                                    sensitivity_names
                                )
                                if sensitivity_names
                                else metabolites,
                            }

                            # Update the summary dictionary
                            if reaction_id not in gapfilling_summary:
                                gapfilling_summary[reaction_id] = []
                            gapfilling_summary[reaction_id].append(
                                f"{media}: {direction}"
                            )

                            # Check if reaction_id is already in dictionary
                            if reaction_id in gapfilling_dict:
                                # Update the media
                                existing_entry = gapfilling_dict[reaction_id]
                                existing_media = existing_entry["media"].split("; ")
                                if media not in existing_media:
                                    existing_media.append(media)
                                    existing_entry["media"] = "; ".join(existing_media)
                            else:
                                gapfilling_dict[reaction_id] = entry

        return list(gapfilling_dict.values()), gapfilling_summary

    # transform data to be used in tabular format to use in build_model_report
    def transform_gapfilling_data(self, gapfilling_data):
        transformed_data = []
        for entry in gapfilling_data:
            row = [
                entry["reaction_id"],
                entry["reaction_name"],
                entry["media"],
                entry["direction"],
                entry["target"],
                entry["gapfilling_sensitivity_id"],
                entry["gapfilling_sensitivity_name"],
            ]
            transformed_data.append(row)
        return transformed_data

    # Extract ATP analysis data
    def extract_atp_analysis_data(self, atp_analysis, atp_expansion_filter):
        entries = []
        if atp_analysis and "core_atp_gapfilling" in atp_analysis:
            for media, data in atp_analysis["core_atp_gapfilling"].items():
                score = data.get("score", None)
                new_reactions = [
                    "{}: {}".format(k, v) for k, v in data.get("new", {}).items()
                ]
                reversed_reactions = [
                    "{}: {}".format(k, v) for k, v in data.get("reversed", {}).items()
                ]
                atp_production = "Not integrated"
                if (
                    "selected_media" in atp_analysis
                    and media in atp_analysis["selected_media"]
                ):
                    atp_production = atp_analysis["selected_media"][media]

                # Extracting the "Filtered Reactions" in the required format
                filtered_reactions = []
                for k, v in atp_expansion_filter.get(media, {}).items():
                    if isinstance(v, dict):
                        for sub_k, sub_v in v.items():
                            if isinstance(sub_v, dict):
                                for reaction, direction_dict in sub_v.items():
                                    direction = list(direction_dict.keys())[0]
                                    filtered_reactions.append(
                                        f"{reaction}: {direction}"
                                    )
                filtered_reactions_str = "; ".join(filtered_reactions)

                if score is not None:
                    entries.append(
                        {
                            "media": media,
                            "no_of_gapfilled_reactions": score,
                            "atp_production": atp_production,
                            "gapfilled_reactions": "; ".join(new_reactions),
                            "reversed_reaction_by_gapfilling": "; ".join(
                                reversed_reactions
                            ),
                            "filtered_reactions": filtered_reactions_str,
                        }
                    )
        # Sorting the entries based on the 'no_of_gapfilled_reactions' column
        entries.sort(key=lambda x: x["no_of_gapfilled_reactions"])
        return entries

    # Extract ATP production data for the ATP Analysis tab
    def extract_atp_production_data(self, atp_analysis):
        atp_production_dict = {}
        if atp_analysis:
            selected_media = atp_analysis.get("selected_media", {})
            core_atp_gapfilling = atp_analysis.get("core_atp_gapfilling", {})

            # First, process selected_media
            for media, value in selected_media.items():
                atp_production_dict[media] = round(value, 2)

            # Next, process core_atp_gapfilling for media not in selected_media
            for media, data in core_atp_gapfilling.items():
                if media not in atp_production_dict:
                    if data.get("failed"):
                        atp_production_dict[media] = "failed"
                    else:
                        # If the media was not processed in selected_media and it's not failed, set as 'Not Integrated'
                        atp_production_dict[media] = "Not Integrated"

        return atp_production_dict

    def build_multitab_report(self, output_path):

        # Build overview data
        overview_data = self.build_overview_data()

        # Get gf_sensitivity attribute from the model
        gf_sensitivity = self.modelutl.attributes.get("gf_sensitivity", None)

        # Extract gapfilling data
        gapfilling_entries, gapfilling_reaction_summary = self.extract_gapfilling_data(
            gf_sensitivity
        )

        # Check if ATP_analysis attribute is present in the model
        atp_analysis = self.modelutl.attributes.get("ATP_analysis", None)
        if atp_analysis:
            atp_expansion_filter = self.modelutl.attributes.get(
                "atp_expansion_filter", {}
            )
            atp_analysis_entries = self.extract_atp_analysis_data(
                atp_analysis, atp_expansion_filter
            )
        else:
            atp_analysis_entries = []

        # Initialize context dictionary
        context = {
            "overview": overview_data,
            "reactions": [],
            "compounds": [],
            "genes": [],
            "biomass": [],
            "gapfilling": gapfilling_entries,  # Populated with gapfilling data
            "atpanalysis": atp_analysis_entries,  # Populated with ATP analysis data
        }

        print("Module Path:", module_path + "/../data/")

        exchanges = {r.id for r in self.model.exchanges}

        # Identify biomass reactions using SBO annotation
        biomass_reactions_ids = {
            rxn.id
            for rxn in self.model.reactions
            if rxn.annotation.get("sbo") == "SBO:0000629"
        }

        # Reactions Tab
        for rxn in self.model.reactions:
            if rxn.id not in exchanges and rxn.id not in biomass_reactions_ids:
                equation = rxn.build_reaction_string(use_metabolite_names=True)
                rxn_data = {
                    "id": rxn.id,
                    "name": rxn.name,
                    "equation": equation,
                    "genes": rxn.gene_reaction_rule,
                    "gapfilling": "; ".join(
                        gapfilling_reaction_summary.get(rxn.id, [])
                    ),  # Empty list results in an empty string
                }
                context["reactions"].append(rxn_data)

        # Compounds Tab
        for cpd in self.model.metabolites:
            cpd_data = {
                "id": cpd.id,
                "name": cpd.name,
                "formula": cpd.formula,
                "charge": cpd.charge,
                "compartment": cpd.compartment,
            }
            context["compounds"].append(cpd_data)

        # Genes Tab
        for gene in self.model.genes:
            gene_data = {
                "gene": gene.id,
                "reactions": "; ".join([rxn.id for rxn in gene.reactions]),
            }
            context["genes"].append(gene_data)

        # Biomass Tab
        if biomass_reactions_ids:
            for biomass_rxn_id in biomass_reactions_ids:
                biomass_rxn = self.model.reactions.get_by_id(biomass_rxn_id)
                for metabolite, coefficient in biomass_rxn.metabolites.items():
                    compound_id = metabolite.id
                    compound_name = metabolite.name.split("_")[0]
                    compartment = compound_id.split("_")[-1]

                    biomass_data = {
                        "biomass_reaction_id": biomass_rxn.id,
                        "biomass_compound_id": compound_id,
                        "name": compound_name,
                        "coefficient": coefficient,
                        "compartment": compartment,
                    }
                    context["biomass"].append(biomass_data)
        else:
            print("No biomass reactions found in the model.")

        # Gapfilling Tab
        gf_sensitivity = self.modelutl.attributes.get("gf_sensitivity", None)
        gapfilling_data = self.extract_gapfilling_data(gf_sensitivity)
        context["gapfilling"] = gapfilling_entries

        # Extract ATP Production Data
        atp_production_data = self.extract_atp_production_data(atp_analysis)

        # Populate the 'atpanalysis' context with ATP production data
        for entry in context["atpanalysis"]:
            media = entry["media"]
            entry["atp_production"] = atp_production_data.get(media, None)

        # Diagnostics
        unique_biomass_rxns = biomass_reactions_ids
        print(f"Unique biomass reactions identified: {len(unique_biomass_rxns)}")
        print(f"Biomass Reaction IDs: {', '.join(unique_biomass_rxns)}")

        print("\nFirst 2 reactions:")
        for rxn in context["reactions"][:2]:
            print(rxn)

        print("\nFirst 2 compounds:")
        for cpd in context["compounds"][:2]:
            print(cpd)

        print("\nFirst 2 genes:")
        for gene in context["genes"][:2]:
            print(gene)

        print("\nFirst 2 biomass compounds:")
        for bm in context["biomass"][:2]:
            print(bm)

        print("\nFirst 2 gapfilling entries:")
        for gf in context["gapfilling"][:2]:
            print(gf)

        print("\nFirst 2 ATP Analysis entries:")
        for entry in context["atpanalysis"][:2]:
            print(entry)

        # Render with template
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(module_path + "/../data/"),
            autoescape=jinja2.select_autoescape(["html", "xml"]),
        )
        html = env.get_template("ModelReportTemplate.html").render(context)
        directory = dirname(output_path)
        os.makedirs(directory, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(html)

    def build_report(self, output_path):
        """Builds model HTML report for the Model Summary table
        Parameters
        ----------
        model : cobra.Model
            Model to use to build the report
        """

        # 1. Utilize the build_overview_data method
        model_summary_data = self.build_overview_data()
        # Remove the unwanted entry
        model_summary_data.pop("Full Gapfilling and ATP Analysis Report", None)
        # 2. Transform the dictionary into a list of tuples
        model_summary_list = [(key, value) for key, value in model_summary_data.items()]
        # 3. Convert to DataFrame
        model_summary_df = pd.DataFrame(model_summary_list, columns=["", ""])

        # Style the DataFrame (as was done previously)
        model_summary_df_styled = model_summary_df.style.hide(
            axis="index"
        ).set_table_styles(
            [
                {
                    "selector": "th",
                    "props": [
                        ("border", "none"),
                        ("background-color", "white"),
                        ("font-family", "Oxygen"),
                        ("font-size", "14px"),
                        ("line-height", "20px"),
                    ],
                },
                {
                    "selector": "td",
                    "props": [
                        ("border", "none"),
                        ("font-family", "Oxygen"),
                        ("font-size", "14px"),
                        ("line-height", "20px"),
                    ],
                },
                {
                    "selector": "tr:nth-child(even)",
                    "props": [("background-color", "white")],
                },
                {
                    "selector": "tr:nth-child(odd)",
                    "props": [("background-color", "#f2f2f2")],
                },
            ]
        )

        # Fetching the gapfilling sensitivity data
        gf_sensitivity = self.modelutl.attributes.get("gf_sensitivity", None)
        gapfilling_data = self.extract_gapfilling_data(gf_sensitivity)
        gapfilling_list = self.transform_gapfilling_data(gapfilling_data[0])

        # Convert the gapfilling_list to a DataFrame
        gapfillings_analysis_df = pd.DataFrame(
            gapfilling_list,
            columns=[
                "Reaction ID",
                "Reaction Name",
                "Media",
                "Direction",
                "Target",
                "Gapfilling Sensitivity ID",
                "Gapfilling Sensitivity Name",
            ],
        )

        # Apply style to Gapfillings Analysis DataFrame
        gapfillings_analysis_df_styled = gapfillings_analysis_df.style.hide(
            axis="index"
        ).set_table_styles(
            [
                {
                    "selector": "th",
                    "props": [
                        ("border", "none"),
                        ("background-color", "white"),
                        ("font-family", "Oxygen"),
                        ("font-size", "14px"),
                        ("line-height", "20px"),
                    ],
                },
                {
                    "selector": "td",
                    "props": [
                        ("border", "none"),
                        ("font-family", "Oxygen"),
                        ("font-size", "14px"),
                        ("line-height", "20px"),
                    ],
                },
                {
                    "selector": "tr:nth-child(even)",
                    "props": [("background-color", "white")],
                },
                {
                    "selector": "tr:nth-child(odd)",
                    "props": [("background-color", "#f2f2f2")],
                },
            ]
        )

        # Legend for Gapfillings Analysis
        annotations_text_gapfillings = """
            <ul>
                <li><b>Reaction ID:</b> The identifier of the reaction.</li>
                <li><b>Reaction Name:</b> The name of the reaction.</li>
                <li><b>Media:</b> The media used by gap filling.</li>
                <li><b>Direction:</b> The direction of the reaction. Can be ">" for forward, "<" for reverse, or "=" for both directions.</li>
                <li><b>Target:</b> The reaction selected as the objective function target for the gapfilling optimization problem. Targets here can be the model’s biomass reaction, commonly named “bio1” for models created by this app.
                Alternatively, “rxn00062” (ATP Production) reaction is shown for cases where gapfilling was applied to guarantee ATP production in a given media.
                When reactions are gapfilled for ATP production, we recommend checking the full Core ATP Analysis in the table below. </li>
                <li><b>Gapfilling Sensitivity ID and Name:</b> Gapfilling is necessary when compounds in the biomass objective function can not be produced by the model.
                For each reaction we list the biomass compound(s) that can not be synthesized by the model without gapfilling.
                In cases where gap filling fails there are two possible scenarios:
                1) FBF (failed before filtering) : the gapfilling immediately failed, even before we filtered out the ATP breaking reactions. This means this objective CANNOT be satisfied with the entire current database.
                2) FAF (failed after filtering): the gapfilling succeeded before filtering, but failed after filtering out reactions that break ATP. This tells you definitively if the ATP filtering caused the gapfilling to fail</li>
            </ul>
        """

        # Extract ATP analysis data
        atp_analysis = self.modelutl.attributes.get("ATP_analysis", None)
        atp_expansion_filter = self.modelutl.attributes.get("atp_expansion_filter", {})
        atp_analysis_entries = self.extract_atp_analysis_data(
            atp_analysis, atp_expansion_filter
        )

        # Convert the atp_analysis_entries list to a DataFrame
        atp_analysis_df = pd.DataFrame(atp_analysis_entries)

        # Apply style to ATP Analysis DataFrame
        atp_analysis_df_styled = atp_analysis_df.style.hide(
            axis="index"
        ).set_table_styles(
            [
                {
                    "selector": "th",
                    "props": [
                        ("border", "none"),
                        ("background-color", "white"),
                        ("font-family", "Oxygen"),
                        ("font-size", "14px"),
                        ("line-height", "20px"),
                    ],
                },
                {
                    "selector": "td",
                    "props": [
                        ("border", "none"),
                        ("font-family", "Oxygen"),
                        ("font-size", "14px"),
                        ("line-height", "20px"),
                    ],
                },
                {
                    "selector": "tr:nth-child(even)",
                    "props": [("background-color", "white")],
                },
                {
                    "selector": "tr:nth-child(odd)",
                    "props": [("background-color", "#f2f2f2")],
                },
            ]
        )

        # Legend for ATP Analysis
        annotations_text_atp_analysis = """
            <ul>
                <li><b>No. of gapfilled reactions:</b> The number of reactions filled by the gapfilling process.</li>
                <li><b>Media:</b> The media in which the reaction takes place.</li>
                <li><b>ATP Production:</b> ATP production by the core metabolism model.</li>
                <li><b>Gapfilled Reactions:</b> Reactions added during the gapfilling process.</li>
                <li><b>Reversed Reaction by Gapfilling:</b> Reactions that have been reversed during the gapfilling process.</li>
                <li><b>Filtered Reactions:</b> Reactions that have been filtered out during the analysis. When a reaction addition would lead to a large increase in ATP production or an infinite energy loop, we filter that reaction out of the gapfilling database and prevent it from being added to the model.</li>
            </ul>
        """

        # ATP analysis explanation text
        explanation_text_atp_analysis = """
        <p>During model reconstruction, we analyze the genome’s core metabolism draft model (model without gapfilling) to assess energy biosynthesis capabilities.
        The goal of this analysis is to ensure the core metabolism model is able to produce ATP before we expand the model to the genome-scale.
        This step is designed to prevent gapfilling from introducing reactions that create energy-generating loops.
        The tests are conducted on a large collection of minimal conditions, with the goal of simulating the model’s capability to produce energy with different electron donor, electron acceptor, and carbon source combinations.</p>
        <p><u>When the draft model of the core metabolism is capable of producing ATP in at least one of the test media, no gapfilling reactions part of this analysis will be added to the model.</u> While we still report the gapfilling requirements for the test media formulations that fail to produce ATP with that draft core model, we only integrate these solutions in the model when no test media succeeds in producing ATP.
        In this case, the integrated gap-filling solution(s) will be displayed in the “Gapfilling Analysis” table above, with the “Target” “rxn00062” (ATP Production) objective function.</p>
        <p>The goal is to display the test results for all media to provide clues for the metabolic capabilities of the genome(s). When many reactions are required for growth on the SO4 testing media conditions, this could be a good indicator that the organism is not capable of performing sulfate reduction.
        On the other hand, when only one gapfill reaction is required  for ATP production in a given media, multiple scenarios can be considered.
        1) Organism(s) can’t grow on test condition, and we correctly did not add the reaction to the model.  2) Possible issue with the source genome annotation missing a specific gene function 3) Possible issue with the model reconstruction database. We hope this data helps make more informed decisions on reactions that may need to be manually curated in the model.
        In cases where is known from the literature or unpublished experimental results that an organism is capable of producing ATP in a given media condition that requires gapfilling in this analysis, you can use the parameter “Force ATP media” in the reconstruction app to ensure those reactions are integrated into the model.
        .</p>
        """

        # Save the data to HTML with the styled DataFrames and the legends
        directory = os.path.dirname(output_path)
        os.makedirs(directory, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write('<meta charset="UTF-8">')
            f.write("<h1>Model Summary</h1>")
            f.write(model_summary_df_styled.render(escape=False))
            f.write("<br><br>")
            f.write("<h1>Gapfillings Analysis</h1>")

            # Check for Gapfillings Analysis data
            if not gapfillings_analysis_df.empty:
                f.write(gapfillings_analysis_df_styled.render(escape=False))
                f.write(f"<br><br><h3>Legend:</h3>{annotations_text_gapfillings}")
            else:
                f.write(
                    "<p><b>Warning:</b> No Gapfillings Analysis data available for this model.</p>"
                )

            f.write("<h1>Core ATP Analysis</h1>")

            # Check for ATP Analysis data
            if not atp_analysis_df.empty:
                f.write(atp_analysis_df_styled.render(escape=False))
                f.write(f"<br><br><h3>Legend:</h3>{annotations_text_atp_analysis}")
                f.write(explanation_text_atp_analysis)
            else:
                f.write(
                    "<p><b>Warning:</b> No Core ATP Analysis data available for this model.</p>"
                )
