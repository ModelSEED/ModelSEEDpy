# -*- coding: utf-8 -*-
import pandas as pd
import logging
import matplotlib.cm as cm
import os
import jinja2
from os.path import dirname, exists
from modelseedpy.core.msmodelutl import MSModelUtil
module_path = dirname(os.path.abspath(__file__))

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO

class MSModelReport:
    def __init__(
        self
    ):
        pass

    def build_multitab_report(self, model_or_mdlutl, output_path):

        # Helper function for extracting gapfilling data
        def extract_gapfilling_data(gf_sensitivity, model):
            gapfilling_entries = []

            if not gf_sensitivity:
                return []

            for media, media_data in gf_sensitivity.items():
                for target, target_data in media_data.items():
                    for reaction_id, reaction_data in target_data.get('success', {}).items():
                        for direction, metabolites in reaction_data.items():
                            entry = {
                                "reaction_id": reaction_id,
                                "reaction_name": model.reactions.get_by_id(reaction_id).name if reaction_id in model.reactions else reaction_id,
                                "media": media,
                                "direction": direction,
                                "target": target,
                                "gapfilling_sensitivity": "; ".join(metabolites) if isinstance(metabolites, (list, tuple)) else str(metabolites)
                            }
                            gapfilling_entries.append(entry)

            return gapfilling_entries

        context = {
            "overview": [{"Model": "Model 3", "Genome": "Genome C"}, {"Model": "Model 4", "Genome": "Genome D"}],
            "reactions": [],
            "compounds": [],
            "genes": [],
            "biomass": [],
            "gapfilling": [],
            "atpanalysis": [{'no_of_gapfilled_reactions': 5, 'media': 'LB', 'atp_production': 'High', 'gapfilled_reactions': 'R005; R006', 'reversed_reaction_by_gapfilling': 'R007', 'filtered_reactions': 'R008; R009'}],
        }

        print("Module Path:", module_path + "/../data/")

        exchanges = {r.id for r in model_or_mdlutl.exchanges}

        # Identify biomass reactions using SBO annotation
        biomass_reactions_ids = {rxn.id for rxn in model_or_mdlutl.reactions if rxn.annotation.get('sbo') == 'SBO:0000629'}

        # Reactions Tab
        for rxn in model_or_mdlutl.reactions:
            if rxn.id not in exchanges and rxn.id not in biomass_reactions_ids:
                equation = rxn.build_reaction_string(use_metabolite_names=True)
                rxn_data = {
                    "id": rxn.id,
                    "name": rxn.name,
                    "equation": equation,
                    "genes": rxn.gene_reaction_rule,
                    "gapfilling": "TBD"
                }
                context["reactions"].append(rxn_data)

        # Compounds Tab
        for cpd in model_or_mdlutl.metabolites:
            cpd_data = {
                "id": cpd.id,
                "name": cpd.name,
                "formula": cpd.formula,
                "charge": cpd.charge,
                "compartment": cpd.compartment
            }
            context["compounds"].append(cpd_data)

        # Genes Tab
        for gene in model_or_mdlutl.genes:
            gene_data = {
                "gene": gene.id,
                "reactions": "; ".join([rxn.id for rxn in gene.reactions])
            }
            context["genes"].append(gene_data)

        # Biomass Tab
        if biomass_reactions_ids:
            for biomass_rxn_id in biomass_reactions_ids:
                biomass_rxn = model_or_mdlutl.reactions.get_by_id(biomass_rxn_id)
                for metabolite, coefficient in biomass_rxn.metabolites.items():
                    compound_id = metabolite.id
                    compound_name = metabolite.name.split('_')[0]
                    compartment = compound_id.split('_')[-1]

                    biomass_data = {
                        "biomass_reaction_id": biomass_rxn.id,
                        "biomass_compound_id": compound_id,
                        "name": compound_name,
                        "coefficient": coefficient,
                        "compartment": compartment
                    }
                    context["biomass"].append(biomass_data)
        else:
            print("No biomass reactions found in the model.")

        # Gapfilling Tab
        gf_sensitivity = model_or_mdlutl.attributes.get('gf_sensitivity', None)
        gapfilling_data = extract_gapfilling_data(gf_sensitivity, model_or_mdlutl)
        context["gapfilling"] = gapfilling_data

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

        # Render with template
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(module_path + "/../data/"),
            autoescape=jinja2.select_autoescape(['html', 'xml'])
        )
        html = env.get_template("ModelReportTemplate.html").render(context)
        directory = dirname(output_path)
        os.makedirs(directory, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(html)


       
    def build_report(
        self,
        model_or_mdlutl,
        output_path
    ):
        """Builds model HTML report
        Parameters
        ----------
        model_or_modelutl : Model | MSModelUtl
            Model to use to run the simulations
        """
        modelutl = model_or_mdlutl
        if not isinstance(model_or_mdlutl, MSModelUtil):
            modelutl = MSModelUtil.get(model_or_mdlutl)
        
        # Process the data
        attributes = modelutl.get_attributes()
        
        selected_media_data = attributes['ATP_analysis']['selected_media']
        core_atp_gapfilling_data = attributes['ATP_analysis']['core_atp_gapfilling']
        gf_filter_data = attributes['gf_filter']
        gf_sensitivity_data = attributes.get('gf_sensitivity')  # Get 'gf_sensitivity_data' if available, otherwise it will be None
        number_gapfills = 0
        if gf_sensitivity_data:
            number_gapfills = len(gf_sensitivity_data)
        
        # Get the names of 'Core Gapfilling Media' and 'Gapfilling Media'
        core_gapfilling_media = [media for media, media_data in (gf_sensitivity_data or {}).items() if 'rxn00062_c0' in media_data]
        gapfilling_media = [media for media, media_data in (gf_sensitivity_data or {}).items() if 'bio1' in media_data]
        core_gapfilling_media_text = ', '.join(core_gapfilling_media)
        gapfilling_media_text = ', '.join(gapfilling_media)
        
        bio_count = 0
        for rxn in modelutl.model.reactions:
            if rxn.id[0:3] == "bio":
                bio_count += 1
        
        # Create the Model Summary table data
        model_summary_data = [
            ('<b>Model ID</b>', modelutl.wsid),
            ('<b>Genome Scale Template</b>',  modelutl.model.template_ref),
            #('<b>Core Template</b>', modelutl.model.core_template_ref),
            ('<b>Core Gapfilling Media</b>', core_gapfilling_media_text),
            ('<b>Gapfilling Media</b>', gapfilling_media_text),
            ('<b>Source Genome</b>',modelutl.model.name),
            ('<b>Total Number of reactions</b>', str(len(modelutl.model.reactions))),
        #    ('<b>Number of reactions in Core</b>', 'TBD - attributes require changes things to support this'),
        #    ('<b>Number of reactions in Genome Scale</b>', 'TBD - attributes require changes things to support this'),
            ('<b>Number compounds</b>', str(len(modelutl.model.metabolites))),
            ('<b>Number compartments</b>', str(len(modelutl.model.compartments))),
            ('<b>Number biomass</b>', str(bio_count)),
            ('<b>Number gapfills</b>', str(number_gapfills)),
        ]
        
        # Create the DataFrame for Model Summary
        model_summary_df = pd.DataFrame(model_summary_data, columns=['', ''])
        
        # Process core_atp_gapfilling_data and gf_filter_data into a list of dictionaries
        gapfilling_list = []
        for media in core_atp_gapfilling_data:
            core_atp_gapfilling_media = core_atp_gapfilling_data[media]
            row = {
                'no of gapfilled reactions': int(core_atp_gapfilling_media['score']),
                'media': media,
                'ATP Production': f"{round(selected_media_data.get(media, 0), 2):.2f}" if media in selected_media_data else '',
                'gapfilled reactions': '',
                'reversed reaction by gapfilling': '',
                'Filtered Reactions': '',
            }
            if 'new' in core_atp_gapfilling_media:
                gapfilled_reactions = core_atp_gapfilling_media['new']
                if gapfilled_reactions:
                    reactions = [f'<a href="https://modelseed.org/biochem/reactions/{rxn[:-3]}">{rxn}</a> : {direction}' if not rxn.startswith("EX") else f'EX_<a href="https://modelseed.org/biochem/compounds/{rxn[3:-3]}">{rxn}</a> : {direction}' for rxn, direction in gapfilled_reactions.items()]
                    row['gapfilled reactions'] = ' | '.join(reactions)
            if 'failed' in core_atp_gapfilling_media and core_atp_gapfilling_media['failed']:
                row['gapfilled reactions'] = 'Failed'
            if 'reversed' in core_atp_gapfilling_media:
                reversed_reactions = core_atp_gapfilling_media['reversed']
                if reversed_reactions:
                    reactions = [f'<a href="https://modelseed.org/biochem/reactions/{rxn[:-3]}">{rxn}</a> : {direction}' if not rxn.startswith("EX") else f'EX_<a href="https://modelseed.org/biochem/compounds/{rxn[3:-3]}">{rxn}</a> : {direction}' for rxn, direction in reversed_reactions.items()]
                    row['reversed reaction by gapfilling'] = ' | '.join(reactions)
            if media in gf_filter_data:
                gf_filter_media_data = gf_filter_data[media]
                atp_production_values = list(gf_filter_media_data.values())
                if atp_production_values:
                    atp_prod_reaction_pairs = list(atp_production_values[0].items())
                    if atp_prod_reaction_pairs:
                        _, reactions = atp_prod_reaction_pairs[0]
                        if reactions:
                            filtered_reactions = ' | '.join([f'<a href="https://modelseed.org/biochem/reactions/{rxn.split("_")[0]}">{rxn}</a> : {list(value.keys())[0]}' if not rxn.startswith("EX") else f'EX_<a href="https://modelseed.org/biochem/compounds/{rxn[3:-3]}">{rxn}</a> : {list(value.keys())[0]}' for rxn, value in reactions.items()])
                            row['Filtered Reactions'] = filtered_reactions if filtered_reactions else ''
                        if not row['reversed reaction by gapfilling']:
                            row['reversed reaction by gapfilling'] = ''
            gapfilling_list.append(row)
        
        
        
        gapfilling_df = pd.DataFrame(gapfilling_list, columns=['no of gapfilled reactions', 'media', 'ATP Production', 'gapfilled reactions', 'reversed reaction by gapfilling', 'Filtered Reactions'])
        gapfilling_df['no of gapfilled reactions'] = pd.to_numeric(gapfilling_df['no of gapfilled reactions'])
        gapfilling_df = gapfilling_df.sort_values('no of gapfilled reactions')
        
        
        reaction_names = {}
        for rxn in modelutl.model.reactions:
            reaction_id = rxn.id
            reaction_name = rxn.name
            reaction_names[reaction_id] = reaction_name
        
        # Gapfillings Analysis DataFrame
        gapfillings_list = []
        if gf_sensitivity_data:
            for media, media_data in gf_sensitivity_data.items():
                for target, target_data in media_data.items():  # Iterate through each target for the current media
                    for status, status_data in target_data.items():
                        if isinstance(status_data, dict):
                            for reaction_id, reaction_directions in status_data.items():
                                for direction, gapfilling_sensitivity in reaction_directions.items():
                                    if status == 'success':
                                        if isinstance(gapfilling_sensitivity, list):
                                            gapfilling_sensitivity = ', '.join(gapfilling_sensitivity)
                                        gapfillings_list.append({
                                            'Reaction ID': reaction_id,
                                            'Reaction Name': reaction_names.get(reaction_id, ''),  # Get reaction name from the dictionary
                                            'Media': media,
                                            'Direction': direction,
                                            'Target': target,
                                            'Gapfilling Sensitivity': gapfilling_sensitivity
                                        })
                        else:
                            # Handle cases where status_data is null
                            gapfillings_list.append({
                                'Reaction ID': '',  # No data available for Reaction ID
                                'Reaction Name': '',  # No data available for Reaction Name
                                'Media': media,
                                'Direction': '',  # No data available for Direction
                                'Target': target,
                                'Gapfilling Sensitivity': 'Failed Before Filtering' if status == 'FBF' else 'Failed After Filtering' if status == 'FAF' else status  # Status is the 'FBF' or other labels in this case
                            })
        
        gapfillings_analysis_df = pd.DataFrame(gapfillings_list, columns=['Reaction ID', 'Reaction Name', 'Media', 'Direction', 'Target', 'Gapfilling Sensitivity'])
        
        
        # Define the custom color mapping function
        def color_gradient(val):
            if val == 0:
                return 'background-color: green'
            else:
                color_map = cm.get_cmap('YlOrRd')  # Choose the color map
                norm_val = val / gapfilling_df['no of gapfilled reactions'].max()  # Normalize the value between 0 and 1
                color = color_map(norm_val)
                r, g, b, _ = color
                return f'background-color: rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})'
        
        # Apply the default style to the Model Summary DataFrame
        model_summary_df_styled = (
            model_summary_df.style
            .hide_index()
            .set_table_styles([
                {'selector': 'th', 'props': [('border', 'none'), ('background-color', 'white'), ('font-family', 'Oxygen'), ('font-size', '14px'), ('line-height', '20px')]},
                {'selector': 'td', 'props': [('border', 'none'), ('font-family', 'Oxygen'), ('font-size', '14px'), ('line-height', '20px')]},
                {'selector': 'tr:nth-child(even)', 'props': [('background-color', 'white')]},
                {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#f2f2f2')]},
            ])
        )
        
            
        # Apply the default style to the Gapfillings Analysis DataFrame
        # Apply the default style to the Gapfillings Analysis DataFrame
        gapfillings_analysis_df_styled = (
            gapfillings_analysis_df.style
            .hide_index()
            .format({
                'Reaction ID': lambda x: f'<a href="https://modelseed.org/biochem/reactions/{x[:-3]}">{x}</a>' if not x.startswith("EX") else f'<a href="https://modelseed.org/biochem/compounds/{x[3:-3]}">{x}</a>',  # Add hyperlink to Reaction ID
                'Gapfilling Sensitivity': lambda x: ', '.join([f'<a href="https://modelseed.org/biochem/compounds/{i[:-3]}">{i}</a>' for i in x.split(', ')]) if x and not x.startswith('Failed') else x  # Add hyperlinks to Gapfilling Sensitivity
            })
            .set_table_styles([
                {'selector': 'th', 'props': [('border', 'none'), ('background-color', 'white'), ('font-family', 'Oxygen'), ('font-size', '14px'), ('line-height', '20px')]},
                {'selector': 'td', 'props': [('border', 'none'), ('font-family', 'Oxygen'), ('font-size', '14px'), ('line-height', '20px')]},
                {'selector': 'tr:nth-child(even)', 'props': [('background-color', 'white')]},
                {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#f2f2f2')]},
            ])
        )
        
        
        # Apply the default style with alternating row colors, Oxygen font, adjusted font size and line height,
        # and switched order of light grey and white backgrounds in the header column for Core ATP Gapfilling Analysis
        gapfilling_df_styled = (
            gapfilling_df.style
            .applymap(color_gradient, subset=['no of gapfilled reactions'])
            .hide_index()
            .format({'Filtered Reactions': lambda x: f'{x}'})
            .format({'gapfilled reactions': lambda x: f'{x}'})
            .set_table_styles([
                {'selector': 'th', 'props': [('border', 'none'), ('background-color', 'white'), ('font-family', 'Oxygen'), ('font-size', '14px'), ('line-height', '20px')]},
                {'selector': 'td', 'props': [('border', 'none'), ('font-family', 'Oxygen'), ('font-size', '14px'), ('line-height', '20px')]},
                {'selector': 'tr:nth-child(even)', 'props': [('background-color', 'white')]},
                {'selector': 'tr:nth-child(odd)', 'props': [('background-color', '#f2f2f2')]},
            ])
        )
        
        
        # Legend text for Table 1
        annotations_text_1 = """
            <ul>
                <li><b>Reaction ID:</b> The identifier of the reaction.</li>
                <li><b>Reaction Name:</b> The name of the reaction.</li>
                <li><b>Media:</b> The media used by gap filling.</li>
                <li><b>Direction:</b> The direction of the reaction. Can be ">" for forward, "<" for reverse, or "=" for both directions.</li>
                <li><b>Target:</b> The reaction selected as the objective function target for the gapfilling optimization problem. Targets here can be the model’s biomass reaction, commonly named “bio1” for models created by this app. 
                Alternatively, “rxn00062” (ATP Production) reaction is shown for cases where gapfilling was applied to guarantee ATP production in a given media. 
                When reactions are gapfilled for ATP production, we recommend checking the full Core ATP Analysis in Table 2 below. </li>
                <li><b>Gapfilling Sensitivity:</b> Gapfilling is necessary when compounds in the biomass objective function can not be produced by the model. 
                For each reaction we list the biomass compound(s) that can not be synthesized by the model without gapfilling.
                In cases where gap filling fails there are two possible scenarios:
                1) FBF (failed before filtering) : the gapfilling immediately failed, even before we filtered out the ATP breaking reactions. This means this objective CANNOT be satisfied with the entire current database. 
                2) FAF (failed after filtering): the gapfilling succeeded before filtering, but failed after filtering out reactions that break ATP. This tells you definitively if the ATP filtering caused the gapfilling to fail</li>
            </ul>
        """
        #table 2 intro text
        introductory_text = """
        <p>During model reconstruction, we analyze the genome’s core metabolism draft model (model without gapfilling) to assess energy biosynthesis capabilities. 
        The goal of this analysis is to ensure the core metabolism model is able to produce ATP before we expand the model to the genome-scale. 
        This step is designed to prevent gapfilling from introducing reactions that create energy-generating loops. 
        The tests are conducted on a large collection of minimal conditions, with the goal of simulating the model’s capability to produce energy with different electron donor, electron acceptor, and carbon source combinations.</p>
        <p><u>When the draft model of the core metabolism is capable of producing ATP in at least one of the test media, no gapfilling reactions part of this analysis will be added to the model.</u> While we still report the gapfilling requirements for the test media formulations that fail to produce ATP with that draft core model, we only integrate these solutions in the model when no test media succeeds in producing ATP. 
        In this case, the integrated gap-filling solution(s) will be displayed in “Table 1 - Gapfilling Analysis” above, with the “Target” “rxn00062” (ATP Production) objective function.</p> 
        <p>The goal is to display the test results for all media to provide clues for the metabolic capabilities of the genome(s). When many reactions are required for growth on the SO4 testing media conditions, this could be a good indicator that the organism is not capable of performing sulfate reduction. 
        On the other hand, when only one gapfill reaction is required  for ATP production in a given media, multiple scenarios can be considered. 
        1) Organism(s) can’t grow on test condition, and we correctly did not add the reaction to the model.  2) Possible issue with the source genome annotation missing a specific gene function 3) Possible issue with the model reconstruction database. We hope this data helps make more informed decisions on reactions that may need to be manually curated in the model. 
        In cases where is known from the literature or unpublished experimental results that an organism is capable of producing ATP in a given media condition that requires gapfilling in this analysis, you can use the parameter “Force ATP media” in the reconstruction app to ensure those reactions are integrated into the model.
        .</p>
        """
        # Legend text for Table 2
        annotations_text_2 = """
            <ul>
                <li><b>No. of gapfilled reactions:</b> The number of reactions filled by the gapfilling process.</li>
                <li><b>Media:</b> The media in which the reaction takes place.</li>
                <li><b>ATP Production:</b> ATP production by the core metabolism model.</li>
                <li><b>Gapfilled Reactions:</b> Reactions added during the gapfilling process.</li>
                <li><b>Reversed Reaction by Gapfilling:</b> Reactions that have been reversed during the gapfilling process.</li>
                <li><b>Filtered Reactions:</b> Reactions that have been filtered out during the analysis. When a reaction addition would lead to a large increase in ATP production or an infinite energy loop, we filter that reaction out of the gapfilling database and prevent it from being added to the model.</li>
            </ul>
        """
        # Save the data to HTML with the styled DataFrames and the legends
        directory = dirname(output_path)
        os.makedirs(directory, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write('<h1>Model Summary</h1>')
            f.write(model_summary_df_styled.render(escape=False))
            f.write('<br><br>')
            if gf_sensitivity_data:
                f.write('<h1>Table 1 - Gapfillings Analysis</h1>')
                f.write(gapfillings_analysis_df_styled.render(escape=False))
                f.write(f'<br><br><h3>Legend:</h3>{annotations_text_1}')
            else:
                f.write('Gapfilling was not selected as a parameter during reconstruction of the model. As a result your model may not grow on your media object when running Flux Balance Analysis. You can gapfill your model after reconstruction by using the bew Gapiflling Metabolic Model app curently in beta')
            f.write('<br><br>')
            f.write('<h1>Table 2 - Core ATP Analysis</h1>')
            f.write(gapfilling_df_styled.render(escape=False))
            f.write(f'<br><br><h3>Legend:</h3>{annotations_text_2}')
            f.write(introductory_text)
