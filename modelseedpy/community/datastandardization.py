# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:44:07 2022

@author: Andrew Freiburger
"""
from modelseedpy.community.commhelper import phenotypes
from modelseedpy.core.exceptions import ParameterError
from modelseedpy.core.optlanghelper import isIterable
from modelseedpy.core.fbahelper import FBAHelper
from optlang import Constraint
from optlang.symbolics import Zero
from scipy.constants import hour
from zipfile import ZipFile, ZIP_LZMA
from itertools import chain
from typing import Union, Iterable
from copy import deepcopy
from icecream import ic
# from cplex import Cplex
import logging, json, os, re
from pandas import read_csv, DataFrame, ExcelFile
import numpy as np


import logging
logger = logging.getLogger(__name__)

def isnumber(string):
    try:
        float(string)
    except:
        return False
    return True

def _findDate(string, numerical=False):
    monthNames = ["January", "February", "March", "April", "May", "June", "July",
                  "August", "September", "October", "November", "December"]
    monthNums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    days = list(range(31, 0, -1))  # [f"{num}-" for num in list(range(31,0,-1))]
    years = list(range(2010, 2025))+list(range(10,25))  # [f"-{num}" for num in list(range(2000, 2100))]
    americanDates = [f"{mon}-{day}-{year}" for mon in monthNums for day in days for year in years]

    for date in americanDates:
        if re.search(date, string):
            month, day, year = date.split("-")
            if numerical:
                return "-".join([day, month, year])
            return f"{monthNames[int(month)-1][:3]} {day}, {year}"
    # # determine the month
    # for monName in monthNames:
    #     if re.search(monName, string):
    #         month = monName
    #         break
    # if not month:
    #     for monNum in monthNums:
    #         if re.search(monNum, string):
    #             month = monNum  # maybe should be converted to the Name for standardization
    # # determine the day
    # for dayNum in days:
    #     if re.search(dayNum, string):
    #         day = dayNum
    #         break
    # # determine the year
    # for yearNum in years:
    #     if re.search(yearNum, string):
    #         year = yearNum
    #         break
    # return day+month+year

def dict_keys_exists(dic, *keys):
    if keys[0] in dic:
        remainingKeys = keys[1:]
        if len(remainingKeys) > 0:
            dict_keys_exists(dic[keys[0]], keys[1:])
        return True
    return False

def find_dic_number(dic):
    for k, v in dic.items():
        if isnumber(v):
            return v
        num = find_dic_number(dic[k])
    return num

def default_dict_values(dic, key, default):
    return default if not key in dic else dic[key]

def trial_contents(short_code, indices_tup, values):
    matches = [ele == short_code for ele in indices_tup]
    return np.array(values)[matches]

def _spreadsheet_extension_load(path):
    if ".csv" in path:
        return read_csv(path)
    elif ".xls" in path:
        return ExcelFile(path)

def _spreadsheet_extension_parse(path, raw_data, org_sheet):
    if ".csv" in path:
        return raw_data
    elif ".xls" in path:
        return raw_data.parse(org_sheet)

def _met_id_parser(met):
    met_id = re.sub('(\_\w\d+)', '', met)
    met_id = met_id.replace('EX_', '', 1)
    met_id = met_id.replace('c_', '', 1)
    return met_id

def _column_reduction(org_df):
    dataframe = org_df.copy()  # this prevents an irrelevant warning from pandas
    dataframe.columns = map(str, dataframe.columns)
    dataframe.index = dataframe['Well']
    dataframe.drop('Well', axis=1, inplace=True)
    for col in dataframe.columns:
        if any([x in col for x in ['Plate', 'Well', 'Cycle']]):
            dataframe.drop(col, axis=1, inplace=True)
    dataframe.columns = list(map(int, list(map(float, dataframe.columns))))
    return dataframe

def _remove_trials(org_df, ignore_trials, signal, name, significant_deviation):
    # refine the ignore_trials parameter
    if isinstance(ignore_trials, dict):
        ignore_trials['columns'] = list(map(str, ignore_trials['columns'])) if 'columns' in ignore_trials else []
        ignore_trials['rows'] = list(map(str, ignore_trials['rows'])) if 'rows' in ignore_trials else []
        ignore_trials['wells'] = ignore_trials['wells'] if 'wells' in ignore_trials else []
    elif isIterable(ignore_trials):
        if ignore_trials[0][0].isalpha() and isnumber(ignore_trials[0][1:]):
            short_code = True  # TODO - drop trials with respect to the short codes, and not the full codes

    dataframe = org_df.copy()  # this prevents an irrelevant warning from pandas
    dropped_trials = []
    for trial in dataframe.index:
        if isinstance(ignore_trials, dict) and any(
                [trial[0] in ignore_trials['rows'], trial[1:] in ignore_trials['columns'], trial in ignore_trials['wells']]
        ) or isIterable(ignore_trials) and trial in ignore_trials:
            dataframe.drop(trial, axis=0, inplace=True)
            dropped_trials.append(trial)
        elif isIterable(ignore_trials) and trial in ignore_trials:
            dataframe.drop(trial, axis=0, inplace=True)
            dropped_trials.append(trial)
    removed_trials = []
    if 'OD' not in signal:
        for trial, row in dataframe.iterrows():
            row_array = np.array(row.to_list())
            ## remove trials for which the biomass growth did not change by the determined minimum deviation
            if row_array[-1] / row_array[0] < significant_deviation:
                dataframe.drop(trial, axis=0, inplace=True)
                removed_trials.append(trial)
        if removed_trials:
            print(f'The {removed_trials} trials were removed from the {name} measurements, '
                  f'with their deviation over time being less than the threshold of {significant_deviation}.')
    if dropped_trials:
        print(f'The {dropped_trials} trials were dropped from the {name} measurements '
              'per the ignore_trials parameter.')
    return dataframe, dropped_trials+removed_trials

def _check_plateau(org_df, signal, name, significant_deviation, timesteps_len):
    significant_deviation = max([2, significant_deviation])
    dataframe = org_df.copy()  # this prevents an irrelevant warning from pandas
    dropped = []
    for trial, row in dataframe.iterrows():
        row_array = np.array(row.to_list())
        values = []
        tracking = False
        ## remove trials for which the biomass growth did not change by the determined minimum deviation
        for index, val in enumerate(row_array):
            if val / row_array[0] >= significant_deviation or tracking:
                tracking = True
                values.append(val)
                if len(values) > timesteps_len:
                    del values[0]
                remaining_values = list(dataframe.columns[index-timesteps_len+1:])
                if all([len(values) == timesteps_len, values[-1] <= values[0],
                        remaining_values[0] <= remaining_values[-1]*1.1]):
                    # the entire plateau, minus the first point of plateau, are removed
                    dropped = remaining_values
                    break
        if dropped:
            break
    if dropped:
        content = f"{name} {signal}" if name != signal else signal
        print(f"The {dropped} timesteps (with {row_array[index-len(values)+1:]} values) were removed "
              f"from the {content} data since the OD plateaued and is no longer valid.")
    return dropped

def _remove_timesteps(org_df, ignore_timesteps, name, signal):
    dataframe = org_df.copy()  # this prevents an irrelevant warning from pandas
    if ignore_timesteps:
        dropped = []
        for col in dataframe:
            if col in ignore_timesteps:
                dataframe.drop(col, axis=1, inplace=True)
                dropped.append(col)
        if dropped == ignore_timesteps:
            print(f"The ignore_timesteps columns were dropped for the {name} {signal} data.")
        else:
            raise ParameterError(f"The ignore_timesteps values {ignore_timesteps} "
                                 f"were unsuccessfully dropped for the {name} {signal} data.")
    return dataframe, ignore_timesteps

def _df_construction(name, df_name, ignore_trials, ignore_timesteps,
                     significant_deviation, dataframe, row_num, buffer_col1=True):
    # refine the DataFrames
    time_df = _column_reduction(dataframe.iloc[0::2])
    values_df = _column_reduction(dataframe.iloc[1::2])
    # display(name, time_df, values_df)

    # remove specified data trials
    if ignore_trials:
        values_df, removed_trials = _remove_trials(
            values_df, ignore_trials, df_name, name, significant_deviation)
        for row in removed_trials:
            time_df.drop(row, axis=0, inplace=True)

    # remove specified data timesteps
    if ignore_timesteps:
        values_df, removed_timesteps = _remove_timesteps(
            values_df, ignore_timesteps, name, df_name)
        for col in list(map(int, removed_timesteps)):
            time_df.drop(col, axis=1, inplace=True)

    # remove undefined trials
    if buffer_col1:
        possible_rows = [chr(ord("A")+row) for row in range(1, row_num+1)]
        for trial_code in values_df.index:
            if trial_code[0] not in possible_rows:
                values_df.drop(trial_code, axis=0, inplace=True)
                time_df.drop(trial_code, axis=0, inplace=True)

    # process the data for subsequent operations and optimal efficiency
    values_df.astype(str); time_df.astype(str)
    return time_df, values_df

def _find_culture(string):
    matches = re.findall(r"([A-Z]{2}\+?[A-Z]*)", string)
    return [m for m in matches if not any([x in m for x in ["BIOLOG", "III"]])]

def reverse_strip_comp(ID):
    return ID.replace("~", "-")

def _process_csv(self, csv_path, index_col):
    self.zipped_output.append(csv_path)
    csv = read_csv(csv_path) ; csv.index = csv[index_col]
    csv.drop(index_col, axis=1, inplace=True)
    csv.astype(str)
    return csv

def add_rel_flux_cons(model, ex, phenoRXN, carbon_ratio, rel_flux=0.2):
    # {ex.id}_uptakeLimit: {net_{carbonous_ex}} >= {net_{carbon_source}}*{rel_flux}*{carbon_ratio}
    #  The negative flux sign of influxes specifies that the carbon_source value must be lesser than the other
    #  carbon influx that is being constrained.
    cons = Constraint(Zero, lb=0, ub=None, name=f"{ex.id}_uptakeLimit")
    model.add_cons_vars(cons)
    cons.set_linear_coefficients({
            ex.forward_variable:1, ex.reverse_variable:-1,
            phenoRXN.forward_variable:-rel_flux*carbon_ratio, phenoRXN.reverse_variable:rel_flux*carbon_ratio})
    return model, cons


class GrowthData:

    @staticmethod
    def process(community_members: dict, base_media=None, solver: str = 'glpk', all_phenotypes=True,
                data_paths: dict = None, species_abundances: str = None, carbon_conc_series: dict = None,
                ignore_trials: Union[dict, list] = None, ignore_timesteps: list = None, species_identities_rows=None,
                significant_deviation: float = 2, extract_zip_path: str = None, determine_requisite_biomass=False):  #, msdb_path:str=None):
        # define the number of rows in the experimental data
        row_num = len(species_identities_rows)
        if "rows" in carbon_conc_series and carbon_conc_series["rows"]:
            row_num = len(list(carbon_conc_series["rows"].values())[0])
        # load and parse data and metadata
        (media_conc, data_timestep_hr, simulation_time, dataframes, trials, fluxes_df
         ) = GrowthData.load_data(
            base_media, community_members, solver, data_paths, ignore_trials, all_phenotypes,
            ignore_timesteps, significant_deviation, row_num, extract_zip_path)
        experimental_metadata, standardized_carbon_conc, trial_name_conversion = GrowthData.metadata(
            base_media, community_members, species_abundances, carbon_conc_series,
            species_identities_rows, row_num, _findDate(data_paths["path"]))
        data_df = GrowthData.data_process(dataframes, trial_name_conversion)
        requisite_biomass = {} if not determine_requisite_biomass else GrowthData.biomass_growth(
            carbon_conc_series, fluxes_df, data_df.index.unique(), trial_name_conversion,
            data_paths, community_members if all_phenotypes else None)
        return (experimental_metadata, data_df, fluxes_df, standardized_carbon_conc, requisite_biomass,
                trial_name_conversion, np.mean(data_timestep_hr), simulation_time, media_conc)

    @staticmethod
    def load_data(base_media, community_members, solver, data_paths, ignore_trials, all_phenotypes,
                  ignore_timesteps, significant_deviation, row_num, extract_zip_path, min_timesteps=False):
        # define default values
        significant_deviation = significant_deviation or 0
        data_paths = data_paths or {}
        ignore_timesteps = ignore_timesteps or "0:0"
        start, end = ignore_timesteps.split(':')
        raw_data = _spreadsheet_extension_load(data_paths['path'])
        for org_sheet, name in data_paths.items():
            if org_sheet == 'path':
                continue
            df = _spreadsheet_extension_parse(data_paths['path'], raw_data, org_sheet)
            df.columns = df.iloc[6]
            df.drop(df.index[:7], inplace=True)
            ## acquire the default start and end indices of ignore_timesteps
            start = int(start or df.columns[0])
            end = int(end or df.columns[-1])
            break
        ignore_timesteps = list(range(start, end+1)) if start != end else None
        if extract_zip_path:
            with ZipFile(extract_zip_path, 'r') as zp:
                zp.extractall()

        # define only species for which data is defined
        fluxes_df, comm_members = phenotypes(community_members, all_phenotypes, solver=solver)
        modeled_species = list(v for v in data_paths.values() if ("OD" not in v and " " not in v))
        removed_phenotypes = [col for col in fluxes_df if not any([species in col for species in modeled_species])]
        fluxes_df.drop(removed_phenotypes, axis=1, inplace=True)
        if removed_phenotypes:
            print(f'The {removed_phenotypes} phenotypes were removed '
                  f'since their species is not among those with data: {modeled_species}.')

        # determine the time range in which all datasets are significant
        data_timestep_hr = []
        dataframes = {}
        max_timestep_cols = []
        if min_timesteps:
            for org_sheet, name in data_paths.items():
                if org_sheet == 'path' or "OD" in sheet:  continue
                ## define the DataFrame
                sheet = org_sheet.replace(' ', '_')
                df_name = f"{name}:{sheet}"
                dataframes[df_name] = _spreadsheet_extension_parse(data_paths['path'], raw_data, org_sheet)
                dataframes[df_name].columns = dataframes[df_name].iloc[6]
                dataframes[df_name].drop(dataframes[df_name].index[:7], inplace=True)
                ## parse the timesteps from the DataFrame
                drop_timestep_range = GrowthData._min_significant_timesteps(
                    dataframes[df_name], ignore_timesteps, significant_deviation, ignore_trials, df_name, name)
                max_timestep_cols.append(drop_timestep_range)
            ## timesteps that must be dropped for the most restrictive dataset is acquired
            max_cols = max(list(map(len, max_timestep_cols)))
            for ignore_timesteps in max_timestep_cols:
                if len(ignore_timesteps) == max_cols:  break

        # remove trials for which the OD has plateaued
        # TODO - this somehow seems to break when the requisite_biomass is ignored
        for org_sheet, name in data_paths.items():
            if "OD" not in name:  continue
            ## load the OD DataFrame
            sheet = org_sheet.replace(' ', '_')
            df_name = f"{name}:{sheet}"
            dataframes[df_name] = _spreadsheet_extension_parse(data_paths['path'], raw_data, org_sheet)
            dataframes[df_name].columns = dataframes[df_name].iloc[6]
            dataframes[df_name].drop(dataframes[df_name].index[:7], inplace=True)
            ## process the OD DataFrame
            data_times_df, data_values_df = _df_construction(
                name, df_name, ignore_trials, ignore_timesteps,
                significant_deviation, dataframes[df_name], row_num)
            plateaued_times = _check_plateau(data_values_df, name, name, significant_deviation, 3)
            ## define and store the final DataFrames
            for col in plateaued_times:
                if col in data_times_df.columns:    data_times_df.drop(col, axis=1, inplace=True)
                if col in data_values_df.columns:   data_values_df.drop(col, axis=1, inplace=True)
            dataframes[df_name] = (data_times_df, data_values_df)
            break

        # refine the non-OD signals
        for org_sheet, name in data_paths.items():
            if org_sheet == 'path' or "OD" in name:  continue
            sheet = org_sheet.replace(' ', '_')
            df_name = f"{name}:{sheet}"
            if df_name not in dataframes:
                dataframes[df_name] = _spreadsheet_extension_parse(
                    data_paths['path'], raw_data, org_sheet)
                dataframes[df_name].columns = dataframes[df_name].iloc[6]
                dataframes[df_name].drop(dataframes[df_name].index[:7], inplace=True)
            # parse the DataFrame for values
            simulation_time = dataframes[df_name].iloc[0, -1] / hour
            data_timestep_hr.append(simulation_time / int(dataframes[df_name].columns[-1]))
            # define the times and data
            data_times_df, data_values_df = _df_construction(
                name, df_name, ignore_trials, ignore_timesteps, significant_deviation,
                dataframes[df_name], row_num)
            # display(data_times_df) ; display(data_values_df)
            for col in plateaued_times:
                if col in data_times_df.columns:  data_times_df.drop(col, axis=1, inplace=True)
                if col in data_values_df.columns:  data_values_df.drop(col, axis=1, inplace=True)
            dataframes[df_name] = (data_times_df, data_values_df)

        # differentiate the phenotypes for each species
        trials = set(chain.from_iterable([list(times.index) for times, values in dataframes.values()]))
        media_conc = {} if not base_media else {cpd.id: cpd.concentration for cpd in base_media.mediacompounds}
        return (media_conc, data_timestep_hr, simulation_time, dataframes, trials, fluxes_df)

    @staticmethod
    def _min_significant_timesteps(full_df, ignore_timesteps, significant_deviation, ignore_trials, df_name, name):
        # refine the DataFrames
        values_df = _column_reduction(full_df.iloc[1::2])
        values_df, removed_trials = _remove_trials(values_df, ignore_trials, df_name, name, significant_deviation)
        timestep_range = list(set(list(values_df.columns)) - set(ignore_timesteps))
        start, end = ignore_timesteps[0], ignore_timesteps[-1]
        start_index = list(values_df.columns).index(start)
        end_index = list(values_df.columns).index(end)
        ## adjust the customized range such that the threshold is reached.
        for trial, row in values_df.iterrows():
            row_array = np.delete(np.array(row.to_list()), list(range(start_index, end_index + 1)))
            ## remove trials for which the biomass growth did not change by the determined minimum deviation
            while all([row_array[-1] / row_array[0] < significant_deviation,
                       end <= values_df.columns[-1], start >= values_df.columns[0]]):
                # print(timestep_range[0], values_df.columns[0], values_df.columns[-1], end, start)
                if timestep_range[0] == values_df.columns[0] and start != values_df.columns[-1]:
                    timestep_range.append(timestep_range[-1] + 1)
                    start += 1
                    print(f"The end boundary for {name} is increased to {timestep_range[-1]}", end="\r")
                elif timestep_range[-1] == values_df.columns[-1] and end != values_df.columns[0]:
                    timestep_range.append(timestep_range[0] - 1)
                    end -= 1
                    print(f"The start boundary for {name} is decreased to {timestep_range[0]}", end="\r")
                else:
                    raise ParameterError(f"All of the timesteps were omitted for {name}.")
                row_array = np.delete(np.array(row.to_list()), list(range(
                    list(values_df.columns).index(start), list(values_df.columns).index(end) + 1)))
            print("\n")
        return list(range(start, end+1))

    @staticmethod
    def metadata(base_media, community_members, species_abundances,
                 carbon_conc, species_identities_rows, row_num, date):
        # define carbon concentrations for each trial
        carbon_conc = carbon_conc or {}
        carbon_conc['columns'] = default_dict_values(carbon_conc, "columns", {})
        carbon_conc['rows'] = default_dict_values(carbon_conc, "rows", {})
        column_num = len(species_abundances)

        # define the metadata DataFrame and a few columns
        constructed_experiments = DataFrame(index = [f"G{x+1}" for x in list(range(column_num*row_num))])
        constructed_experiments.index.name = "short_code"
        base_media_path = "minimal components media" if not base_media else base_media.path[0]
        constructed_experiments["base_media"] = [base_media_path] * (column_num*row_num)

        # define community content
        # species_mets = {mem["name"]: np.array([mets["consumed"] for mets in mem["phenotypes"].values()]).flatten()
        #                 for mem in community_members.values()}
        # define the strains column
        strains, additional_compounds, experiment_ids = [], [], []
        trial_name_conversion = {}
        count = 1
        ## apply universal values to all trials
        base_row_conc = [] if '*' not in carbon_conc else [
            ':'.join([met, str(carbon_conc['*'][met][0]), str(carbon_conc['*'][met][1])]) for met in carbon_conc['*']]
        members = list(mem["name"] for mem in community_members.values())
        for row in range(1, row_num+1):
            row_conc = base_row_conc[:]
            trial_letter = chr(ord("A") + row)
            trial_name_conversion[trial_letter] = {}
            ## add rows where the initial concentration in the first trial is non-zero
            for met, conc_dict in carbon_conc["rows"].items():
                if conc_dict[sorted(list(conc_dict.keys()))[row-1]] > 0:
                    row_conc.append(':'.join([
                        met, str(conc_dict[sorted(list(conc_dict.keys()))[row-1]]),
                        str(conc_dict[sorted(list(conc_dict.keys()), reverse=True)[-row]])]))

            row_concentration = ';'.join(row_conc)
            composition = {}
            for col in range(1, column_num+1):
                ## construct the columns of information
                additional_compounds.append(row_concentration)
                experiment_id = []
                for member in members:
                    ### define the relative community abundances
                    composition[member] = [member, f"r{species_abundances[col][member]}"]
                    ### define the member strain, where it is appropriate
                    if member in species_identities_rows[row]:
                        composition[member][0] += f"_{species_identities_rows[row][member]}"
                    ### the experimental ID is abundance+memberID
                    if int(composition[member][1][1:]) != 0:
                        experiment_id.append(f"{composition[member][1]}_{composition[member][0]}")
                    composition[member] = ':'.join(composition[member])
                strains.append(';'.join(composition[member] for member in members))
                # for row2 in row_conc:
                    # metID, init, end = row2.split(':')
                    # ### get the met_name for the corresponding match in values
                    # met_name = None
                    # for index, mets in enumerate(species_mets.values()):
                    #     if metID in mets:
                    #         met_name = list(species_mets.keys())[index]
                    #         break
                    # if "met_name" not in locals() or not met_name:
                    #     logger.critical(f"The specified phenotypes {species_mets} for the {members} members"
                    #                     f" does not include the consumption of the available sources"
                    #                     f" {row_conc}; hence, the model cannot grow.")
                    #     content = ""
                    # else:
                    #     content = f"{init}_{met_name}"
                    # experiment_id.append(content)
                experiment_id.extend([":".join(row.split(":")[:2]) for row in row_conc])
                experiment_id = '-'.join(experiment_id)
                experiment_ids.append(experiment_id)
                trial_name_conversion[trial_letter][str(col+1)] = ("G"+str(count), experiment_id)
                count += 1

        # convert the variable concentrations to short codes
        standardized_carbon_conc = {}
        for met, conc in carbon_conc["rows"].items():
            standardized_carbon_conc[met] = {}
            for row, val in conc.items():
                standardized_carbon_conc[met].update({short_code:val for (
                    short_code, expID) in trial_name_conversion[row].values()})
        for met, conc in carbon_conc["columns"].items():
            standardized_carbon_conc[met] = default_dict_values(standardized_carbon_conc, met, {})
            for col, val in conc.items():
                for row in trial_name_conversion:
                    standardized_carbon_conc[met][trial_name_conversion[row][str(col)][0]] = val

        # add columns to the exported dataframe
        constructed_experiments.insert(0, "trial_IDs", experiment_ids)
        constructed_experiments["additional_compounds"] = additional_compounds
        constructed_experiments["strains"] = strains
        constructed_experiments["date"] = [date] * (column_num*row_num)
        constructed_experiments.to_csv("growth_metadata.tsv", sep="\t")
        return constructed_experiments, standardized_carbon_conc, trial_name_conversion

    @staticmethod
    def biomass_growth(carbon_conc, fluxes_df, data_df_trials, trial_name_conversion,
                       data_paths, community_members=None, pheno_info=None):
        # TODO - leverage cFBA to partition metabolite consumption between the defined phenotypes
        pheno_info = pheno_info or {f"{content['name']}_{pheno}": mets
                                    for model, content in community_members.items()
                                    for pheno, mets in content["phenotypes"].items()}
        # invert the trial_name_conversion and data_paths keys and values
        short_code_trials = {contents[0]: row+col for row in trial_name_conversion
                             for col, contents in trial_name_conversion[row].items()}
        # short_code_trials = {contents[0]:contents[1] for contents in trial_name_conversion[row].values()}
        name_signal = {name: signal for signal, name in data_paths.items()}

        # calculate the 90% concentration for each carbon source
        requisite_fluxes = {}
        for trial in [short_code_trials[ID] for ID in data_df_trials]:
            row_letter = trial[0] ; col_number = trial[1:]
            ## add rows where the initial concentration in the first trial is non-zero
            utilized_phenos = {}
            food_gradient = carbon_conc.copy()
            for dimension, content in food_gradient.items():
                for met, conc_dict in content.items():
                    source_conc = conc_dict[row_letter if dimension == "rows" else int(col_number)]
                    # print(met, source_conc)
                    if source_conc == 0 or f"EX_{met}_e0" not in fluxes_df.index:  continue
                    for pheno, val in fluxes_df.loc[f"EX_{met}_e0"].items():
                        # print(pheno, val)
                        if val < 0:  utilized_phenos[pheno] = source_conc*0.9 / val
            total_consumed = sum(list(utilized_phenos.values()))
            # print(utilized_phenos)

            display(fluxes_df)
            short_code = trial_name_conversion[row_letter][col_number][0]
            requisite_fluxes[short_code] = {}
            excreta = {}
            for pheno, flux_conversion in utilized_phenos.items():
                species, phenotype = pheno.split("_", 1)
                fluxes = fluxes_df.loc[:, pheno]*abs(flux_conversion) * abs(flux_conversion/total_consumed)
                requisite_fluxes[short_code][f"{species}|{name_signal[species]}"] = fluxes[fluxes != 0]
                pheno = reverse_strip_comp(pheno)
                if "excreted" in pheno_info[pheno]:
                    # print(pheno_info[pheno]["excreted"])
                    excreta.update({met:fluxes.loc[met] for met in pheno_info[pheno]["excreted"]})
            ## determine the fluxes for the other members of the community through cross-feeding
            participated_species = []
            for pheno, mets in pheno_info.items():
                species, phenotype = pheno.split("_", 1)
                if any([species in ph for ph in utilized_phenos]) or species in participated_species:  continue
                for met in mets["consumed"]:
                    exMet = f"EX_{met}_e0"
                    if exMet not in excreta:  continue
                    fluxes = abs(excreta[exMet] * 0.99 / fluxes_df.loc[exMet, pheno]) * fluxes_df.loc[:, pheno]
                    requisite_fluxes[short_code][f"{species}|{name_signal[species]}"] = fluxes[fluxes != 0]
                    participated_species.append(species)
        # print(requisite_fluxes)
        return requisite_fluxes

    @staticmethod
    def data_process(dataframes, trial_name_conversion):
        short_codes, trials_list = [], []
        values, times = {}, {}  # The times must capture upstream
        first = True
        for df_name, (times_df, values_df) in dataframes.items():
            # print(df_name)
            # display(times_df) ; display(values_df)
            times_tup = FBAHelper.parse_df(times_df)
            average_times = np.mean(times_tup.values, axis=0)
            values[df_name], times[df_name] = [], []
            for trial_code in values_df.index:
                row_let, col_num = trial_code[0], trial_code[1:]
                # print(trial_code, row_let, col_num)
                for trial_row_values in trial_contents(trial_code, values_df.index, values_df.values):
                    if first:
                        short_code, experimentalID = trial_name_conversion[row_let][col_num]
                        trials_list.extend([experimentalID] * len(values_df.columns))
                        short_codes.extend([short_code] * len(values_df.columns))
                    values[df_name].extend(trial_row_values)
                    times[df_name].extend(average_times)
            first = False
        # process the data to the smallest dataset, to accommodate heterogeneous data sizes
        minVal = min(list(map(len, values.values())))
        for df_name, data in values.items():
            values[df_name] = data[:minVal]
        times2 = times.copy()
        for df_name, data in times2.items():
            times[df_name] = data[:minVal]
        # construct the growth DataFrame
        df_data = {"trial_IDs": trials_list[:minVal], "short_codes": short_codes[:minVal]}
        df_data.update({"Time (s)": np.mean(list(times.values()), axis=0)})  # element-wise average
        df_data.update({df_name:vals for df_name, vals in values.items()})
        data_df = DataFrame(df_data)
        data_df.index = data_df["short_codes"]
        data_df = data_df.drop(["short_codes"], axis=1)
        data_df.to_csv("growth_spectra.tsv", sep="\t")
        return data_df


class BiologData:

    @staticmethod
    def process(data_paths, trial_conditions_path, community_members, col_row_num, member_conversions,
                culture=None, date=None, significant_deviation=None, solver="glpk", msdb_path:str=None):
        row_num = 8 ; column_num = 12
        (zipped_output, data_timestep_hr, simulation_time, dataframes, trials, culture, date, fluxes_df
         ) = BiologData.load_data(data_paths, significant_deviation, community_members,
                                  col_row_num, row_num, culture, date, solver)
        experimental_metadata, standardized_carbon_conc, trial_name_conversion = BiologData.metadata(
            trial_conditions_path, row_num, column_num, culture, date)
        biolog_df = BiologData.data_process(dataframes, trial_name_conversion)
        requisite_biomass = BiologData.biomass_growth(biolog_df, member_conversions)
        return (experimental_metadata, biolog_df, fluxes_df, standardized_carbon_conc, requisite_biomass,
                trial_name_conversion, np.mean(data_timestep_hr), simulation_time)

    @staticmethod
    def load_data(data_paths, significant_deviation, community_members, col_row_num,
                  row_num, culture, date, solver):
        zipped_output = [data_paths['path'], "fluxes.tsv"]
        # determine the metabolic fluxes for each member and phenotype
        # import and parse the raw CSV data
        # TODO - this may be capable of emulating leveraged functions from the GrowthData object
        fluxes_df = phenotypes(community_members, solver=solver)
        # fluxes_df = None
        data_timestep_hr = []
        dataframes = {}
        raw_data = _spreadsheet_extension_load(data_paths['path'])
        significant_deviation = significant_deviation or 2
        # culture = culture or _find_culture(data_paths['path'])
        culture = culture or ",".join([x for x in data_paths.values() if (x not in ["OD"] and not re.search(r"\w\.\w", x))])
        date = date or _findDate(data_paths['path'])
        for org_sheet, name in data_paths.items():
            if org_sheet == 'path':
                continue
            sheet = org_sheet.replace(" ", "_")
            df_name = f"{name}:{sheet}"
            if df_name not in dataframes:
                dataframes[df_name] = _spreadsheet_extension_parse(
                    data_paths['path'], raw_data, org_sheet)
                dataframes[df_name].columns = dataframes[df_name].iloc[col_row_num]
                dataframes[df_name].drop(dataframes[df_name].index[:col_row_num+1], inplace=True)
                dataframes[df_name].dropna(inplace=True)
            # parse the DataFrame for values
            dataframes[df_name].columns = [str(x).strip() for x in dataframes[df_name].columns]
            simulation_time = dataframes[df_name].iloc[0, -1] / hour
            # display(dataframes[df_name])
            data_timestep_hr.append(simulation_time / int(float(dataframes[df_name].columns[-1])))
            # define the times and data
            data_times_df, data_values_df = _df_construction(
                name, df_name, None, None, significant_deviation,
                dataframes[df_name], row_num, False)
            # display(data_times_df) ; display(data_values_df)
            dataframes[df_name] = (data_times_df, data_values_df)

        # differentiate the phenotypes for each species
        trials = set(chain.from_iterable([list(df.index) for df, times in dataframes.values()]))
        return (zipped_output, data_timestep_hr, simulation_time, dataframes, trials, culture, date, fluxes_df)

    @staticmethod
    def metadata(trial_conditions_path, row_num, column_num, culture, date):
        # define the conditions for each trial
        with open(trial_conditions_path) as trials:
            trial_conditions = json.load(trials)

        # define the metadata DataFrame and a few columns
        constructed_experiments = DataFrame()
        ex_prefix = "B"
        constructed_experiments.index = [f"{ex_prefix}{x+1}" for x in list(range(row_num*column_num))]
        constructed_experiments.index.name = "short_code"

        # define the strains column
        experiment_ids, trial_names = [], []
        trial_name_conversion, trial_mets = {}, {}
        count = 1
        ## apply universal values to all trials
        for row in range(row_num):
            trial_letter = chr(ord("A") + row)
            trial_name_conversion[trial_letter] = {}
            ## add rows where the initial concentration in the first trial is non-zero
            for col in range(1, column_num+1):
                ## construct the columns of information
                dataID = trial_letter+str(col)
                MSID = trial_conditions[dataID]["ModelSEED_ID"]
                short_code = ex_prefix+str(count)

                experiment_ids.append(MSID)
                trial_names.append(trial_conditions[dataID]["name"])
                trial_name_conversion[trial_letter][str(col)] = (short_code, MSID)
                trial_mets[MSID] = {short_code:trial_conditions[dataID]["mM"]}
                count += 1

        # add columns to the exported dataframe
        constructed_experiments.insert(0, "ModelSEED_ID", experiment_ids)
        constructed_experiments.insert(0, "condition", trial_names)
        constructed_experiments["strain"] = [culture] * (column_num*row_num)
        constructed_experiments["date"] = [date] * (column_num*row_num)
        constructed_experiments.to_csv("growth_metadata.tsv", sep="\t")
        return constructed_experiments, trial_mets, trial_name_conversion

    @staticmethod
    def data_process(dataframes, trial_name_conversion):
        short_codes, trials_list = [], []
        values, times = {}, {}  # The times must capture upstream
        first = True
        for df_name, (times_df, values_df) in dataframes.items():
            # display(df_name, times_df, values_df)
            times_tup = FBAHelper.parse_df(times_df)
            # display(DataFrame(times_tup.values))
            average_times = list(np.mean(times_tup.values, axis=0))
            # print(average_times)
            # print(len(average_times))
            values[df_name], times[df_name] = [], []
            for exprID in values_df.index:
                row_let, col_num = exprID[0], exprID[1:]
                for trial_row_values in trial_contents(exprID, values_df.index, values_df.values):
                    if first:
                        short_code, experimentalID = trial_name_conversion[row_let][col_num]
                        trials_list.extend([experimentalID] * len(values_df.columns))
                        short_codes.extend([short_code] * len(values_df.columns))
                    if len(trial_row_values) != len(average_times):
                        print(f"The length of the trial data {len(trial_row_values)} "
                              f"exceeds that of the timesteps {len(average_times)} "
                              f"which creates an incompatible DataFrame.")
                    values[df_name].extend(trial_row_values)
                    times[df_name].extend(average_times)
            first = False
        # process the data to the smallest dataset, to accommodate heterogeneous data sizes
        minVal = min(list(map(len, values.values())))
        for df_name, data in values.items():
            values[df_name] = data[:minVal]
        times2 = times.copy()
        for df_name, data in times2.items():
            times[df_name] = data[:minVal]
        df_data = {"trial_IDs": trials_list, "short_codes": short_codes}
        df_data.update({"Time (s)": list(np.mean(list(times.values()), axis=0))})  # element-wise average
        df_data.update({df_name:vals for df_name, vals in values.items()})
        biolog_df = DataFrame(df_data)
        biolog_df.index = biolog_df["short_codes"]
        del biolog_df["short_codes"]
        biolog_df.to_csv("growth_spectra.tsv", sep="\t")

        return biolog_df

    @staticmethod
    def biomass_growth(biolog_df, member_conversions):
        requisite_biomass = {}
        for short_code in biolog_df.index.unique():
            requisite_biomass[short_code] = {}
            for signal, conversion in member_conversions.items():
                short_code_df = biolog_df[biolog_df.index == short_code]
                requisite_biomass[short_code][signal] = conversion * short_code_df[
                    signal.replace("|", ":").replace(" ", "_")].iloc[-1]
        return requisite_biomass
