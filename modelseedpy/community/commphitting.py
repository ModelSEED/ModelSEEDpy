# -*- coding: utf-8 -*-
# from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.exceptions import FeasibilityError, ParameterError, ObjectAlreadyDefinedError, NoFluxError
from modelseedpy.core.optlanghelper import OptlangHelper, Bounds, tupVariable, tupConstraint, tupObjective, isIterable, define_term
from modelseedpy.community.datastandardization import GrowthData
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.biochem import from_local
from scipy.constants import hour, minute
from zipfile import ZipFile, ZIP_LZMA
from optlang import Model, Objective
from time import sleep, process_time
from typing import Union, Iterable
from optlang.symbolics import Zero
from scipy.optimize import newton
from matplotlib import pyplot
from math import inf, isclose
from deepdiff import DeepDiff
from pandas import DataFrame
from itertools import chain
from pprint import pprint
from h5py import File
from icecream import ic
import numpy as np
import cobra.io
# from cplex import Cplex
import warnings, logging, json, os, re

logger = logging.getLogger(__name__)

def dict_keys_exists(dic, *keys):
    result = keys[0] in dic
    if keys[0] in dic:
        remainingKeys = keys[1:]
        if len(remainingKeys) > 0:  result = dict_keys_exists(dic[keys[0]], *remainingKeys)
        return result
    return result

def find_dic_number(dic):
    for k, v in dic.items():
        if FBAHelper.isnumber(v):  return v
        num = find_dic_number(dic[k])
    return num

def trial_contents(short_code, indices_tup, values):
    matches = [ele == short_code for ele in indices_tup]
    return np.array(values)[matches]

def dic_keys(dic):
    keys = []
    if isinstance(dic, dict):
        for key, value in dic.items():
            keys.append(key)
            keys.extend(dic_keys(value))
    return keys

# define data objects
def _name(name, suffix, short_code, timestep, names):
    name = '-'.join([x for x in list(map(str, [name + suffix, short_code, timestep])) if x])
    if name not in names:  names.append(name) ; return name
    else:  pprint(names) ; raise ObjectAlreadyDefinedError(f"The object {name} is already defined for the problem.")

def _export_model_json(json_model, path):
    with open(path, 'w') as lp:  json.dump(json_model, lp, indent=3)

def _met_id_parser(met):
    met_id = re.sub('(\_\w\d+)', '', met)
    met_id = met_id.replace('EX_', '', 1)
    met_id = met_id.replace('c_', '', 1)
    return met_id

# define an entity as a variable or a constant
def _obj_val(primal, name, pheno, short_code, timestep, bounds, data_timestep_hr, names):
    time_hr = int(timestep) * data_timestep_hr
    return tupVariable(_name(name, pheno, short_code, timestep, names),
                       Bounds=bounds) if not primal else primal[short_code][name+pheno][time_hr]

def _michaelis_menten(conc, vmax, km):
    return (conc*vmax)/(km+conc)

def clamp(val, minimum, maximum):
    return min(max(val, minimum), maximum)

# parse primal values for use in the optimization loops
def parse_primals(primal_values, entity_labels=None, coefs=None, kcat_vals=None):
    if kcat_vals:
        kcat_primal = {}
        for trial, content in primal_values.items():
            for primal, time_value in content.items():
                if "bin" not in primal:  continue
                name, trial = primal.split("-")
                number = re.search(r"(\d)", name).group()
                species, pheno = re.sub(r"(bin\d_)", "", name).split("_")
                if "stationary" in pheno:  continue
                if species not in kcat_primal:  kcat_primal[species] = {}
                if pheno not in kcat_primal[species]:  kcat_primal[species][pheno] = 0
                # kcat_(k,new) = sum_z^Z ( kcat_z * bin_k^z ) * kcat_(k,old) < 10
                if time_value == 0 and kcat_primal[species][pheno] < 10:
                    kcat_primal[species][pheno] += coefs[int(number)-1]*kcat_vals[species][pheno]
                kcat_primal[species][pheno] = clamp(kcat_primal[species][pheno], 1e-4, 10)
        return kcat_primal
    select_primals = {}
    for trial, entities in primal_values.items():
        select_primals[trial] = {}
        for entity, times in entities.items():
            # a poor man's dictionary copy
            if any([label in entity for label in entity_labels]):  select_primals[trial][entity] = dict(list(times.items()))
    return select_primals

def signal_species(signal):
    return signal.split(":")[0].replace(" ", "_")

def _partition_coefs(initial_val, divisor):
    return (initial_val, initial_val/divisor, initial_val/divisor**2, initial_val/divisor**3, initial_val/divisor**4)


biomass_partition_coefs = [_partition_coefs(10, 10), _partition_coefs(2, 2), _partition_coefs(1, 3)]


class CommPhitting:

    def __init__(self, msdb_path, community_members: dict=None, fluxes_df=None, data_df=None, carbon_conc=None,
                 media_conc=None, experimental_metadata=None, base_media=None, solver: str = 'glpk', all_phenotypes=True,
                 data_paths: dict = None, species_abundances: str = None, ignore_trials: Union[dict, list] = None,
                 ignore_timesteps: list = None, species_identities_rows=None, significant_deviation: float = 2,
                 extract_zip_path: str = None, determine_requisite_biomass:bool = True, consumed_mets:iter=None):
        self.msdb = from_local(msdb_path) ; self.msdb_path = msdb_path
        self.solver = solver ; self.all_phenotypes = all_phenotypes ; self.data_paths = data_paths
        self.species_abundances = species_abundances ; self.ignore_trials = ignore_trials
        self.ignore_timesteps = ignore_timesteps ; self.species_identities_rows = species_identities_rows
        self.significant_deviation = significant_deviation ; self.extract_zip_path = extract_zip_path

        self.community_members = community_members
        self.consumed_mets = consumed_mets or set([
            met for content in community_members.values() for met in content["phenotypes"]])
        if community_members is not None or any([x is None for x in [fluxes_df, data_df]]):
            (self.experimental_metadata, data_df, fluxes_df, carbon_conc, self.requisite_biomass,
             self.trial_name_conversion, self.data_timestep_hr, simulation_timestep, media_conc
             ) = GrowthData.process(community_members, base_media, solver, all_phenotypes, data_paths,
                                    species_abundances, carbon_conc, ignore_trials, ignore_timesteps,
                                    species_identities_rows, significant_deviation, extract_zip_path,
                                    determine_requisite_biomass)
            #                       for content in community_members.values() for met in content["phenotypes"]]
        self.fluxes_tup = FBAHelper.parse_df(fluxes_df)
        self.fluxes_df = fluxes_df ; self.data_df = data_df
        self.default_excreta = [index for index, row in fluxes_df.iterrows() if any(row > 1)]
        self.parameters, self.variables, self.constraints = {}, {}, {}
        self.zipped_output, self.plots, self.names = [], [], []
        self.experimental_metadata = experimental_metadata
        self.carbon_conc = carbon_conc; self.media_conc = media_conc

    #################### FITTING PHASE METHODS ####################

    def fit_kcat(self, parameters: dict = None, mets_to_track: list = None, rel_final_conc: dict = None,
                 zero_start: list = None, abs_final_conc: dict = None, graphs: list = None, data_timesteps: dict = None,
                 export_zip_name: str = None, export_parameters: bool = True, requisite_biomass: dict = None,
                 export_lp:str = f'solveKcat.lp', figures_zip_name:str=None, publishing=True, primals_export_path=None):
        if export_zip_name and os.path.exists(export_zip_name):  os.remove(export_zip_name)
        kcat_primal = None
        requisite_biomass = requisite_biomass or self.requisite_biomass
        for index, coefs in enumerate(biomass_partition_coefs):
            # solve for growth rate constants with the previously solved biomasses
            newSim = CommPhitting(self.msdb_path, None, self.fluxes_df, self.data_df, self.carbon_conc,
                                  self.media_conc, self.experimental_metadata, None, self.solver, self.all_phenotypes,
                                  self.data_paths, self.species_abundances, self.ignore_trials, self.ignore_timesteps,
                                  self.species_identities_rows, self.significant_deviation, self.extract_zip_path,
                                  True, self.consumed_mets)
            newSim.define_problem(parameters, mets_to_track, rel_final_conc, zero_start, abs_final_conc,
                                  data_timesteps, export_zip_name, export_parameters, export_lp,
                                  kcat_primal, coefs, requisite_biomass)
            newSim.compute(graphs, export_zip_name, figures_zip_name, publishing,
                           primals_export_path or re.sub(r"(.lp)", ".json", export_lp))
            kcat_primal = parse_primals(newSim.values, coefs=coefs, kcat_vals=newSim.parameters["kcat"])
            pprint(kcat_primal)
            print(f"Interation {index+1} is complete\n")
        kcats = {k: val for k, val in newSim.values.items() if "kcat" in k}
        DataFrame(kcats).T.to_csv("pheno_growth_kcat.tsv", sep="\t")
        return kcats

    def fit(self, parameters:dict=None, mets_to_track: list = None, rel_final_conc:dict=None, zero_start:list=None,
            abs_final_conc:dict=None, graphs: list = None, data_timesteps: dict = None,
            export_zip_name: str = None, export_parameters: bool = True, requisite_biomass: dict = None,
            export_lp: str = 'CommPhitting.lp', figures_zip_name:str=None, publishing:bool=False, primals_export_path=None):
        if hasattr(self, "requisite_biomass"):  requisite_biomass = self.requisite_biomass
        self.define_problem(parameters, mets_to_track, rel_final_conc, zero_start, abs_final_conc,
                            data_timesteps, export_zip_name, export_parameters, export_lp,
                            None, None, requisite_biomass)
        self.compute(graphs, export_zip_name, figures_zip_name, publishing,
                     primals_export_path or re.sub(r"(.lp)", ".json", export_lp))

    def define_b_vars(self, pheno, short_code, timestep, variables):
        self.variables['b_' + pheno][short_code][timestep] = tupVariable(
            _name("b_", pheno, short_code, timestep, self.names), Bounds(0, 1000))
        self.variables['b1_' + pheno][short_code][timestep] = tupVariable(
            _name("b1_", pheno, short_code, timestep, self.names), Bounds(0, 1000))
        self.variables['b2_' + pheno][short_code][timestep] = tupVariable(
            _name("b2_", pheno, short_code, timestep, self.names), Bounds(0, 1000))
        self.variables['b3_' + pheno][short_code][timestep] = tupVariable(
            _name("b3_", pheno, short_code, timestep, self.names), Bounds(0, 1000))
        self.variables['b4_' + pheno][short_code][timestep] = tupVariable(
            _name("b4_", pheno, short_code, timestep, self.names), Bounds(0, 1000))
        self.variables['b5_' + pheno][short_code][timestep] = tupVariable(
            _name("b5_", pheno, short_code, timestep, self.names), Bounds(0, 1000))
        variables.extend([self.variables['b_' + pheno][short_code][timestep],
                          self.variables['b1_' + pheno][short_code][timestep],
                          self.variables['b2_' + pheno][short_code][timestep],
                          self.variables['b3_' + pheno][short_code][timestep],
                          self.variables['b4_' + pheno][short_code][timestep],
                          self.variables['b5_' + pheno][short_code][timestep]])
        if short_code not in self.variables[f"bin1_{pheno}"]:
            self.variables[f"bin1_{pheno}"][short_code] = tupVariable(
                _name("bin1_", pheno, short_code, "", self.names), Bounds(0, 1), "binary")
            self.variables[f"bin2_{pheno}"][short_code] = tupVariable(
                _name("bin2_", pheno, short_code, "", self.names), Bounds(0, 1), "binary")
            self.variables[f"bin3_{pheno}"][short_code] = tupVariable(
                _name("bin3_", pheno, short_code, "", self.names), Bounds(0, 1), "binary")
            self.variables[f"bin4_{pheno}"][short_code] = tupVariable(
                _name("bin4_", pheno, short_code, "", self.names), Bounds(0, 1), "binary")
            self.variables[f"bin5_{pheno}"][short_code] = tupVariable(
                _name("bin5_", pheno, short_code, "", self.names), Bounds(0, 1), "binary")
            variables.extend([self.variables[f"bin1_{pheno}"][short_code], self.variables[f"bin2_{pheno}"][short_code],
                              self.variables[f"bin3_{pheno}"][short_code], self.variables[f"bin4_{pheno}"][short_code],
                              self.variables[f"bin5_{pheno}"][short_code]])
        return variables

    def define_b_cons(self, pheno, short_code, timestep, biomass_coefs):
        biomass_coefs = biomass_coefs or biomass_partition_coefs[-1]
        # define the partitioned biomass groups
        ## b_n{pheno,t} <= coef*b_tot{pheno,t}
        self.constraints['b1c_' + pheno][short_code][timestep] = tupConstraint(
            _name("b1c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    {"elements": [biomass_coefs[0], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b1_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b2c_' + pheno][short_code][timestep] = tupConstraint(
            _name("b2c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    {"elements": [biomass_coefs[1], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b2_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b3c_' + pheno][short_code][timestep] = tupConstraint(
            _name("b3c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    {"elements": [biomass_coefs[2], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b3_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b4c_' + pheno][short_code][timestep] = tupConstraint(
            _name("b4c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    {"elements": [biomass_coefs[3], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b4_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b5c_' + pheno][short_code][timestep] = tupConstraint(
            _name("b5c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    {"elements": [biomass_coefs[4], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b5_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })

        # define the comprehensive biomass constraints
        ## coef*b{pheno,t} - b_n{pheno,t} - 1000*bin_n{pheno} <= 0
        self.constraints['b1c_control_' + pheno][short_code][timestep] = tupConstraint(
            _name("b1c_control_", pheno, short_code, timestep, self.names), Bounds(None, 0), {
                "elements": [
                    {"elements": [biomass_coefs[0], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b1_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin1_{pheno}"][short_code].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b2c_control_' + pheno][short_code][timestep] = tupConstraint(
            _name("b2c_control_", pheno, short_code, timestep, self.names), Bounds(None, 0), {
                "elements": [
                    {"elements": [biomass_coefs[1], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b2_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin2_{pheno}"][short_code].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b3c_control_' + pheno][short_code][timestep] = tupConstraint(
            _name("b3c_control_", pheno, short_code, timestep, self.names), Bounds(None, 0), {
                "elements": [
                    {"elements": [biomass_coefs[2], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b3_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin3_{pheno}"][short_code].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b4c_control_' + pheno][short_code][timestep] = tupConstraint(
            _name("b4c_control_", pheno, short_code, timestep, self.names), Bounds(None, 0), {
                "elements": [
                    {"elements": [biomass_coefs[3], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b4_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin4_{pheno}"][short_code].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })
        self.constraints['b5c_control_' + pheno][short_code][timestep] = tupConstraint(
            _name("b5c_control_", pheno, short_code, timestep, self.names), Bounds(None, 0), {
                "elements": [
                    {"elements": [biomass_coefs[4], self.variables['b_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1, self.variables['b5_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin5_{pheno}"][short_code].name],
                     "operation": "Mul"},
                ],
                "operation": "Add"
            })

        # define the binary constraints
        ## b_n{pheno,t} <= 1000 - 1000*bin_n{pheno}
        self.constraints['bin1c_' + pheno][short_code][timestep] = tupConstraint(
            _name("bin1c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    1000,
                    {"elements": [-1, self.variables['b1_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin1_{pheno}"][short_code].name],
                     "operation": "Mul"}
                ],
                "operation": "Add"
            })
        self.constraints['bin2c_' + pheno][short_code][timestep] = tupConstraint(
            _name("bin2c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    1000,
                    {"elements": [-1, self.variables['b2_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin2_{pheno}"][short_code].name],
                     "operation": "Mul"}
                ],
                "operation": "Add"
            })
        self.constraints['bin3c_' + pheno][short_code][timestep] = tupConstraint(
            _name("bin3c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    1000,
                    {"elements": [-1, self.variables['b3_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin3_{pheno}"][short_code].name],
                     "operation": "Mul"}
                ],
                "operation": "Add"
            })
        self.constraints['bin4c_' + pheno][short_code][timestep] = tupConstraint(
            _name("bin4c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    1000,
                    {"elements": [-1, self.variables['b4_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin4_{pheno}"][short_code].name],
                     "operation": "Mul"}
                ],
                "operation": "Add"
            })
        self.constraints['bin5c_' + pheno][short_code][timestep] = tupConstraint(
            _name("bin5c_", pheno, short_code, timestep, self.names), Bounds(0, None), {
                "elements": [
                    1000,
                    {"elements": [-1, self.variables['b5_' + pheno][short_code][timestep].name],
                     "operation": "Mul"},
                    {"elements": [-1000, self.variables[f"bin5_{pheno}"][short_code].name],
                     "operation": "Mul"}
                ],
                "operation": "Add"
            })

        # load the constraints to the model
        return [self.constraints['b1c_' + pheno][short_code][timestep],
                self.constraints['b2c_' + pheno][short_code][timestep],
                self.constraints['b3c_' + pheno][short_code][timestep],
                self.constraints['b4c_' + pheno][short_code][timestep],
                self.constraints['b5c_' + pheno][short_code][timestep],
                self.constraints['b1c_control_' + pheno][short_code][timestep],
                self.constraints['b2c_control_' + pheno][short_code][timestep],
                self.constraints['b3c_control_' + pheno][short_code][timestep],
                self.constraints['b4c_control_' + pheno][short_code][timestep],
                self.constraints['b5c_control_' + pheno][short_code][timestep],
                self.constraints['bin1c_' + pheno][short_code][timestep],
                self.constraints['bin2c_' + pheno][short_code][timestep],
                self.constraints['bin3c_' + pheno][short_code][timestep],
                self.constraints['bin4c_' + pheno][short_code][timestep],
                self.constraints['bin5c_' + pheno][short_code][timestep]]

    def initialize_vars_cons(self, pheno, short_code):
        # cvt and cvf
        self.variables['cvt_' + pheno] = {}; self.variables['cvf_' + pheno] = {}
        self.variables['cvt_' + pheno][short_code] = {}; self.variables['cvf_' + pheno][short_code] = {}
        # total biomass and growth
        self.variables['b_' + pheno] = {}; self.variables['g_' + pheno] = {}
        self.variables['b_' + pheno][short_code] = {}; self.variables['g_' + pheno][short_code] = {}
        self.constraints['gc_' + pheno] = {}; self.constraints['cvc_' + pheno] = {}
        self.constraints['gc_' + pheno][short_code] = {}; self.constraints['cvc_' + pheno][short_code] = {}
        # partitioned biomasses
        self.variables['b1_' + pheno] = {}; self.variables['b2_' + pheno] = {};  self.variables['b3_' + pheno] = {}
        self.variables['b4_' + pheno] = {}; self.variables['b5_' + pheno] = {}
        self.variables['b1_' + pheno][short_code] = {}; self.variables['b2_' + pheno][short_code] = {}
        self.variables['b3_' + pheno][short_code] = {}; self.variables['b4_' + pheno][short_code] = {}
        self.variables['b5_' + pheno][short_code] = {}
        ## biomass binary variables
        self.variables[f'bin1_{pheno}'] = {}; self.variables[f'bin2_{pheno}'] = {}; self.variables[f'bin3_{pheno}'] = {}
        self.variables[f'bin4_{pheno}'] = {}; self.variables[f'bin5_{pheno}'] = {}
        self.variables[f"bin1_{pheno}"][short_code] = {}; self.variables[f"bin2_{pheno}"][short_code] = {}
        self.variables[f"bin3_{pheno}"][short_code] = {}; self.variables[f"bin4_{pheno}"][short_code] = {}
        self.variables[f"bin5_{pheno}"][short_code] = {}
        ## biomass partition constraints
        self.constraints['b1c_' + pheno] = {}; self.constraints['b2c_' + pheno] = {}; self.constraints['b3c_' + pheno] = {}
        self.constraints['b4c_' + pheno] = {}; self.constraints['b5c_' + pheno] = {}
        self.constraints['b1c_' + pheno][short_code] = {}; self.constraints['b2c_' + pheno][short_code] = {}
        self.constraints['b3c_' + pheno][short_code] = {}; self.constraints['b4c_' + pheno][short_code] = {}
        self.constraints['b5c_' + pheno][short_code] = {}
        self.constraints['b1c_control_' + pheno] = {}; self.constraints['b2c_control_' + pheno] = {}
        self.constraints['b3c_control_' + pheno] = {}; self.constraints['b4c_control_' + pheno] = {}
        self.constraints['b5c_control_' + pheno] = {}
        self.constraints['b1c_control_' + pheno][short_code] = {}; self.constraints['b2c_control_' + pheno][short_code] = {}
        self.constraints['b3c_control_' + pheno][short_code] = {}; self.constraints['b4c_control_' + pheno][short_code] = {}
        self.constraints['b5c_control_' + pheno][short_code] = {}
        self.constraints[f'binc_{pheno}'] = {}; self.constraints[f'binc_{pheno}'][short_code] = {}
        self.constraints['bin1c_' + pheno] = {}; self.constraints['bin2c_' + pheno] = {}
        self.constraints['bin3c_' + pheno] = {}; self.constraints['bin4c_' + pheno] = {}; self.constraints['bin5c_' + pheno] = {}
        self.constraints['bin1c_' + pheno][short_code] = {}; self.constraints['bin2c_' + pheno][short_code] = {}
        self.constraints['bin3c_' + pheno][short_code] = {}; self.constraints['bin4c_' + pheno][short_code] = {}
        self.constraints['bin5c_' + pheno][short_code] = {}

    def get_timestep_bin(self, timestep):
        if timestep < self.first: return 0
        elif timestep < self.second: return 1
        elif timestep < self.third: return 2
        elif timestep < self.fourth: return 3
        return 4

    def define_problem(self, parameters=None, mets_to_track=None, rel_final_conc=None, zero_start=None,
                       abs_final_conc=None, data_timesteps=None, export_zip_name: str=None,
                       export_parameters: bool=True, export_lp: str='CommPhitting.lp', primal_values=None,
                       biomass_coefs=None, requisite_biomass:dict=None, biolog_simulation=False,
                       export_phenotype_profiles=True):
        # parse the growth data
        growth_tup = FBAHelper.parse_df(self.data_df, False)
        self.phenotypes = list(self.fluxes_tup.columns)
        self.phenotypes.extend([signal_species(signal)+"_stationary" for signal in growth_tup.columns if (
                ":" in signal and "OD" not in signal)])
        self.species_list = [signal_species(signal) for signal in growth_tup.columns if ":" in signal]
        num_sorted = np.sort(np.array([int(obj[1:]) for obj in set(growth_tup.index)]))
        # TODO - short_codes must be distinguished for different conditions
        unique_short_codes = [f"{growth_tup.index[0][0]}{num}" for num in map(str, num_sorted)]
        full_times = growth_tup.values[:, growth_tup.columns.index("Time (s)")]
        self.times = {short_code: trial_contents(short_code, growth_tup.index, full_times)
                      for short_code in unique_short_codes}
        average_time_series = np.mean(list(self.times.values()), axis=0)  ;  points = len(average_time_series)
        self.first, self.second, self.third, self.fourth = int(points*0.1), int(points*0.25), int(points*0.45), int(points*0.7)
        self.time_ranges = {0: average_time_series[:self.first], 1: average_time_series[self.first:self.second],
                            2: average_time_series[self.second:self.third], 3: average_time_series[self.third:self.fourth],
                            4: average_time_series[self.fourth:]}

        # define default values
        # TODO render bcv and cvmin dependent upon temperature, and possibly trained on Carlson's data
        parameters, data_timesteps = parameters or {}, data_timesteps or {}
        self.parameters["data_timestep_hr"] = np.mean(np.diff(np.array(list(
            self.times.values())).flatten()))/hour if not hasattr(self, "data_timestep_hr") else self.data_timestep_hr
        self.parameters.update({
            "timestep_hr": self.parameters['data_timestep_hr'],
            "cvct": 0.01, "cvcf": 0.01,
            "bcv": 0.01, "cvmin": 0.01,
            "kcat": 0.33,
            'diffpos': 1, 'diffneg': 1,  # coefficients that weight difference between experimental and predicted biomass
            "stationary": 10,  # the penalty coefficient for the stationary phenotype
        })
        self.parameters.update(parameters)
        # distribute kcat values to all phenotypes of all species and update from previous simulations where necessary
        self.parameters.update(self._universalize(self.parameters, "kcat", exclude=["stationary"]))
        if primal_values is not None:
            for species, content in self.parameters["kcat"].items():
                if species not in primal_values:  continue
                for pheno, content2 in content.items():
                    if pheno not in primal_values[species]:  continue
                    for time, val in content2.items():
                        if time not in primal_values[species][pheno]:  continue
                        self.parameters["kcat"][species][pheno][time] = val
        print(self.parameters["kcat"])
        # define the metabolites that are tracked, exchanged, and not available in the media
        # TODO the default zero_start logic appears to be incorrect
        self.zero_start = zero_start or [met for met in self.consumed_mets
                                         if (met not in self.carbon_conc or self.carbon_conc[met] == 0)]
        self.rel_final_conc = rel_final_conc or {
            met:0.1 for met, concs in self.carbon_conc.items() if any(
                [concs[short_code] > 0 for short_code in self.data_df.index.unique()]
            ) and met not in self.zero_start}
        self.abs_final_conc = abs_final_conc or {}
        if mets_to_track:  self.mets_to_track = mets_to_track
        elif not isinstance(rel_final_conc, dict):  self.mets_to_track = self.fluxes_tup.index
        else:  self.mets_to_track = list(self.rel_final_conc.keys()) + self.zero_start
        print(self.mets_to_track)

        ts_to_delete = {}  # {short_code: full_times for short_code in unique_short_codes}
        if data_timesteps:  # {short_code:[times]}
            for short_code, times in data_timesteps.items():
                ts_to_delete[short_code] = set(list(range(len(full_times)))) - set(times)
                self.times[short_code] = np.delete(self.times[short_code], list(ts_to_delete[short_code]))

        # construct the problem
        objective = tupObjective("minimize variance and phenotypic transitions", [], "min")
        constraints, variables, simulated_mets = [], [], []
        time_1 = process_time()
        for exID in self.fluxes_tup.index:
            if exID == "bio":  continue
            met_id = re.search(r"(cpd\d{5})", exID).group()
            met = self.msdb.compounds.get_by_id(met_id)
            if "C" not in met.elements:  continue
            concID = f"c_{met_id}_e0"
            simulated_mets.append(met_id)
            self.variables[concID] = {}; self.constraints['dcc_' + met_id] = {}

            # define the growth rate for each metabolite and concentrations
            # TODO the MM parameters may be deletable once the binned kcat method is refined
            if "Vmax" and "Km" in self.parameters:
                self.parameters["Vmax"].update(self._universalize(self.parameters["Vmax"], met_id))
                self.parameters["Km"].update(self._universalize(self.parameters["Km"], met_id))
            for short_code in unique_short_codes:
                self.variables[concID][short_code] = {}; self.constraints['dcc_' + met_id][short_code] = {}
                timesteps = list(range(1, len(self.times[short_code]) + 1))
                for timestep in timesteps:
                    ## define the concentration variables
                    conc_var = tupVariable(_name(concID, "", short_code, timestep, self.names))
                    ## constrain initial time concentrations to the media or a large default
                    if timestep == timesteps[0]:
                        initial_val = None
                        if met_id in self.media_conc:  initial_val = self.media_conc[met_id]
                        if met_id in self.zero_start:  initial_val = 0
                        if dict_keys_exists(self.carbon_conc, met_id, short_code):
                            initial_val = self.carbon_conc[met_id][short_code]
                        if initial_val is not None:
                            conc_var = conc_var._replace(bounds=Bounds(initial_val, initial_val))
                            if biolog_simulation:  conc_var = conc_var._replace(bounds=Bounds(1, None))
                    ## mandate complete carbon consumption
                    elif timestep == timesteps[-1] and (met_id in self.rel_final_conc or met_id in self.abs_final_conc):
                        if met_id in self.rel_final_conc:
                            final_bound = self.variables[concID][short_code][1].bounds.lb * self.rel_final_conc[met_id]
                        if met_id in self.abs_final_conc:  # this intentionally overwrites rel_final_conc
                            final_bound = self.abs_final_conc[met_id]
                        conc_var = conc_var._replace(bounds=Bounds(0, final_bound))
                        if met_id in self.zero_start:
                            conc_var = conc_var._replace(bounds=Bounds(final_bound, final_bound))
                    self.variables[concID][short_code][timestep] = conc_var
                    variables.append(self.variables[concID][short_code][timestep])
        for pheno in self.phenotypes:
            self.constraints['dbc_' + pheno] = {short_code: {} for short_code in unique_short_codes}

        # define growth and biomass variables and constraints
        for pheno in self.phenotypes:
            for short_code in unique_short_codes:
                self.initialize_vars_cons(pheno, short_code)
                timesteps = list(range(1, len(self.times[short_code]) + 1))
                nth_percentile_timestep = timesteps[int(0.90*len(timesteps))]
                penalty_range = np.linspace(self.parameters['stationary'], self.parameters['stationary']/10,
                                            len(timesteps[nth_percentile_timestep:]))
                timestep_excess_count = 0
                for timestep in map(int, timesteps):
                    variables = self.define_b_vars(pheno, short_code, timestep, variables)
                    if short_code not in self.constraints[f"binc_{pheno}"]:
                        self.constraints[f"binc_{pheno}"][short_code] = tupConstraint(
                            _name("binc_", pheno, short_code, "", self.names), Bounds(0, 4), {
                                "elements": [self.variables[f"bin1_{pheno}"][short_code].name,
                                             self.variables[f"bin2_{pheno}"][short_code].name,
                                             self.variables[f"bin3_{pheno}"][short_code].name,
                                             self.variables[f"bin4_{pheno}"][short_code].name,
                                             self.variables[f"bin5_{pheno}"][short_code].name],
                                "operation": "Add"})
                        constraints.append(self.constraints[f'binc_{pheno}'][short_code])
                    constraints.extend(self.define_b_cons(pheno, short_code, timestep, biomass_coefs))

                    ## define the growth rate variable or primal value
                    species, phenotype = pheno.split("_")
                    self.variables['g_' + pheno][short_code][timestep] = tupVariable(
                        _name("g_", pheno, short_code, timestep, self.names))
                    variables.append(self.variables['g_' + pheno][short_code][timestep])

                    if 'stationary' in pheno:
                        weight = self.parameters['stationary']
                        if timestep > nth_percentile_timestep:
                            weight = penalty_range[timestep_excess_count]
                            timestep_excess_count += 1
                        objective.expr.extend([{
                            "elements": [{"elements": [weight, self.variables['b_' + pheno][short_code][timestep].name],
                                          "operation": "Mul"}],
                            "operation": "Add"}])
                        continue
                    # the conversion rates to and from the stationary phase
                    self.variables['cvt_' + pheno][short_code][timestep] = tupVariable(
                        _name("cvt_", pheno, short_code, timestep, self.names), Bounds(0, 100))
                    self.variables['cvf_' + pheno][short_code][timestep] = tupVariable(
                        _name("cvf_", pheno, short_code, timestep, self.names), Bounds(0, 100))
                    variables.extend([self.variables['cvf_' + pheno][short_code][timestep],
                                      self.variables['cvt_' + pheno][short_code][timestep]])

                    # cvt <= bcv*b_{pheno} + cvmin
                    self.constraints['cvc_' + pheno][short_code][timestep] = tupConstraint(
                        _name('cvc_', pheno, short_code, timestep, self.names), (0, None), {
                            "elements": [{"elements": [-1, self.variables['cvt_' + pheno][short_code][timestep].name],
                                          "operation": "Mul"}],
                            "operation": "Add"})
                    # biomass_term = [self.parameters['bcv']*b_value + self.parameters['cvmin']] if FBAHelper.isnumber(b_value) else [
                    biomass_term = [self.parameters['cvmin'],
                                    {"elements": [self.parameters['bcv'],
                                                  self.variables["b_"+pheno][short_code][timestep].name],
                                     "operation": "Mul"}]
                    self.constraints['cvc_' + pheno][short_code][timestep].expr["elements"].extend(biomass_term)

                    # g_{pheno} = b_{pheno}*v_{pheno}
                    b_values = [self.variables['b1_' + pheno][short_code][timestep].name,
                                self.variables['b2_' + pheno][short_code][timestep].name,
                                self.variables['b3_' + pheno][short_code][timestep].name,
                                self.variables['b4_' + pheno][short_code][timestep].name,
                                self.variables['b5_' + pheno][short_code][timestep].name]
                    self.constraints['gc_' + pheno][short_code][timestep] = tupConstraint(
                        name=_name('gc_', pheno, short_code, timestep, self.names),
                        expr={"elements": [*[{"elements": [-self.parameters["kcat"][species][phenotype], b],
                                              "operation": "Mul"} for b in b_values],
                                           self.variables['g_' + pheno][short_code][timestep].name],
                              "operation": "Add"})

                    constraints.extend([self.constraints['cvc_' + pheno][short_code][timestep],
                                        self.constraints['gc_' + pheno][short_code][timestep]])
                                        # self.constraints["binTot_" + pheno][short_code]])

        # define the concentration constraint
        half_dt = self.parameters['data_timestep_hr'] / 2
        time_2 = process_time()
        print(f'Done with concentrations and biomass loops: {(time_2 - time_1) / 60} min')
        for r_index, met in enumerate(self.fluxes_tup.index):
            met_id = _met_id_parser(met)
            if met_id not in simulated_mets:  continue
            concID = f"c_{met_id}_e0"
            for short_code in unique_short_codes:
                timesteps = list(range(1, len(self.times[short_code]) + 1))
                for timestep in timesteps[:-1]:
                    # c_{met} + dt/2*sum_k^K(n_{k,met} * (g_{pheno}+g+1_{pheno})) = c+1_{met}
                    next_timestep = timestep + 1
                    growth_phenos = [[self.variables['g_' + pheno][short_code][next_timestep].name,
                                      self.variables['g_' + pheno][short_code][timestep].name]
                                     for pheno in self.fluxes_tup.columns]
                    self.constraints['dcc_' + met_id][short_code][timestep] = tupConstraint(
                        name=_name("dcc_", met_id, short_code, timestep, self.names),
                        expr={
                            "elements": [
                                self.variables[concID][short_code][timestep].name,
                                {"elements": [-1, self.variables[concID][short_code][next_timestep].name],
                                 "operation": "Mul"},
                                *OptlangHelper.dot_product(
                                    growth_phenos, heuns_coefs=half_dt * self.fluxes_tup.values[r_index])],
                            "operation": "Add"})
                    constraints.append(self.constraints['dcc_' + met_id][short_code][timestep])

        #   define the conversion variables of every signal for every phenotype
        # for signal in growth_tup.columns[2:]:
        #     for pheno in self.fluxes_tup.columns:
        #         conversion_name = "_".join([signal, pheno, "__conversion"])
        #         self.variables[conversion_name] = tupVariable(conversion_name)
        #         variables.append(self.variables[conversion_name])

        time_3 = process_time()
        print(f'Done with DCC loop: {(time_3 - time_2) / 60} min')
        species_phenos = {}
        self.conversion_bounds = [5e-6, 50]
        for index, org_signal in enumerate(growth_tup.columns[2:]):
            # signal = org_signal.split(":")[1]
            signal = org_signal.replace(":", "|")
            species = signal_species(org_signal)
            species_phenos[species] = {None if "OD" in species else f"{species}_stationary"}
            signal_column_index = index + 2
            data_timestep = 1
            self.variables[signal + '|conversion'] = tupVariable(
                signal + '|conversion', bounds=Bounds(*self.conversion_bounds))
            variables.append(self.variables[signal + '|conversion'])

            self.variables[signal + '|bio'] = {}; self.variables[signal + '|diffpos'] = {}
            self.variables[signal + '|diffneg'] = {}; self.variables['g_' + species] = {}
            self.constraints[signal + '|bioc'] = {}; self.constraints[signal + '|diffc'] = {}
            self.constraints["gc_" + species] = {}; self.constraints["totVc_" + species] = {}
            self.constraints["totGc_" + species] = {}; self.constraints[signal + '|bio_finalc'] = {}
            for short_code in unique_short_codes:
                self.variables[signal + '|bio'][short_code] = {}
                self.variables[signal + '|diffpos'][short_code] = {}
                self.variables[signal + '|diffneg'][short_code] = {}
                self.variables['g_' + species][short_code] = {}
                self.constraints[signal + '|bioc'][short_code] = {}
                self.constraints[signal + '|diffc'][short_code] = {}
                self.constraints["gc_" + species][short_code] = {}
                self.constraints["totVc_" + species][short_code] = {}
                self.constraints["totGc_" + species][short_code] = {}
                # self.constraints[signal + '|bio_finalc'][short_code] = {}
                # the value entries are matched to only the timesteps that are condoned by data_timesteps
                values_slice = trial_contents(short_code, growth_tup.index, growth_tup.values)
                if ts_to_delete:  values_slice = np.delete(values_slice, list(ts_to_delete[short_code]), axis=0)
                timesteps = list(range(1, len(values_slice) + 1))
                # the last timestep is omitted since Heun's method in the modelled biomass
                ## requires a future timestep, which does not exist for the last timestep
                for timestep in timesteps[:-1]:
                    ## the user timestep and data timestep must be synchronized
                    if (int(timestep)*self.parameters['timestep_hr']
                            < data_timestep*self.parameters['data_timestep_hr']):
                        print(f"Skipping timestep {timestep} that does not align with the user's timestep") ; continue
                    data_timestep += 1
                    if data_timestep > int(self.times[short_code][-1] / self.parameters["data_timestep_hr"]):
                        print(f"The user-defined time exceeds the simulation time, so the DBC & diff loop is broken.")
                        break
                    next_timestep = int(timestep) + 1
                    ## the phenotype transition terms are aggregated
                    total_biomass, signal_sum, from_sum, to_sum = [], [], [], []
                    for pheno_index, pheno in enumerate(self.phenotypes):
                        ### define the collections of signal and pheno terms
                        if species in pheno or "OD" in signal:
                            # if not FBAHelper.isnumber(b_values[pheno][short_code][timestep]):
                            signal_sum.append({"operation": "Mul", "elements": [
                                -1, self.variables['b_' + pheno][short_code][timestep].name]})
                            # else:
                            #     signal_sum.append(-b_values[pheno][short_code][timestep])
                            ### total_biomass.append(self.variables["b_"+pheno][short_code][timestep].name)
                            if all(['OD' not in signal, species in pheno, 'stationary' not in pheno]):
                                species_phenos[species].add(pheno)
                                from_sum.append({"operation": "Mul", "elements": [
                                    -1, self.variables["cvf_" + pheno][short_code][timestep].name]})
                                to_sum.append(self.variables["cvt_" + pheno][short_code][timestep].name)
                    for pheno in species_phenos[species]:
                        if "OD" in signal:  continue
                        # print(pheno, timestep, b_values[pheno][short_code][timestep], b_values[pheno][short_code][next_timestep])
                        if "stationary" in pheno:
                            # b_{phenotype} - sum_k^K(es_k*cvf) + sum_k^K(pheno_bool*cvt) = b+1_{phenotype}
                            self.constraints['dbc_' + pheno][short_code][timestep] = tupConstraint(
                                name=_name("dbc_", pheno, short_code, timestep, self.names),
                                expr={"elements": [*from_sum, *to_sum], "operation": "Add"})
                        else:
                            # b_{phenotype} + dt/2*(g_{phenotype} + g+1_{phenotype}) + cvf-cvt = b+1_{phenotype}
                            self.constraints['dbc_' + pheno][short_code][timestep] = tupConstraint(
                                name=_name("dbc_", pheno, short_code, timestep, self.names),
                                expr={
                                    "elements": [
                                        self.variables['cvf_' + pheno][short_code][timestep].name,
                                        {"elements": [half_dt, self.variables['g_' + pheno][short_code][timestep].name],
                                         "operation": "Mul"},
                                        {"elements": [half_dt, self.variables['g_' + pheno][short_code][next_timestep].name],
                                         "operation": "Mul"},
                                        {"elements": [-1, self.variables['cvt_' + pheno][short_code][timestep].name],
                                         "operation": "Mul"}],
                                    "operation": "Add"})
                        # if not FBAHelper.isnumber(self.variables['b_' + pheno][short_code][timestep]):
                        biomass_term = [self.variables['b_' + pheno][short_code][timestep].name, {
                            "elements": [-1, self.variables['b_' + pheno][short_code][next_timestep].name],
                            "operation": "Mul"}]
                        # else:
                        #     biomass_term = [b_values[pheno][short_code][timestep]-b_values[pheno][short_code][next_timestep]]
                        self.constraints['dbc_' + pheno][short_code][timestep].expr["elements"].extend(biomass_term)
                        constraints.append(self.constraints['dbc_' + pheno][short_code][timestep])

                    if not requisite_biomass or any([timestep != timesteps[-2], signal not in requisite_biomass[short_code]]):
                        self.variables[signal + '|bio'][short_code][timestep] = tupVariable(
                            _name(signal, '|bio', short_code, timestep, self.names))
                    else:
                        biomass_flux = requisite_biomass[short_code][signal]["bio"]
                        estimated_biomass = biomass_flux #* int(timestep)*self.parameters['data_timestep_hr']
                        self.variables[signal + '|bio'][short_code][timestep] = tupVariable(
                            _name(signal, '|bio', short_code, timestep, self.names),
                            Bounds(estimated_biomass, None))
                    self.variables[signal + '|diffpos'][short_code][timestep] = tupVariable(
                        _name(signal, '|diffpos', short_code, timestep, self.names), Bounds(0, 100))
                    self.variables[signal + '|diffneg'][short_code][timestep] = tupVariable(
                        _name(signal, '|diffneg', short_code, timestep, self.names), Bounds(0, 100))
                    variables.extend([self.variables[signal + '|bio'][short_code][timestep],
                                      self.variables[signal + '|diffpos'][short_code][timestep],
                                      self.variables[signal + '|diffneg'][short_code][timestep]])

                    # {signal}__conversion*datum = {signal}__bio
                    # TODO - the conversion variable must be a constant for BIOLOG conditions
                    self.constraints[signal + '|bioc'][short_code][timestep] = tupConstraint(
                        name=_name(signal, '|bioc', short_code, timestep, self.names),
                        expr={
                            "elements": [
                                {"elements": [-1, self.variables[signal + '|bio'][short_code][timestep].name],
                                 "operation": "Mul"},
                                {"elements": [self.variables[signal + '|conversion'].name,
                                              values_slice[timestep, signal_column_index]],
                                 "operation": "Mul"}],
                            "operation": "Add"})
                    constraints.append(self.constraints[signal + '|bioc'][short_code][timestep])

                    # {speces}_bio + {signal}_diffneg-{signal}_diffpos = sum_k^K(es_k*b_{phenotype})
                    self.constraints[signal + '|diffc'][short_code][timestep] = tupConstraint(
                        name=_name(signal, '|diffc', short_code, timestep, self.names),
                        expr={
                            "elements": [
                                self.variables[signal + '|bio'][short_code][timestep].name,
                                self.variables[signal + '|diffneg'][short_code][timestep].name,
                                {"elements": [-1, self.variables[signal + '|diffpos'][short_code][timestep].name],
                                 "operation": "Mul"}],
                            "operation": "Add"})
                    if all([isinstance(val, dict) for val in signal_sum]):
                        self.constraints[signal + "|diffc"][short_code][timestep].expr["elements"].extend(signal_sum)
                    else:  raise ValueError(f"The {signal_sum} value has unexpected contents.")
                    constraints.append(self.constraints[signal + '|diffc'][short_code][timestep])

                    objective.expr.extend([{
                        "elements": [
                            {"elements": [self.parameters['diffpos'],
                                          self.variables[f'{signal}|diffpos'][short_code][timestep].name],
                             "operation": "Mul"},
                            {"elements": [self.parameters['diffneg'],
                                          self.variables[f'{signal}|diffneg'][short_code][timestep].name],
                             "operation": "Mul"}],
                        "operation": "Add"}])

        time_4 = process_time()
        print(f'Done with the DBC & diffc loop: {(time_4 - time_3) / 60} min')

        # construct the problem
        self.problem = OptlangHelper.define_model("CommPhitting model", variables, constraints, objective, True)
        self.hdf5_name = export_lp.replace(".lp", ".h5")
        self.hdf5_file = File(self.hdf5_name, 'w')
        time_5 = process_time()
        print(f'Done with constructing the {type(self.problem)} model: {(time_5 - time_4) / 60} min')

        # export contents
        if export_phenotype_profiles:
            phenotype_profiles_name = 'phenotype_profiles.tsv'
            self.fluxes_df.to_csv(phenotype_profiles_name, sep="\t")
            self.zipped_output.append(phenotype_profiles_name)
        if export_parameters:
            parameter_name = 'parameters.tsv'
            DataFrame(data=list(self.parameters.values()), index=list(self.parameters.keys()),
                      columns=['values']).to_csv(parameter_name, sep="\t")
            self.zipped_output.append(parameter_name)
        if export_lp:
            if re.search(r"(\\\\/)", export_lp):  os.makedirs(os.path.dirname(export_lp), exist_ok=True)
            with open(export_lp, 'w') as lp:  lp.write(self.problem.to_lp())
            model_name = 'CommPhitting.json'
            _export_model_json(self.problem.to_json(), model_name)
            self.zipped_output.extend([export_lp, model_name])
        if export_zip_name:
            self.zip_name = export_zip_name
            sleep(2)
            with ZipFile(self.zip_name, 'a', compression=ZIP_LZMA) as zp:
                for file in self.zipped_output:
                    zp.write(file) ; os.remove(file) ; self.zipped_output.remove(file)
        time_6 = process_time()
        print(f'Done exporting the content: {(time_6 - time_5) / 60} min')

    def compute(self, graphs: list = None, export_zip_name=None, figures_zip_name=None, publishing=False,
                primals_export_path:str = "primal_values.json", remove_empty_plots=False):
        print("starting optimization")
        time1 = process_time()
        self.values = {}
        solution = self.problem.optimize()
        timesteps = min(list(map(len, self.times.values())))
        fit_quality = self.problem.objective.value/timesteps
        print(f"The optimization fit quality is {fit_quality}")
        if "parameters.tsv" in self.zipped_output:
            self.parameters["fit"] = fit_quality
            parameter_name = 'parameters.tsv'
            DataFrame(data=list(self.parameters.values()), index=list(self.parameters.keys()),
                      columns=['values']).to_csv(parameter_name, sep="\t")
            with ZipFile(self.zip_name, 'a', compression=ZIP_LZMA) as zp:
                for file in self.zipped_output:
                    zp.write(file) ; os.remove(file)

        # TODO approximate a threshold of good fits, and trigger black box optimization for bad fits
        ## that iteratively adjust parameters until the fit metric surmounts the threshold.

        # categorize the primal values by trial and time
        if "optimal" not in solution:
            raise FeasibilityError(f'The solution is sub-optimal, with a(n) {solution} status.')
        if all(np.array(list(self.problem.primal_values.values())) == 0):
            raise NoFluxError("The simulation lacks any flux.")
        for variable, value in self.problem.primal_values.items():
            if "v_" in variable:  self.values[variable] = value
            elif 'conversion' in variable or re.search(r"(bin\d)", variable):
                self.values[short_code].update({variable: value})
                if value in self.conversion_bounds:
                    warnings.warn(f"The conversion factor {value} optimized to a bound, which may be "
                                  f"indicative of an error, such as improper kinetic rates.")
            else:
                basename, short_code, timestep = variable.split('-')
                time_hr = int(timestep) * self.parameters['data_timestep_hr']
                self.values[short_code] = self.values.get(short_code, {})
                self.values[short_code][basename] = self.values[short_code].get(basename, {})
                self.values[short_code][basename][time_hr] = value

        # export the processed primal values for graphing
        # with open(primals_export_path, 'w') as out:
        #     json.dump(self.values, out, indent=3)
        # if not export_zip_name and hasattr(self, 'zip_name'):
        #     export_zip_name = self.zip_name
        # if export_zip_name:
        #     with ZipFile(export_zip_name, 'a', compression=ZIP_LZMA) as zp:
        #         zp.write(primals_export_path)
        #         os.remove(primals_export_path)
        # visualize the specified information
        time2 = process_time()
        if graphs:  self.graph(graphs, export_zip_name=figures_zip_name or export_zip_name,
                               publishing=publishing, remove_empty_plots=remove_empty_plots)

        # parse the primal values
        values_df = DataFrame(self.values)
        values_index = values_df.index.tolist()
        for col in values_df.columns:
            trial_values = values_df[col].tolist()
            ## process the times
            times = [list(ele.keys()) for ele in trial_values if isinstance(ele, dict)]
            max_time = max(list(map(len, times)))
            for max_time_series in times:
                if len(max_time_series) == max_time:  break
            trial_path = f'results/primals/{col}/'
            self.hdf5_file.create_dataset(f'{trial_path}/times', data=max_time_series)
            ## process the data values
            for index, ele in enumerate(trial_values):
                dataset_name = f'{trial_path}/{values_index[index]}'
                if FBAHelper.isnumber(ele):  self.hdf5_file.create_dataset(dataset_name, data=[float(ele)])
                elif isinstance(ele, dict):
                    self.hdf5_file.create_dataset(dataset_name, data=list(map(float, ele.values())))
                    self.hdf5_file[dataset_name].attrs["full_time"] = (len(ele.values()) == max_time)

        self.hdf5_file.close()
        with ZipFile(self.zip_name, 'a', compression=ZIP_LZMA) as zp:
            zp.write(self.hdf5_name)  ;  os.remove(self.hdf5_name)

        time3 = process_time()
        print(f"Optimization completed in {(time2-time1)/60} minutes")
        print(f"Graphing completed in {(time3-time2)/60} minutes")

    def load_model(self, mscomfit_json_path: str = None, zip_name: str = None, model_to_load: dict = None):
        if zip_name:
            with ZipFile(zip_name, 'r') as zp:  zp.extract(mscomfit_json_path)
        if mscomfit_json_path:
            with open(mscomfit_json_path, 'r') as mscmft:  return json.load(mscmft)
        if model_to_load:  self.problem = Model.from_json(model_to_load)

    @staticmethod
    def assign_values(param, var, next_dimension, kcat=True):
        dic = {var: {}}
        for dim1, dim2_list in next_dimension.items():
            if isinstance(dim2_list, dict):  dic[var].update(CommPhitting.assign_values(param, dim1, dim2_list))
            else:
                if kcat:  dic[var][dim1] = param
                else:  dic[var][dim1] = {dim2: param for dim2 in dim2_list}
        return dic

    def _universalize(self, param, var, next_dimension=None, exclude=None, tsBin=False):
        if not next_dimension:
            next_dimension = {}
            for organism in self.fluxes_tup.columns:
                species, pheno = organism.split("_")
                if pheno in exclude:  continue
                if not tsBin:
                    if species in next_dimension:  next_dimension[species].append(pheno)
                    else:  next_dimension[species] = [pheno]
                else:
                    if species in next_dimension:  next_dimension[species].update({pheno: self.time_ranges})
                    else:  next_dimension[species] = {pheno: self.time_ranges}
        if FBAHelper.isnumber(param):  return CommPhitting.assign_values(param, var, next_dimension)
        elif FBAHelper.isnumber(param[var]):  return CommPhitting.assign_values(param[var], var, next_dimension)
        elif isinstance(param[var], dict):
            return {var: {dim1: {dim2: param[var][dim1] for dim2 in dim2_list}
                          for dim1, dim2_list in next_dimension.items()}}
        else:  logger.critical(f"The param (with keys {dic_keys(param)}) and var {var} are not amenable"
                               " with the parameterizing a universal value.")
                    # {short_code: {list(timestep_info.keys())[0]: find_dic_number(param)} for short_code, timestep_info in variable.items()}}

    def adjust_color(self, color, amount=0.5):
        """
        Lightens the given color by multiplying (1-luminosity) by the given amount.
        Input can be matplotlib color string, hex string, or RGB tuple.

        Examples:
        >> lighten_color('g', 0.3)
        >> lighten_color('#F034A3', 0.6)
        >> lighten_color((.3,.55,.1), 0.5)
        """
        import colorsys
        import matplotlib.colors as mc
        try:  c = mc.cnames[color]
        except:  c = color
        c = colorsys.rgb_to_hls(*mc.to_rgb(c))
        return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

    def _add_plot(self, ax, labels, label, basename, trial, x_axis_split, linestyle="solid",
                  scatter=False, color=None, xs=None, ys=None):
        labels.append(label or basename.split('-')[-1])
        xs = xs if xs is not None else list(map(float, self.values[trial][basename].keys()))
        ys = ys if ys is not None else list(map(float, self.values[trial][basename].values()))
        if scatter:  ax.scatter(xs, ys, s=10, label=labels[-1], color=color or None)
        else:  ax.plot(xs, ys, label=labels[-1], linestyle=linestyle, color=color or None)
        ax.set_xticks(list(map(int, xs))[::x_axis_split])
        return ax, labels

    def graph(self, graphs, primal_values_filename: str = None, primal_values_zip_path: str = None,
              export_zip_name: str = None, data_timestep_hr: float = 0.163, publishing: bool = False,
              title: str = None, remove_empty_plots:bool = False):
        print(export_zip_name)
        # define the default timestep ratio as 1
        data_timestep_hr = self.parameters.get('data_timestep_hr', data_timestep_hr)
        timestep_ratio = data_timestep_hr / self.parameters.get('timestep_hr', data_timestep_hr)
        if primal_values_filename:
            if primal_values_zip_path:
                with ZipFile(primal_values_zip_path, 'r') as zp:  zp.extract(primal_values_filename)
            with open(primal_values_filename, 'r', encoding='utf-8') as primal:  self.values = json.load(primal)

        # plot the content for desired trials
        x_axis_split = int(3 / data_timestep_hr / timestep_ratio)
        self.plots = set()
        contents = {"biomass": 'b_', "all_biomass": 'b_', "growth": 'g_', "conc": "c_"}
        mM_threshold = 1e-3
        for graph_index, graph in enumerate(graphs):
            content = contents.get(graph['content'], graph['content'])
            y_label = 'Variable value'; x_label = r'Time ($hr$)'
            if any([x in graph['content'] for x in ['biomass', 'OD']]):
                total_biomasses = {name: [] for name in self.species_list}
                total_biomasses.update({"OD":[]})
                if "species" not in graph:  graph['species'] = self.species_list
            if "biomass" in graph['content']:  y_label = r'Biomass ($\frac{g}{L}$)'
            elif 'growth' in graph['content']:  y_label = r'Biomass growth ($\frac{g}{hr}$)'
            graph["experimental_data"] = graph.get("experimental_data", False)
            if "painting" not in graph:
                graph["painting"] = {
                    "OD": {
                        "color": "blue",
                        "linestyle": "solid",
                        "name": "Total biomass"
                    },
                    "ecoli": {
                        "color": "red",
                        "linestyle": "dashed",
                        "name": "E. coli"
                    },
                    "pf": {
                        "color": "green",
                        "linestyle": "dotted",
                        "name": "P. fluorescens"
                    }}
            graph["parsed"] = graph.get("parsed", False)
            if 'phenotype' in graph and graph['phenotype'] == '*':
                if "species" not in graph:  graph['species'] = self.species_list
                graph['phenotype'] = set([pheno.split("_")[-1] for pheno in self.phenotypes
                                          if pheno.split("_")[0] in graph["species"]])
            # TODO - a species-resolved option must be developed for the paper figure
            if 'species' in graph and graph['species'] == '*':  graph['species'] = self.species_list
            elif content == "c_" and 'mets' not in graph:
                print(self.mets_to_track)
                graph["mets"] = self.mets_to_track
            elif not any(["species" in graph, "mets" in graph]):
                raise ValueError(f"The specified graph {graph} must define species for which data will be plotted.")
            print(f"graph_{graph_index}") ; pprint(graph)

            # define figure specifications
            if publishing:
                pyplot.rc('axes', titlesize=22, labelsize=28)
                pyplot.rc('xtick', labelsize=24)
                pyplot.rc('ytick', labelsize=24)
                pyplot.rc('legend', fontsize=18)
            if graph["parsed"]:
                parsed_graphs = {}
                for species in graph["species"]:
                    parsed_graphs[species] = pyplot.subplots(dpi=200, figsize=(11, 7))
            else:  fig, ax = pyplot.subplots(dpi=200, figsize=(11, 7))
            yscale = "linear"

            # populate the figures
            for trial, basenames in self.values.items():
                if trial not in graph['trial']:  continue
                labels = []
                for basename, values in basenames.items():
                    # graph experimental and total simulated biomasses
                    if any([x in graph['content'] for x in ['biomass', 'OD']]):
                        if 'b_' in basename:
                            vals = list(map(float, values.values()))
                            var_name, species, phenotype = basename.split('_')
                            # ic(basename)
                            label = f'{species}_biomass (model)'
                            if publishing:
                                species_name = graph["painting"][species]["name"]
                                label = f'{species_name} total (model)'
                            labels.append({species: label})
                            if remove_empty_plots and all([v == 0 for v in vals]):
                                print(f"The {basename} is empty and thus is removed.")
                                continue
                            if (any([x in graph['content'] for x in ["total", "biomass", 'OD']]) or
                                    graph['species'] == self.species_list): # and not graph["parsed"]:
                                total_biomasses['OD'].append(vals)
                                if "OD" not in graph['content']:  total_biomasses[species].append(vals)
                        if all([graph['experimental_data'], '|bio' in basename, ]):
                            # any([content in basename])]):  # TODO - any() must include all_biomass and total
                            species, signal, phenotype = basename.split('|')
                            label = basename
                            if publishing:
                                species_name = "total" if "OD" in signal else graph["painting"][species]["name"]
                                label = f'Experimental {species_name} (from {signal})'
                            # print(basename, label, self.values[trial][basename].values())
                            if remove_empty_plots and all(self.values[trial][basename].values() == 0):
                                print(f"The {basename} is empty and thus is removed.")
                                continue
                            ax, labels = self._add_plot(ax, labels, label, basename, trial, x_axis_split, scatter=True,
                                                        color=self.adjust_color(graph["painting"][species]["color"], 1.5))

                    if content not in basename:  continue
                    # graph individual phenotypes
                    if "phenotype" in graph:
                        # print(graph['phenotype'])
                        for specie in graph["species"]:
                            if specie not in basename:  continue
                            if not any([p in basename for p in graph['phenotype']]):
                                print(f"{basename} data with unknown phenotype.")
                                continue
                            if remove_empty_plots and all(self.values[trial][basename].values() == 0):
                                print(f"The {specie} is empty and thus is removed.")
                                continue
                            if graph["parsed"]:  fig, ax = parsed_graphs[specie]
                            ## define graph characteristics
                            label = basename.split("_")[-1]
                            style = "solid"
                            if len(graph["species"]) > 1:
                                label = re.sub(r"(^[a-b]+\_)", "", basename)
                                style = graph["painting"][specie]["linestyle"]
                            ax, labels = self._add_plot(ax, labels, label, basename, trial, x_axis_split, style)
                            if graph["parsed"]:  parsed_graphs[specie] = (fig, ax)
                    # graph media concentration plots
                    elif "mets" in graph and all([any([x in basename for x in graph["mets"]]), 'c_cpd' in basename]):
                        if not any(np.array(list(self.values[trial][basename].values())) > mM_threshold):  continue
                        if remove_empty_plots and all(self.values[trial][basename].values() == 0):  continue
                        label=self.msdb.compounds.get_by_id(re.search(r"(cpd\d+)", basename).group()).name
                        ax, labels = self._add_plot(ax, labels, label, basename, trial, x_axis_split)
                        yscale = "log"
                        y_label = r'Concentration ($mM$)'

                if labels:  # assesses whether graph(s) were created
                    ## graph all of the total biomasses
                    if any([x in graph['content'] for x in ['OD', 'biomass', 'total']]):
                        labeled_species = [label for label in labels if isinstance(label, dict)]
                        for name, vals in total_biomasses.items():
                            # ic(name)
                            if not vals or (len(total_biomasses) == 2 and "OD" not in name):  continue
                            if len(total_biomasses) == 2:
                                specie_label = [graph["painting"][name]["name"] for name in total_biomasses
                                                if "OD" not in name][0]
                                label = f"{graph['painting'][name]['name']} ({specie_label})"
                            else:
                                label = f'{name}_biomass (model)'
                                if labeled_species:
                                    for label_specie in labeled_species:
                                        if name in label_specie:  label = label_specie[name] ; break
                            style = "solid" if (len(graph["species"]) < 1 or name not in graph["painting"]
                                                ) else graph["painting"][name]["linestyle"]
                            style = "dashdot" if "model" in label else style
                            style = "solid" if ("OD" in name and not graph["experimental_data"]
                                                or "total" in graph["content"]) else style
                            total_biomass = sum(np.array(vals))[:-1]
                            xs = list(map(float, values.keys()))
                            if graph["parsed"]:  fig, ax = parsed_graphs[name]
                            self._add_plot(ax, labels, label, None, None, x_axis_split, style, False,
                                           graph["painting"][name]["color"], xs, total_biomass)
                            if graph["parsed"]:
                                ## process and export the parsed figures
                                ax.set_xlabel(x_label) ; ax.set_ylabel(y_label) ; ax.grid(axis="y")
                                ax.set_yscale(yscale) ; ax.legend()
                                phenotype_id = graph.get('phenotype', "")
                                if "phenotype" in graph and not isinstance(graph['phenotype'], str):
                                    phenotype_id = f"{','.join(graph['phenotype'])} phenotypes"
                                fig_name = f'{"_".join([trial, name, phenotype_id, content])}.jpg'
                                fig.savefig(fig_name, bbox_inches="tight", transparent=True)
                                self.plots.add(fig_name)

                    if graph["parsed"]:  continue
                    ## process and export the non-parsed figures
                    phenotype_id = graph.get('phenotype', "")
                    if "phenotype" in graph and not isinstance(graph['phenotype'], str):
                        phenotype_id = f"{','.join(graph['phenotype'])} phenotypes"

                    species_id = ""
                    if "mets" not in graph and content != "c_":
                        species_id = graph["species"] if isinstance(graph["species"], str) else ",".join(graph["species"])
                        if "species" in graph and graph['species'] == self.species_list:  species_id = 'all species'
                        else:  phenotype_id = f"{','.join(graph['species'])} species"
                        if species_id == "all species" and not phenotype_id: phenotype_id = ','.join(graph['species'])

                    ax.set_xlabel(x_label) ; ax.set_ylabel(y_label)
                    if "mets" in graph:  ax.set_ylim(mM_threshold)
                    ax.grid(axis="y")
                    if len(labels) > 1:  ax.legend()
                    else:  yscale = "linear"
                    ax.set_yscale(yscale)
                    if not publishing:
                        if not title:
                            org_content = content if content not in contents.values() else list(
                                contents.keys())[list(contents.values()).index(content)]
                            this_title = f'{org_content} of {species_id} ({phenotype_id}) in the {trial} trial'
                            if content == "c_":  this_title = f"{org_content} in the {trial} trial"
                            ax.set_title(this_title)
                        else:  ax.set_title(title)
                    fig_name = f'{"_".join([trial, species_id, phenotype_id, content])}.jpg'
                    if "mets" in graph:  fig_name = f"{trial}_{','.join(graph['mets'])}_c.jpg"
                    fig.savefig(fig_name, bbox_inches="tight", transparent=True)

                    self.plots.add(fig_name)

        # export the figures with other simulation content
        if export_zip_name:
            with ZipFile(export_zip_name, 'a', compression=ZIP_LZMA) as zp:
                for plot in self.plots:
                    zp.write(plot)  ;  os.remove(plot)


    #################### ENGINEERING PHASE METHODS ####################

    def engineering(self):
        if not hasattr(self, "problem"):
            self.fit()  # TODO - accommodate both fitting a new model and loading an existing model

        # This will capture biomass variables at all times and trials, which seems undesirable
        self.problem.objective = Objective(sum([x for x in self.problem.variables if "bio" in x.name]))

        # Use a community COBRA model and CommKinetics with the fitted kinetic parameters?

    def _add_phenotypes(self):
        pass



    def _change_obj(self):
        pass


class BIOLOGPhitting(CommPhitting):
    def __init__(self, carbon_conc, media_conc, biolog_df, fluxes_df,
                 experimental_metadata, msdb_path, community_members):
        self.biolog_df = biolog_df; self.experimental_metadata = experimental_metadata
        self.carbon_conc = carbon_conc; self.media_conc = media_conc or []
        self.fluxes_df = fluxes_df ; self.phenotypes = list(self.fluxes_df.columns)
        self.phenotypes.extend([signal_species(signal)+"_stationary"
                                for signal in self.biolog_df if ":" in signal])
        self.community_members = community_members
        # import os
        from modelseedpy.biochem import from_local
        self.msdb_path = msdb_path ; self.msdb = from_local(msdb_path)

    def fitAll(self, parameters: dict = None, rel_final_conc: float = None,
               abs_final_conc: dict = None, graphs: list = None, data_timesteps: dict = None,
               export_zip_name: str = None, export_parameters: bool = True, requisite_biomass: dict = None,
               figures_zip_name: str = None, publishing: bool = False):
        # simulate each condition
        if export_zip_name and os.path.exists(export_zip_name):
            os.remove(export_zip_name)
        org_rel_final_conc = rel_final_conc
        # total_reactions = set(list(chain.from_iterable([model.reactions for model in models_dict.values()])))
        model_abbreviations = ','.join([content["name"] for content in self.community_members.values()])
        for exp_index, experiment in self.experimental_metadata.iterrows():
            print(f"\n{exp_index} {experiment}")
            display(experiment)
            pheno = experiment["ModelSEED_ID"]
            if not pheno:
                print("The BIOLOG condition is not defined.")
                continue
            for model in self.community_members:
                cpd = self.msdb.compounds.get_by_id(pheno)
                if "C" not in cpd.elements or not any([re.search(pheno, rxn.id) for rxn in model.reactions]):
                    if "valid_condition" not in locals():
                        valid_condition = False
                    continue
                exp_list = [pheno] if isinstance(pheno, str) else pheno
                self.community_members[model].update({"phenotypes": {
                    re.sub(r"(-|\s)", "", experiment["condition"]): {"consumed": exp_list} }})
                # determine the requisite biomass for each condition based on which member consumes the compound
                valid_condition = True
            # proceed if none of the members can utilize the phenotype condition
            if not valid_condition:
                print(f"The BIOLOG condition with {experiment['ModelSEED_ID']} is not"
                      f" absorbed by the {model_abbreviations} model(s).")
                continue
            print(f"The {experiment['ModelSEED_ID']} ({cpd.formula}) metabolite of the "
                  f"{experiment['condition']} condition may feed the {model_abbreviations} model(s).")
            if not any([experiment["ModelSEED_ID"] in pheno for pheno in self.phenotypes]):
                print(e)
                print(f"The {experiment['ModelSEED_ID']} ({cpd.formula}) metabolite of the "
                      f"{experiment['condition']} condition is not a suitable phenotype for "
                      f"the {model_abbreviations} model(s).")
                continue

            # for exp_index, experiment in self.experimental_metadata.iterrows():
            # the model(s) for which the condition is a suitable carbon source must be defined here
            # simulate through the kinetics ranges with conditions that can be used by one of members
            rel_final_conc = {experiment["ModelSEED_ID"]: org_rel_final_conc}
            export_path = os.path.join(os.getcwd(), "BIOLOG_LPs", f"{exp_index}_{','.join(exp_list)}.lp")
            kcat_primal = None
            for coef_index, coefs in enumerate(biomass_partition_coefs):
                # solve for growth rate constants with the previously solved biomasses
                new_simulation = CommPhitting(self.fluxes_df, self.carbon_conc, self.media_conc,
                                              self.msdb_path, self.biolog_df.loc[exp_index,:],
                                              self.experimental_metadata)
                new_simulation.define_problem(
                    parameters, exp_list, rel_final_conc,
                    set(list(chain.from_iterable([
                        content["excretions"] for content in self.community_members.values()]))),
                    abs_final_conc, data_timesteps, export_zip_name, export_parameters, export_path,
                    kcat_primal, coefs, requisite_biomass, True)
                time1 = process_time()
                primals_export_path = primals_export_path or f"BIOLOG_{experiment['ModelSEED_ID']}.json"
                try:
                    new_simulation.compute(graphs, export_zip_name, None, publishing, primals_export_path, True)
                except (NoFluxError) as e:
                    print(e)
                kcat_primal = parse_primals(new_simulation.values, coefs=coefs,
                                            kcat_vals=new_simulation.parameters["kcat"])
                time2 = process_time()
                print(f"Done simulating with the coefficients for biomass partitions: {coef_index}"
                      f"\n{(time2 - time1) / 60} minutes")
                pprint(kcat_primal)
            print("\n\n\n")
        return {k: val for k, val in new_simulation.values.items() if "kcat" in k}
