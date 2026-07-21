# -*- coding: utf-8 -*-

from scipy.constants import milli, hour, minute, day, femto
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy import MSModelUtil
from optlang import Constraint
from modelseedpy.core.fbahelper import FBAHelper
from collections import OrderedDict
from optlang.symbolics import Zero
from numpy import log10, nan, mean
from warnings import warn
from matplotlib import pyplot
from pprint import pprint
from datetime import date
from math import inf
import pandas
import json, re, os


def _x_axis_determination(total_time):
    time = total_time * minute
    if time <= 600:
        return minute, "s"
    if time > 600:
        return 1, "min"
    if time > 7200:
        return 1 / hour, "hr"
    return 1 / day, "days"


def _check_datum(datum):
    if "substituted_rate_law" not in datum:
        print(f"RateLawError: The {datum} datum lacks a rate law.")
        return False
    remainder = re.sub("([0-9A-Za-z/()e\-\+\.\*\_])", "", datum["substituted_rate_law"])
    if remainder != "":
        print(
            f'RateLawError: The {datum["substituted_rate_law"]}'
            f" rate law contains unknown characters: {remainder}"
        )
        return False
    return True


class MSKineticsFBA:
    def __init__(
        self,
        model,
        warnings: bool = True,
        verbose: bool = False,
        printing: bool = False,
        jupyter: bool = False,
    ):
        self.warnings, self.verbose, self.printing, self.jupyter = (
            warnings,
            verbose,
            printing,
            jupyter,
        )
        self.model_util = MSModelUtil(model)
        self.met_ids = OrderedDict(
            {met.id: met.id for met in self.model_util.model.metabolites}
        )

    def baseKinFBA(
        self,
        kinetics_path: str = None,
        kinetics_data: dict = None,
        initial_M: dict = None,  # a dictionary of the initial metabolic concentrations, which supplants concentrations from the defined kinetics data
        total_min: float = 200,
        ts_min: float = 20,
        export_name=None,
        export_directory=None,
        chemostat_L: float = None,
        feed_profile: dict = None,
        chemostat_L_hr: float = None,
        temperature: float = 25,
        p_h: float = 7,
        cell_dry_g: float = 1.44e-13,
        cellular_L: float = 1e-18,
        conc_figure_title="Metabolic perturbation",
        included_mets: list = None,
        labeled_plots=True,
        visualize=True,
        export=True,
    ):
        # define the dataframe for the time series content
        feed_profile, constrained, self.constraints = feed_profile or {}, {}, {}
        included_mets, self.sols = included_mets or [], []
        self.parameters = {
            "timesteps": int(total_min / ts_min),
            "pH": p_h,
            "temperature": temperature,
        }
        self.variables = {"elapsed_time": 0}
        self.ts_min, self.minimum = ts_min, inf
        timestep_hr = self.ts_min / (hour / minute)
        self.constrained = OrderedDict()
        cell_g_L = (
            cell_dry_g / cellular_L
        )  # https://journals.asm.org/doi/full/10.1128/AEM.64.2.688-694.1998

        # define reaction kinetics and initial concentrations
        assert (
            kinetics_path or kinetics_data
        ), "Either < kinetics_path > or < kinetics_data > must be provided"
        if kinetics_path:
            with open(kinetics_path) as data:
                self.kinetics_data = json.load(data)
        elif kinetics_data:
            self.kinetics_data = kinetics_data.copy()
        ## define the concentration, moles, and fluxes DataFrames
        self.time = "0 min"
        self.conc = pandas.DataFrame(
            [0] * len(self.met_ids),
            index=list(self.met_ids.keys()),
            columns=[self.time],
        )
        self.conc.index.name = "metabolite (mM)"
        self.moles = self.conc.copy(deep=True)
        self.fluxes = pandas.DataFrame(
            index=[rxn.id for rxn in self.model_util.model.reactions],
            columns=[self.time],
        )
        self.fluxes.index.name = "reaction (\u0394mmol/hr*g_(dw)))"  # Delta
        ## parse the kinetics data
        for content in self.kinetics_data.values():
            for condition, datum in content.items():
                if "initial_M" not in datum:
                    continue
                for var, conc in datum["initial_M"].items():
                    met_id = datum["met_id"][var]
                    if met_id in self.met_ids:
                        self.conc.at[met_id, self.time] += conc / milli
                    elif self.warnings:
                        warn(
                            f"KineticsError: The {met_id} reagent ({var}) in the"
                            f" {datum['substituted_rate_law']} rate law is not defined by the model."
                        )
        ## incorporate custom initial concentrations, which overwrites values from the kinetics data
        for met_id in initial_M:
            self.conc.at[met_id, self.time] = initial_M[met_id] / milli
        defined_concs = self.conc[self.conc[self.time] != 0][self.time].to_dict()
        chemostat_requirements = [
            chemostat_L is not None,
            feed_profile != {},
            chemostat_L_hr is not None,
        ]
        # execute FBA for each timestep, then calculate custom fluxes, constrain the model, and update concentrations
        model_rxns = [rxn.id for rxn in self.model_util.model.reactions]
        newTime = 0
        for timestep in range(1, self.parameters["timesteps"] + 1):
            oldTime = newTime
            newTime = timestep * self.ts_min
            t = timestep * timestep_hr
            self.previous_time = f"{oldTime} min"
            self.time = f"{newTime} min"
            self.conc[self.time] = [float(0)] * len(self.conc.index)
            self.fluxes[self.time] = [0] * len(self.fluxes.index)
            ## create a metabolite variable that prevents negative concentrations
            for met in self.model_util.model.metabolites:
                if met.id not in defined_concs:
                    continue
                if met.id not in self.constraints:
                    self.constraints[met.id] = {}
                coef = {}
                for rxn in met.reactions:
                    ### The product of the reaction stoichiometry and the timestep
                    stoich = abs(timestep_hr * rxn.metabolites[met])
                    coef[rxn.forward_variable], coef[rxn.reverse_variable] = (
                        stoich,
                        -stoich,
                    )
                ### build the metabolite constraint
                if newTime - self.ts_min in self.constraints[met.id]:
                    self.model_util.remove_cons_vars(
                        [self.constraints[met.id][newTime - self.ts_min]]
                    )
                self.constraints[met.id][newTime] = Constraint(
                    Zero, lb=0, ub=None, name=f"{met.id}_conc"
                )
                self.model_util.create_constraint(
                    self.constraints[met.id][newTime], coef
                )
            ## calculate the flux
            display(self.conc[self.conc["0 min"] != 0], self.fluxes)
            for rxnID in self.kinetics_data:
                # TODO allocate the following code into a function and recusively reduce the timestep until
                ## the concentration becomes not negative, following the model of microBialSim. This may require
                ## time dependency in the kinetics expression to achieve the desired behavior.
                if rxnID not in model_rxns and self.warnings:
                    warn(f"ReactionError: {rxnID} is not in the model.")
                    continue
                fluxes = []
                for source in self.kinetics_data[rxnID]:
                    datum = self.kinetics_data[rxnID][source]
                    if not _check_datum(datum):
                        continue
                    ### define rate law variables; calculate flux; average or overwrite the flux based on data criteria
                    locals().update(
                        {
                            metID: self.conc.at[metID, self.previous_time] * milli
                            for metID in datum["mets"]
                        }
                    )
                    flux = eval(datum["substituted_rate_law"])
                    print(datum["substituted_rate_law"], flux)
                    if (
                        "metadata" not in self.kinetics_data[rxnID][source]
                        or self.__find_data_match(rxnID, source) == "a"
                    ):
                        fluxes.append(flux)
                    else:
                        fluxes = [flux]

                flux = mean(fluxes)
                rxn = self.model_util.model.reactions.get_by_id(rxnID)
                rxn.lb = rxn.ub = flux
                self.fluxes.at[rxnID, self.time] = flux
            ## execute the COBRA model
            sol = self.model_util.model.optimize()
            self.sols.append(sol)
            ## add previously undefined fluxes and concentrations
            for rxnID in self.fluxes.index:
                if self.fluxes.at[rxnID, self.time] == 0:
                    self.fluxes.at[rxnID, self.time] = sol.fluxes[rxnID]
            for met in self.model_util.model.metabolites:
                self.conc.at[met.id, self.time] = 0
                for rxn in met.reactions:
                    flux = self.fluxes.at[rxn.id, self.time]
                    if flux == 0:
                        continue
                    # print(rxn.metabolites[met], flux, timestep_hr, cell_g_L)
                    self.conc.at[met.id, self.time] += (
                        rxn.metabolites[met] * flux * timestep_hr * cell_g_L
                    )
            if all(chemostat_requirements):
                self.moles[self.time] = self.conc[self.time] * milli * chemostat_L
                self._chemostat(feed_profile, chemostat_L_hr, chemostat_L)
            elif any(chemostat_requirements):
                warn(
                    "The < chemostat_L > , < feed_profile >, and < chemostat_L_hr >"
                    " parameters must all be defined to simulate a chemostat."
                )
            self.variables["elapsed_time"] += self.ts_min
            if self.printing:
                print(
                    f"\nObjective value (\u0394t{self.ts_min}): ",
                    self.sols[-1].objective_value,
                )

        # identify the chemicals that dynamically changed in concentrations
        self.changed = set(
            [
                met_id
                for met_id in self.met_ids
                if self.conc.at[met_id, "0 min"] != self.conc.at[met_id, self.time]
            ]
        )
        self.unchanged = set(self.met_ids.keys()) - self.changed

        # visualize concentration changes over time
        if visualize:
            self._visualize(conc_figure_title, included_mets, labeled_plots)
        if export:
            self._export(export_name, export_directory, total_min)
        if self.verbose:
            print(
                f"\nChanged concentrations:\t{self.changed}",
                f"\nConstrained reactions:\t{constrained.keys()}",
            )
        elif self.printing:
            if self.jupyter:
                pandas.set_option("max_rows", None)
                display(self.conc, self.fluxes)
            if self.unchanged == set():
                print(
                    "All of the metabolites changed concentration over the simulation"
                )
            else:
                print(f"\nUnchanged metabolite concentrations\t{self.unchanged}")
        return self.conc, self.fluxes

    def _chemostat(self, feed_profile: dict, chemostat_L_hr, chemostat_L):
        L_changed = chemostat_L_hr * self.ts_min
        # chemostat addition
        for met_id, conc in feed_profile.items():
            self.moles.at[met_id, self.time] += conc * L_changed
            self.conc.at[met_id, self.time] = (
                self.moles.at[met_id, self.time] / milli / chemostat_L
            )  # normalize to the chemostat volume
        # chemostat subtraction
        for met in self.model_util.model.metabolites:
            if met.compartment[0] != "e":
                continue
            ## update the chemical moles
            self.moles.at[met.id, self.time] -= (
                self.conc.at[met.id, self.time] * L_changed
            )
            ## define the chemical concentration
            self.conc.at[met.id, self.time] = (
                self.moles.at[met.id, self.time] / milli / chemostat_L
            )

    # nested functions
    def __find_data_match(self, rxnID: str, source: str):
        # identifies the datum whose experimental conditions most closely matches the simulation conditions
        temperature_deviation = ph_deviation = 0
        if FBAHelper.isnumber(
            self.kinetics_data[rxnID][source]["metadata"]["Temperature"]
        ):
            temp = float(self.kinetics_data[rxnID][source]["metadata"]["Temperature"])
            temperature_deviation = (
                abs(self.parameters["temperature"] - temp)
                / self.parameters["temperature"]
            )
        if FBAHelper.isnumber(self.kinetics_data[rxnID][source]["metadata"]["pH"]):
            pH = float(self.kinetics_data[rxnID][source]["metadata"]["pH"])
            ph_deviation = abs(self.parameters["pH"] - pH) / self.parameters["pH"]

        # equally weight between temperature and pH deviation from the simulation conditions
        old_minimum = self.minimum
        deviation = mean(temperature_deviation, ph_deviation)
        self.minimum = min(deviation, self.minimum)
        return (
            "a" if old_minimum == self.minimum else "w"
        )  # append or write a list of data

    def _visualize(self, conc_fig_title, included_mets, labeled_plots):
        # TODO construct a Vega visualization with a range bind that permits scanning over a time series
        ## and accordingly adjusting arrowhead widths to reflect flux at the particularly timestep.
        ## The heatmap may likewise be dynamic for each timestep over a bind range.

        # define the figure
        pyplot.rcParams["figure.figsize"] = (11, 7)
        pyplot.rcParams["figure.dpi"] = 150
        self.figure, ax = pyplot.subplots()
        ax.set_title(conc_fig_title)
        ax.set_ylabel("Concentrations (mM)")

        x_axis_scalar, unit = _x_axis_determination(self.total_min)
        ax.set_xlabel("Time " + unit)
        legend_list = []
        times = [
            t * self.ts_min * x_axis_scalar
            for t in range(self.parameters["timesteps"] + 1)
        ]

        # determine the plotted metabolites and the scale of the figure axis
        bbox = (1, 1)
        if not included_mets:
            bbox = (1.7, 1)
            # 1e-2 is an arbitrary concentration threshold for plotting on the figure
            included_mets = [
                chem
                for chem in self.changed
                if max(self.conc.loc[[chem]].values[0].tolist()) > 1e-2
            ]

        log_axis = False
        minimum, maximum = inf, -inf
        printed_concentrations = {}
        for chem in self.changed:
            if chem not in included_mets:
                continue
            concentrations = self.conc.loc[[chem]].values[0].tolist()
            maximum = max(maximum, max([x if x > 1e-9 else 0 for x in concentrations]))
            minimum = min(minimum, min([x if x > 1e-9 else 0 for x in concentrations]))
            # plot chemicals with perturbed concentrations
            ax.plot(times, concentrations)
            if len(chem) > 25:
                chem = list(self.met_ids.keys())[self.met_ids.index(chem)]
            if not concentrations[0] < 1e-9:
                legend_list.append(chem)
            else:
                legend_list.append(f"(rel) {chem}")

            # design the proper location of the overlaid labels in the figure
            if not labeled_plots:
                continue
            for i, conc in enumerate(concentrations):
                if conc <= 1e-9:
                    continue
                x_value = i * self.ts_min
                vertical_adjustment = 0
                if x_value in printed_concentrations:
                    vertical_adjustment = (maximum - minimum) * 0.05
                    if log_axis:
                        vertical_adjustment = log10(maximum - minimum) / 3
                ax.text(
                    x_value,
                    conc + vertical_adjustment,
                    f"{chem} - {round(conc, 4)}",
                    ha="left",
                )
                printed_concentrations[x_value] = conc
                break

        # finalize figure details
        if maximum > 10 * minimum:
            ax.set_yscale("log")
        ax.set_xticks(times)
        ax.grid(True)
        ax.legend(
            legend_list,
            title="Changed chemicals",
            loc="upper right",
            bbox_to_anchor=bbox,
            title_fontsize="x-large",
            fontsize="large",
        )

    def _export(self, export_name="kineticsFBA", export_directory: str = None):
        # define a unique simulation name
        directory = (
            os.path.dirname(export_directory) if export_directory else os.getcwd()
        )
        self.parameters["simulation_path"] = self.simulation_path = os.path.join(
            directory, export_name
        )
        # export simulation content
        self.fluxes.to_csv(os.path.join(self.simulation_path, "fluxes.csv"))
        self.conc.to_csv(os.path.join(self.simulation_path, "concentrations.csv"))
        obj_vals_df = pandas.DataFrame(
            [
                (self.fluxes.columns[index].replace(" min", ""), sol.objective_value)
                for index, sol in enumerate(self.sols)
            ],
            columns=["min", "objective_value"],
        )
        obj_vals_df.index = obj_vals_df["min"]
        obj_vals_df.drop(["min"], axis=1, inplace=True)
        obj_vals_df.to_csv(os.path.join(self.simulation_path, "objective_values.csv"))
        # export the parameters
        parameters_table = pandas.DataFrame(
            self.parameters, columns=["parameter", "value"]
        )
        parameters_table.to_csv(os.path.join(self.simulation_path, "parameters.csv"))
        # export the figure
        self.figure.savefig(
            os.path.join(self.simulation_path, "changed_concentrations.svg")
        )
        if self.verbose and not self.jupyter:
            self.figure.show()
