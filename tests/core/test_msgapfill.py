# -*- coding: utf-8 -*-
"""
from glob import glob
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = 'xx'
kbase = cobrakbase.KBaseAPI(token)
import re

# define the example individual model and associated API media package
model = kbase.get_from_ws('e_coli_core.kb', 95098)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
import modelseedpy
from modelseedpy.core.msgapfill import MSGapfill
gapfill = MSGapfill(model)

def test_init():
    assert type(gapfill.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert type(gapfill.blacklist) is list
    assert type(gapfill.solutions) is dict

def test_run_gapfilling_and_integrate_gapfill_solution():
    solutions = gapfill.run_gapfilling()

    # test that the objective expression is correctly set
    if solutions is not None:
        assert type(solutions) is dict

        # verify the integrate_gapfill_solution function
        model_2 = gapfill.integrate_gapfill_solution(solutions)
        assert type(model_2) is cobrakbase.core.kbasefba.fbamodel.FBAModel

        for reaction in solutions['reversed']:
            if solution["reversed"][reaction] == ">":
                assert reaction.upper_bound == 100
            else:
                assert reaction.lower_bound == -100

        for reaction in solutions['new']:
            if solution["new"][reaction] == ">":
                assert reaction.upper_bound == 100
                assert reaction.lower_bound == 0
            else:
                assert reaction.upper_bound == 0
                assert reaction.lower_bound == -100

def test_gapfill():
    pass
"""
import os
import pytest
import json
import cobra
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy import MSGapfill, MSMedia


@pytest.fixture
def template():
    with open(
        os.path.join(
            os.path.dirname(__file__), "..", "test_data", "template_core_bigg.json"
        ),
        "r",
    ) as fh:
        return MSTemplateBuilder.from_dict(json.load(fh)).build()


@pytest.fixture
def get_model():
    def _method(ko=None):
        if ko is None:
            ko = []
        with open(
            os.path.join(
                os.path.dirname(__file__), "..", "test_data", "e_coli_core.json"
            ),
            "r",
        ) as fh:
            model_json = json.load(fh)
            model_json["compartments"] = {
                k + "0": v for (k, v) in model_json["compartments"].items()
            }
            metabolites = {}
            for m in model_json["metabolites"]:
                m["id"] += "0"
                m["compartment"] += "0"
                metabolites[m["id"]] = m
            for r in model_json["reactions"]:
                r["metabolites"] = {i + "0": v for (i, v) in r["metabolites"].items()}
                compartments = set(
                    [metabolites[k]["compartment"] for k in r["metabolites"].keys()]
                )
                if r["id"].endswith("_e"):
                    r["id"] += "0"
                elif len(compartments) == 1:
                    r["id"] += "_" + list(compartments)[0]
                else:
                    r["id"] += (
                        "_" + "c0"
                    )  # hack cause there is only combo between e0 and c0

            model_json["reactions"] = [
                x for x in model_json["reactions"] if x["id"] not in ko
            ]
            return cobra.io.from_json(json.dumps(model_json))

    return _method


@pytest.fixture
def media_glucose_aerobic():

    return MSMedia.from_dict(
        {
            "glc__D": (-10, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
            "pi": (-1000, 1000),
            "co2": (-1000, 1000),
            "nh4": (-1000, 1000),
        }
    )


def test_model_default(get_model):
    solution = get_model([]).optimize()
    assert solution.status == "optimal"
    assert solution.objective_value > 0.8


def test_model_no_solution(get_model):
    solution = get_model(["GLCpts_c0"]).optimize()
    assert solution.status == "infeasible"


def test_ms_gap_fill1(template, get_model, media_glucose_aerobic):
    """
    Test gap filling with glucose aerobic requesting minimal growth value
    """
    model = get_model(["CYTBD_c0", "NADH16_c0", "GLCpts_c0", "O2t_c0"])
    gap_fill = MSGapfill(model, [template])
    result = gap_fill.run_gapfilling(
        media_glucose_aerobic, "BIOMASS_Ecoli_core_w_GAM_c0", minimum_obj=0.1
    )
    assert result
    assert "new" in result
    assert "GLCpts_c0" in result["new"]
    assert result["new"]["GLCpts_c0"] == ">"


def test_ms_gap_fill2(template, get_model, media_glucose_aerobic):
    """
    Test gap filling with glucose aerobic requesting higher growth value
    """
    model = get_model(["CYTBD_c0", "NADH16_c0", "GLCpts_c0", "O2t_c0"])
    gap_fill = MSGapfill(model, [template])
    result = gap_fill.run_gapfilling(
        media_glucose_aerobic, "BIOMASS_Ecoli_core_w_GAM_c0", minimum_obj=0.8
    )
    assert result
    assert "new" in result
    assert "GLCpts_c0" in result["new"] and result["new"]["GLCpts_c0"] == ">"
    assert "CYTBD_c0" in result["new"] and result["new"]["CYTBD_c0"] == ">"
    assert "NADH16_c0" in result["new"] and result["new"]["NADH16_c0"] == ">"
    assert "O2t_c0" in result["new"] and result["new"]["O2t_c0"] == ">"
