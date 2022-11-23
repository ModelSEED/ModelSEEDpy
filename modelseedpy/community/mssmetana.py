import numpy
from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.msminimalmedia import MSMinimalMedia
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.community.mscommunity import MSCommunity
from itertools import combinations, permutations, chain
from optlang import Variable, Constraint, Objective
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.exceptions import ObjectiveError, ParameterError
from modelseedpy.core.msmodelutl import MSModelUtil
from cobra.medium import minimal_medium
from collections import Counter
from deepdiff import DeepDiff  # (old, new)
from typing import Iterable
from pprint import pprint
from numpy import array
from icecream import ic
from numpy import mean
# from math import prod


def _compatibilize(member_models: Iterable, printing=False):
    # return member_models
    models = MSCompatibility.standardize(member_models, conflicts_file_name='exchanges_conflicts.json', printing=printing)
    if not isinstance(member_models, (set, list, tuple)):
        return models[0]
    return models

def _load_models(member_models: Iterable, com_model=None, compatibilize=True, printing=False):
    ic(member_models, com_model, compatibilize)
    if not com_model and member_models:
        model, names, abundances = build_from_species_models(member_models, name="SMETANA_example")
        return member_models, model  # (model, names=names, abundances=abundances)
    # models = PARSING_FUNCTION(community_model) # TODO the individual models of a community model can be parsed
    if compatibilize:
        # return _compatibilize(member_models, printing), MSCommunity(_compatibilize([com_model], printing)[0])
        return _compatibilize(member_models, printing), _compatibilize([com_model], printing)[0]
    return member_models, com_model  # MSCommunity(com_model)

def _get_media(media=None, com_model=None, model_s_=None, min_growth=None,
               interacting=True, printing=False, minimization_method="minComponents"):
    ic(media, com_model, model_s_)
    if not com_model and not model_s_:
        raise TypeError("Either the com_model or model_s_ arguments must be parameterized.")
    if media:
        if model_s_ and not isinstance(model_s_, (list,set,tuple)):
            return media["members"][model_s_.id]["media"]
        elif com_model:
            return media["community_media"]
        return media
    # model_s_ is either a singular model or a list of models
    if com_model:
        com_media = MSMinimalMedia.determine_min_media(com_model, minimization_method, min_growth, interacting, printing)
    if model_s_:
        if not isinstance(model_s_, (list,set,tuple)):
            return MSMinimalMedia.determine_min_media(model_s_, minimization_method, min_growth, interacting, printing)
        if isinstance(model_s_, (list,set,tuple)):
            members_media = {}
            for model in model_s_:
                members_media[model.id] = {"media":MSMinimalMedia.determine_min_media(
                    model, minimization_method, min_growth, interacting, printing)}
            if not com_model:
                return members_media
    else:
        return com_media
    return {"community_media":com_media, "members":members_media}


class MSSmetana:
    def __init__(self, member_models: Iterable, com_model, min_growth=0.1, n_solutions=100, environment=None,
                 abstol=1e-3, media_dict=None, printing=True, raw_content=False, antismash_json_path:str=None,
                 antismash_zip_path:str=None, minimal_media_method="minComponents"):
        self.min_growth = min_growth ; self.abstol = abstol ; self.n_solutions = n_solutions
        self.printing = printing ; self.raw_content = raw_content
        self.antismash_json_path = antismash_json_path ; self.antismash_zip_path = antismash_zip_path

        # process the models
        self.models = _compatibilize(member_models)
        self.community = MSModelUtil(com_model or build_from_species_models(self.models, cobra_model=True))
        ## define the environment
        if environment:
            if hasattr(environment, "get_media_constraints"):
                ### standardize modelseed media into COBRApy media
                environment = {"EX_" + exID: -bound[0] for exID, bound in environment.get_media_constraints().items()}
            self.community.add_medium(environment)
        self.environment = environment
        ## test growth
        for model in self.models:
            if model.slim_optimize() == 0:
                raise ObjectiveError(f"The model {model.id} possesses an objective value of 0 in complete media, "
                                     "which is incompatible with minimal media computations and hence SMETANA.")
        if self.community.model.slim_optimize() == 0:
            raise ObjectiveError(f"The community model {self.community.model.id} possesses an objective value of 0 in complete "
                                 f"media, which is incompatible with minimal media computations and hence SMETANA.")
        ## determine the minimal media for each model, including the community
        self.media = media_dict or _get_media(
                None, self.community.model, self.models, min_growth, True, printing, minimal_media_method)

    def all_scores(self):
        mro = self.mro_score()
        mip = self.mip_score(interacting_media=self.media)
        mp = self.mp_score()
        mu = self.mu_score()
        sc = self.sc_score()
        smetana = self.smetana_score()
        return (mro, mip, mp, mu, sc, smetana)

    def mro_score(self):
        self.mro = MSSmetana.mro(self.models, self.media["members"], self.min_growth,
                                 self.media, self.raw_content, self.printing, True)
        if not self.printing:
            return self.mro
        if self.raw_content:
            for pair, (interaction, media) in self.mro.items():
                newcomer, established = pair.split('---')
                print(f"\n(MRO) The {newcomer} media {media} possesses {interaction} shared "
                      f"requirements with the {established} established member.")
                return self.mro
        for pair, mro in self.mro.items():
            newcomer, established = pair.split('---')
            print(f"\nThe {newcomer} on {established} MRO score: {mro[0]} ({mro[0]*100:.2f}%). "
                  f"This is the percent of nutritional requirements in {newcomer} "
                  f"that overlap with {established} ({mro[1]}/{mro[2]}).")
        return self.mro

    def mip_score(self, interacting_media:dict=None, noninteracting_media:dict=None):
        diff, self.mip = MSSmetana.mip(self.community.model, self.models, self.min_growth, interacting_media,
                                       noninteracting_media, self.printing, True)
        if not self.printing:
            return self.mip
        print(f"\nMIP score: {self.mip}\t\t\t{self.mip} required compound(s) can be sourced via syntrophy:")
        if self.raw_content:
            pprint(diff)
        return self.mip

    def mp_score(self):
        self.mp = MSSmetana.mp(self.models, self.environment, self.community.model, self.abstol, self.printing)
        if not self.printing:
            return self.mp
        if self.raw_content:
            print("\n(MP) The possible contributions of each member in the member media include:\n")
            pprint(self.mp)
        else:
            print("\nMP score:\t\t\tEach member can possibly contribute the following to the community:\n")
            for member, contributions in self.mp.items():
                print(member, "\t", len(contributions))
        return self.mp

    def mu_score(self):
        member_excreta = self.mp_score() if not hasattr(self, "mp") else self.mp
        self.mu = MSSmetana.mu(self.models, self.n_solutions, member_excreta, self.n_solutions,
                               self.abstol, True, self.printing)
        if not self.printing:
            return self.mu
        print("\nMU score:\t\t\tThe fraction of solutions in which each member is the "
              "syntrophic receiver that contain a respective metabolite:\n")
        pprint(self.mu)
        return self.mu

    def sc_score(self):
        self.sc = MSSmetana.sc(self.models, self.community.model, self.min_growth,
                               self.n_solutions, self.abstol, True, self.printing)
        if not self.printing:
            return self.sc
        print("\nSC score:\t\t\tThe fraction of community members who syntrophically contribute to each species:\n")
        pprint(self.sc)
        return self.sc

    def smetana_score(self):
        if not hasattr(self, "sc"):
            self.sc = self.sc_score()
        sc_coupling = all(array(list(self.sc.values())) != None)
        if not hasattr(self, "mu"):
            self.mu = self.mu_score()
        if not hasattr(self, "mp"):
            self.mp = self.mp_score()

        self.smetana = MSSmetana.smetana(
            self.models, self.community.model, self.min_growth, self.n_solutions, self.abstol,
            (self.sc, self.mu, self.mp), True, sc_coupling, self.printing)
        if self.printing:
            print("\nsmetana score:\n")
            pprint(self.smetana)
        return self.smetana

    def antiSMASH_scores(self):
        if not hasattr(self, "antismash_json_path"):
            raise TypeError("The antismash_json_path argument must be specified to conduct the antiSMASH scores.")
        self.antismash = MSSmetana.antiSMASH(self.antismash_json_path)
        if not self.printing:
            return self.antismash
        if self.raw_content:
            print("\n(antismash) The biosynthetic_areas, BGCs, protein_annotations, clusterBlast, and "
                  "num_clusterBlast from the provided antiSMASH results:\n")
            print("The 'areas' that antiSMASH determines produce biosynthetic products:")
            pprint(self.antismash[0])
            print("The set of biosynthetic gene clusters:")
            pprint(self.antismash[1])
            print("The set of clusterblast protein annotations:")
            pprint(self.antismash[2])
            print("Resistance information from clusterblast")
            pprint(self.antismash[3])
            print("The number of proteins associated with resistance")
            pprint(self.antismash[4])
            return self.antismash
        print("\nantiSMASH scores:\n")
        print("The community exhibited:"
              f"- {len(self.antismash[0])}'areas' that antiSMASH determines produce biosynthetic products."
              f"- {len(self.antismash[1])} biosynthetic gene clusters."
              f"- {len(self.antismash[2])} clusterblast protein annotations."
              f"- {len(self.antismash[3])} parcels of resistance information from clusterblast."
              f"- {self.antismash[4]} proteins associated with resistance.")
        return list(map(len, self.antismash[:4]))+[self.antismash[4]]


    ###### STATIC METHODS OF THE SMETANA SCORES, WHICH ARE APPLIED IN THE ABOVE CLASS OBJECT ######

    @staticmethod
    def mro(member_models:Iterable=None, mem_media:dict=None, min_growth=0.1, media_dict=None,
            raw_content=False, printing=False, compatibilized=False):
        """Determine the overlap of nutritional requirements (minimal media) between member organisms."""
        # determine the member minimal media if they are not parameterized
        if not mem_media:
            if not member_models:
                raise ParameterError("The either member_models or minimal_media parameter must be defined.")
            member_models = _compatibilize(member_models, printing)
            mem_media = _get_media(media_dict, None, member_models, min_growth, printing=printing)
            if "community_media" in mem_media:
                mem_media = mem_media["members"]
        # MROs = array(list(map(len, pairs.values()))) / array(list(map(len, mem_media.values())))
        mro_values = {}
        for model1, model2 in permutations(member_models, 2):
            intersection = set(mem_media[model1.id]["media"].keys()) & set(mem_media[model2.id]["media"].keys())
            member_media = mem_media[model1.id]["media"]
            if raw_content:
                mro_values[f"{model1.id}---{model2.id})"] = (intersection, member_media)
            else:
                mro_values[f"{model1.id}---{model2.id})"] = (
                    len(intersection)/len(member_media), len(intersection), len(member_media))
        return mro_values
        # return mean(list(map(len, pairs.values()))) / mean(list(map(len, mem_media.values())))

    @staticmethod
    def mip(com_model=None, member_models:Iterable=None, min_growth=0.1, interacting_media_dict=None,
            noninteracting_media_dict=None, printing=True, compatibilized=False):
        """Determine the maximum quantity of nutrients that can be sourced through syntrophy"""
        member_models, community = _load_models(
            member_models, com_model, compatibilized==False, printing=printing)
        # determine the interacting and non-interacting media for the specified community  .util.model
        noninteracting_medium = _get_media(noninteracting_media_dict, community, None, min_growth, False)
        if "community_media" in noninteracting_medium:
            noninteracting_medium = noninteracting_medium["community_media"]
        interacting_medium = _get_media(interacting_media_dict, community, None, min_growth, True)
        if "community_media" in interacting_medium:
            interacting_medium = interacting_medium["community_media"]
        # differentiate the community media
        interact_diff = DeepDiff(noninteracting_medium, interacting_medium)
        if "dictionary_item_removed" in interact_diff:
            return interact_diff["dictionary_item_removed"], len(interact_diff["dictionary_item_removed"])
        return None, 0

    @staticmethod
    def contributions(org_possible_contributions, scores, model_util, abstol):
        # identify and log excreta from the solution
        FBAHelper.add_objective(model_util.model, sum(
            ex_rxn.flux_expression for ex_rxn in org_possible_contributions))
        sol = model_util.model.optimize()
        if sol.status != "optimal":
            ## exit the while loop by returning the original possible_contributions, hence DeepDiff == {} and the while loop terminates
            return scores, org_possible_contributions
        # identify and log excreta from the solution
        possible_contributions = org_possible_contributions[:]
        for ex in org_possible_contributions:
            if ex.id in sol.fluxes.keys() and sol.fluxes[ex.id] >= abstol:
                possible_contributions.remove(ex)
                scores[model_util.model.id].update([met.id for met in ex.metabolites])
        return scores, possible_contributions

    @staticmethod
    def mp(member_models:Iterable, environment, com_model=None, minimal_media=None, abstol=1e-3, printing=False):
        """Discover the metabolites that each species can contribute to a community"""
        community = _compatibilize(com_model) or build_from_species_models(member_models, cobra_model=True, standardize=True)
        # TODO support parsing the individual members through the MSCommunity object
        scores = {}
        for org_model in member_models:
            model_util = MSModelUtil(org_model)
            model_util.compatibilize(printing=printing)
            if environment:
                model_util.add_medium(environment)
            # TODO leverage extant minimal media as the default instead of the community complete media
            scores[model_util.model.id] = set()
            # determines possible member contributions in the community environment, where the excretion of media compounds is irrelevant
            org_possible_contributions = [ex_rxn for ex_rxn in model_util.exchange_list()
                                          if ex_rxn.id not in community.medium and ex_rxn.upper_bound > 0]
            ic(org_possible_contributions, len(model_util.exchange_list()), len(community.medium))
            scores, possible_contributions = MSSmetana.contributions(
                org_possible_contributions, scores, model_util, abstol)
            while DeepDiff(org_possible_contributions, possible_contributions):
                print("remaining possible_contributions", len(possible_contributions), end="\r")
                ## optimize the sum of the remaining exchanges that have not surpassed the abstol
                org_possible_contributions = possible_contributions[:]
                scores, possible_contributions = MSSmetana.contributions(
                    org_possible_contributions, scores, model_util, abstol)

            ## individually checks the remaining possible contributions
            for ex_rxn in possible_contributions:
                model_util.model.objective = Objective(ex_rxn.flux_expression)
                sol = model_util.model.optimize()
                if sol.status == 'optimal' or sol.objective_value > abstol:
                    for met in ex_rxn.metabolites:
                        if met.id in scores[model_util.model.id]:
                            print("removing", met.id)
                            scores[model_util.model.id].remove(met.id)
        return scores

    @staticmethod
    def mu(member_models:Iterable, environment=None, member_excreta=None, n_solutions=100, abstol=1e-3,
           compatibilized=False, printing=True):
        """the fractional frequency of each received metabolite amongst all possible alternative syntrophic solutions"""
        # member_solutions = member_solutions if member_solutions else {model.id: model.optimize() for model in member_models}
        scores = {}
        member_models = member_models if compatibilized else _compatibilize(member_models, printing)
        member_excreta = member_excreta or MSSmetana.mp(member_models, environment, None, abstol, printing)
        missing_members = [model for model in member_models if model.id not in member_excreta]
        if missing_members:
            member_excreta.update(MSSmetana.mp([missing_members], environment))
        ic(member_models)
        for org_model in member_models:
            other_excreta = set(numpy.array([excreta for model, excreta in member_excreta.items()
                                            if model.id != org_model.id]).flatten())
            model_util = MSModelUtil(org_model)
            if environment:
                model_util.add_medium(environment)
            ex_rxns = {ex_rxn: met for ex_rxn in model_util.exchange_list() for met in ex_rxn.metabolites}
            variables = {ex_rxn.id: Variable('___'.join([model_util.model.id, ex_rxn.id]),
                                             lb=0, ub=1, type="binary") for ex_rxn in ex_rxns}
            FBAHelper.add_cons_vars(model_util.model, [list(variables.values())])
            media, solutions = [], []
            sol = model_util.model.optimize()
            while sol.status == "optimal" and len(solutions) < n_solutions:
                solutions.append(sol)
                medium = set([ex for ex in ex_rxns if sol.fluxes[ex.id] < -abstol and ex in other_excreta])
                constraint = Constraint(sum([variables[ex.id] for ex in medium]),
                                        ub=len(medium)-1, name=f"iteration_{len(solutions)}")
                FBAHelper.add_cons_vars(model_util.model, [constraint])
                media.append(medium)
                sol = model_util.model.optimize()
            counter = Counter(chain(*media))
            scores[model_util.model.id] = {met.id: counter[ex] / len(media)
                                           for ex, met in ex_rxns.items() if counter[ex] > 0}
        return scores

    @staticmethod
    def sc(member_models:Iterable=None, com_model=None, min_growth=0.1, n_solutions=100, abstol=1e-6,
           compatibilized=True, printing=False):
        """Calculate the frequency of interspecies dependency in a community"""
        member_models, community = _load_models(
            member_models, com_model, compatibilized==False, printing=printing)
        for rxn in com_model.reactions:
            rxn.lower_bound = 0 if 'bio' in rxn.id else rxn.lower_bound

        # c_{rxn.id}_lb: rxn < 1000*y_{species_id}
        # c_{rxn.id}_ub: rxn > -1000*y_{species_id}
        variables = {}
        constraints = []
        # TODO this can be converted to an MSCommunity object by looping through each index
        for org_model in member_models:
            model = org_model.copy()
            variables[model.id] = Variable(name=f'y_{model.id}', lb=0, ub=1, type='binary')
            FBAHelper.add_cons_vars(com_model, [variables[model.id]])
            for rxn in model.reactions:
                if "bio" not in rxn.id:
                    # print(rxn.flux_expression)
                    lb = Constraint(rxn.flux_expression + 1000*variables[model.id], name="_".join(["c", model.id, rxn.id, "lb"]), lb=0)
                    ub = Constraint(rxn.flux_expression - 1000*variables[model.id], name="_".join(["c", model.id, rxn.id, "ub"]), ub=0)
                    constraints.extend([lb, ub])
        FBAHelper.add_cons_vars(com_model, constraints, sloppy=True)

        # calculate the SCS
        scores = {}
        for model in member_models:
            com_model = com_model.copy()
            other_members = [other for other in member_models if other.id != model.id]
            # model growth is guaranteed while minimizing the growing members of the community
            ## SMETANA_Biomass: {biomass_reactions} > {min_growth}
            smetana_biomass = Constraint(sum(rxn.flux_expression for rxn in model.reactions if "bio" in rxn.id),
                                         name='SMETANA_Biomass', lb=min_growth)
            FBAHelper.add_cons_vars(com_model, [smetana_biomass], sloppy=True)
            FBAHelper.add_objective(com_model, sum([variables[other.id] for other in other_members]), "min")
            previous_constraints, donors_list = [], []
            for i in range(n_solutions):
                sol = com_model.optimize()  # FIXME The solution is not optimal
                if sol.status != 'optimal':
                    scores[model.id] = None
                    break
                donors = [o for o in other_members if com_model.solver.primal_values[f"y_{o.id}"] > abstol]
                donors_list.append(donors)
                previous_con = f'iteration_{i}'
                previous_constraints.append(previous_con)
                FBAHelper.add_cons_vars(com_model, [Constraint(
                    sum(variables[o.id] for o in donors), name=previous_con, ub=len(previous_constraints)-1)], sloppy=True)
            if i != 0:
                donors_counter = Counter(chain(*donors_list))
                scores[model.id] = {o.id: donors_counter[o] / len(donors_list) for o in other_members}
        return scores

    @staticmethod
    def smetana(member_models: Iterable, environment, com_model=None, min_growth=0.1, n_solutions=100, abstol=1e-6,
                prior_values=None, compatibilized=False, sc_coupling=False, printing=False):
        """Quantifies the extent of syntrophy as the sum of all exchanges in a given nutritional environment"""
        member_models, community = _load_models(
            member_models, com_model, compatibilized==False, printing=printing)
        sc = None
        if not prior_values:
            mp = MSSmetana.mp(member_models, environment, com_model, abstol)
            mu = MSSmetana.mu(member_models, environment, mp, n_solutions, abstol, compatibilized)
            if sc_coupling:
                sc = MSSmetana.sc(member_models, com_model, min_growth, n_solutions, abstol, compatibilized)
        elif len(prior_values) == 3:
            sc, mu, mp = prior_values
        else:
            mu, mp = prior_values

        smetana_scores = {}
        for pairs in combinations(member_models, 2):
            for model1, model2 in permutations(pairs):
                if model1.id not in smetana_scores:
                    smetana_scores[model1.id] = {}
                if not any([not mu[model1.id], not mp[model1.id]]):
                    sc_score = 1 if not sc_coupling else sc[model1.id][model2.id]
                    models_mets = list(model1.metabolites)+list(model2.metabolites)
                    unique_mets = set([met.id for met in models_mets])
                    smetana_scores[model1.id][model2.id] = 0
                    for met in models_mets:
                        if met.id in unique_mets:
                            mp_score = 0 if met.id not in mp[model1.id] else 1
                            smetana_scores[model1.id][model2.id] += mu[model1.id].get(met.id,0)*sc_score*mp_score
        return smetana_scores

    @staticmethod
    def antiSMASH(json_path=None, zip_path=None):
        # TODO Scores 2, 4, and 5 are being explored for relevance to community formation and reveal specific member interactions/targets
        # load the antiSMASH report from either the JSON or the raw ZIP, or both
        from os import mkdir, listdir, path
        from zipfile import ZipFile
        from json import load
        if json_path:
            cwd_files = listdir()
            if json_path not in cwd_files and zip_path:
                with ZipFile(zip_path, "r") as zip_file:
                    zip_file.extract(json_path)
            with open(json_path, "r") as json_file:
                data = load(json_file)
        elif zip_path:
            mkdir("extracted_antiSMASH")
            with ZipFile(zip_path, "r") as zip_file:
                zip_file.extractall("extracted_antiSMASH")
            json_files = [x for x in listdir("extracted_antiSMASH") if x.endswith("json")]
            if len(json_files) > 1:
                print(f"The antiSMASH report describes {len(json_files)} JSON files, the first of which is selected "
                      f"{json_files[0]} for analysis, otherwise explicitly identify the desired JSON file in the json_path parameter.")
            with open(path.join("extracted_antiSMASH", json_files[0]), "r") as json_file:
                data = load(json_file)
        else:
            raise ParameterError("Either the json_path or zip_path from the antiSMASH analysis must be provided,"
                                 " for these scores to be determined.")
        # Parse data and scores from the antiSMASH report
        biosynthetic_areas = data["records"][0]['areas']
        BGCs = set(numpy.array([data["records"][0]['areas'][i]['products'] for i in range(biosynthetic_areas)]).flatten())
        len_proteins = len(data["records"][0]['modules']['antismash.modules.clusterblast']['knowncluster']['proteins'])
        protein_annotations = [data["records"][0]['modules']['antismash.modules.clusterblast']['knowncluster']['proteins'][i]['annotations']
                           for i in range(len_proteins)]
        clusterBlast = [s for s in protein_annotations if "resistance" in s]
        num_clusterBlast = sum([item.count("resistance") for item in protein_annotations])

        return biosynthetic_areas, BGCs, protein_annotations, clusterBlast, num_clusterBlast
