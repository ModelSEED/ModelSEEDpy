from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.minimalmediapkg import MinimalMediaPkg
from modelseedpy.community.mscommunity import MSCommunity
from itertools import combinations, permutations, chain
from optlang import Variable, Constraint, Objective
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.exceptions import ObjectiveError
from cobra.medium import minimal_medium
from collections import Counter
from deepdiff import DeepDiff  # (old, new)
from typing import Iterable
from pprint import pprint
from numpy import array
from numpy import mean
from math import prod


def _compatibilize_models(cobra_models: Iterable, printing=False):
    return cobra_models
    # return MSCompatibility.standardize(cobra_models, conflicts_file_name='exchanges_conflicts.json', printing=printing)


class MSSmetana:
    def __init__(self, cobra_models: Iterable, com_model, min_growth=0.1, n_solutions=100,
                 abstol=1e-3, media_dict=None, printing=True, compatibilize=True, minimal_media_method="jenga"):

        self.min_growth = min_growth ; self.abstol = abstol ; self.n_solutions = n_solutions
        self.printing = printing ; self.compatibilize = compatibilize

        self.models, self.community = MSSmetana._load_models(cobra_models, com_model, compatibilize)
        for model in self.models:
            if model.slim_optimize() == 0:
                raise ObjectiveError(f"The model {model.id} possesses an objective value of 0 in complete media, "
                                     "which is incompatible with minimal media computations and hence SMETANA.")
        if self.community.slim_optimize() == 0:
            raise ObjectiveError(f"The community model {self.community.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations and hence SMETANA.")
        if not media_dict:
            # if cobra_models:
            #     self.media = MinimalMediaPkg.interacting_comm_media(self.models, printing=self.printing)
            # else:
            if minimal_media_method == "minFlux":
                self.media = MinimalMediaPkg.minimize_flux(self.community, min_growth)
            elif minimal_media_method == "minComponent":
                self.media = MinimalMediaPkg.minimize_components(self.community, min_growth)
            elif minimal_media_method == "jenga":
                self.media = MinimalMediaPkg.jenga_method(self.community, self.models)
        else:
            self.media = media_dict
        self.media = {"community_media": self.media, "members": {}}

    def mro_score(self):
        self.mro = MSSmetana.mro(self.models, self.min_growth, self.media, self.compatibilize, self.printing)
        if self.printing:
            print(f"\nMRO score: {self.mro}\t\t\t{self.mro*100:.2f}% of member requirements, on average, overlap with other members "
                  "requirements.")
        return self.mro

    def mip_score(self, interacting_media:dict=None, noninteracting_media:dict=None):
        diff, self.mip = MSSmetana.mip(self.models, self.min_growth, interacting_media,
                                       noninteracting_media, compatibilize=self.compatibilize)
        if self.printing:
            print(f"\nMIP score: {self.mip}\t\t\t{self.mip} required compound(s) can be sourced via syntrophy:")
            pprint(diff)
        return self.mip

    def mu_score(self):
        self.mu = MSSmetana.mu(self.models, self.n_solutions, self.abstol, compatibilize=self.compatibilize)
        if self.printing:
            print("\nMU score:\t\t\tThe fraction of solutions in which each member is the "
                  "syntrophic receiver that contain a respective metabolite:\n")
            pprint(self.mu)
        return self.mu

    def mp_score(self):
        self.mp = MSSmetana.mp(self.models, self.community, self.abstol, compatibilize=self.compatibilize)
        if self.printing:
            print("\nMP score:\t\t\tThe possible contributions of each member in the member media include:\n")
            pprint(self.mp)
        return self.mp

    def sc_score(self):
        self.sc = MSSmetana.sc(self.models, self.community, self.min_growth,
                               self.n_solutions, self.abstol, compatibilize=self.compatibilize)
        if self.printing:
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
            self.models, self.community, self.min_growth, self.n_solutions, self.abstol,
            (self.sc, self.mu, self.mp), self.compatibilize, sc_coupling)
        if self.printing:
            print("\nsmetana score:\n")
            pprint(self.smetana)
        return self.smetana

    ###### STATIC METHODS OF THE SMETANA SCORES, WHICH ARE APPLIED IN THE ABOVE CLASS OBJECT ######

    @staticmethod
    def _load_models(cobra_models: Iterable, com_model=None, compatibilize=True):
        if not com_model and cobra_models:
            return cobra_models, MSCommunity.build_from_species_models(cobra_models, name="SMETANA_example", cobra_model=True)
        # models = PARSING_FUNCTION(community_model) # TODO the individual models of a community model must be parsed
        if compatibilize:
            return _compatibilize_models(cobra_models), _compatibilize_models([com_model])[0]
        return cobra_models, com_model

    @staticmethod
    def _get_media(media=None, com_model=None, model_s_=None, min_growth=None,
                   syntrophy=True, printing=False, minimization_method="jenga"):
        if not media:  # model_s_ is either a singular model or a list of models
            if com_model or not isinstance(model_s_, (list,set,tuple)):
                return MinimalMediaPkg.jenga_method(
                    com_model or model_s_, minimal_growth=min_growth, syntrophy=syntrophy, printing=printing)
            elif isinstance(model_s_, (list,set,tuple)):
                if syntrophy:
                    return MinimalMediaPkg.interacting_comm_media(
                        model_s_, minimization_method, None, printing)
                return MinimalMediaPkg.comm_media_est(model_s_, 0.1, minimization_method, printing)
            raise TypeError("Either the com_model or model_s_ arguments must be parameterized.")
        # and's here to quit after a false, unlike all()
        if model_s_ and not isinstance(model_s_, (list,set,tuple)) and model_s_.id not in media["members"]:
            media["members"][model_s_.id] = {"media": MinimalMediaPkg.jenga_method(
                    model_s_, minimal_growth=min_growth, syntrophy=syntrophy, printing=printing)}
            return media["members"][model_s_.id]["media"]
        return media["community_media"]

    @staticmethod
    def mro(cobra_models:Iterable, min_growth=0.1, media_dict=None, compatibilize=True, printing=False):
        """Determine the maximal overlap of minimal media between member organisms."""
        cobra_models = _compatibilize_models(cobra_models) if compatibilize else cobra_models
        ind_media = {model.id: set(MSSmetana._get_media(media_dict, None, model, min_growth, printing)) for model in cobra_models}
        pairs = {(model1, model2): ind_media[model1.id] & ind_media[model2.id] for model1, model2 in combinations(cobra_models, 2)}
        # ratio of the average size of intersecting minimal media between any two members and the minimal media of all members
        return mean(list(map(len, pairs.values()))) / mean(list(map(len, ind_media.values())))

    @staticmethod
    def mip(cobra_models:Iterable, min_growth=0.1, interacting_media_dict=None,
            noninteracting_media_dict=None, compatibilize=True):
        """Determine the maximum quantity of nutrients that can be sourced through syntrophy"""
        cobra_models = _compatibilize_models(cobra_models) if compatibilize else cobra_models
        noninteracting_medium = MSSmetana._get_media(noninteracting_media_dict, None, cobra_models, min_growth, False)
        interacting_medium = MSSmetana._get_media(interacting_media_dict, None, cobra_models, min_growth, True)
        nonint_med = noninteracting_medium if "community_media" not in noninteracting_medium else noninteracting_medium["community_media"]
        int_med = interacting_medium if "community_media" not in interacting_medium else interacting_medium["community_media"]
        interact_diff = DeepDiff(nonint_med, int_med)
        if "dictionary_item_removed" in interact_diff:
            return interact_diff["dictionary_item_removed"], len(interact_diff["dictionary_item_removed"])
        return None, 0

    @staticmethod
    def mu(cobra_models:Iterable, n_solutions=100, abstol=1e-3, compatibilize=True):
        """the fractional frequency of each received metabolite amongst all possible alternative syntrophic solutions"""
        # determine the solutions for each member
        # member_solutions = member_solutions if member_solutions else {model.id: model.optimize() for model in cobra_models}
        scores = {}
        if compatibilize:
            cobra_models = _compatibilize_models(cobra_models)
        for org_model in cobra_models:
            model = org_model.copy()
            ex_rxns = {ex_rxn: met for ex_rxn in FBAHelper.exchange_reactions(model) for met in ex_rxn.metabolites}
            variables = {ex_rxn.id: Variable('___'.join([model.id, ex_rxn.id]), lb=0, ub=1, type="binary") for ex_rxn in ex_rxns}
            FBAHelper.add_cons_vars(model, [list(variables.values())])
            media, solutions = [], []
            for i in range(0, n_solutions):
                if i > 0:
                    constraint = Constraint(sum([variables[ex.id] for ex in medium]), ub=len(medium)-1, name=f"iteration_{i}")
                    FBAHelper.add_cons_vars(model, [constraint])
                sol = model.optimize()
                if sol.status != 'optimal':
                    break
                # determine the de facto medium for this simulated growth
                solutions.append(sol)
                medium = set([ex_rxn for ex_rxn in ex_rxns if sol.fluxes[ex_rxn.id] < -abstol])
                media.append(medium)
            counter = Counter(chain(*media))
            scores[model.id] = {met.id: counter[ex] / len(media) for ex, met in ex_rxns.items() if counter[ex] > 0}
        return scores

    @staticmethod
    def mp(cobra_models:Iterable, com_model=None, abstol=1e-3, compatibilize=True):
        """Discover the metabolites that each species can contribute to a community"""
        cobra_models, community = MSSmetana._load_models(cobra_models, com_model, compatibilize)
        com_model = community.copy()
        scores = {}
        for org_model in cobra_models:
            model = org_model.copy()
            scores[model.id] = []
            # determines possible member contributions in the community environment, where the excretion of media compounds is irrelevant
            approximate_minimal_media = MinimalMediaPkg.comm_media_est(cobra_models)
            possible_contributions = [ex_rxn for ex_rxn in FBAHelper.exchange_reactions(model) if ex_rxn.id not in approximate_minimal_media]
            while len(possible_contributions) > 0:
                print("remaining possible_contributions", len(possible_contributions), end="\r")
                FBAHelper.add_objective(com_model, sum(ex_rxn.flux_expression for ex_rxn in possible_contributions))
                sol = com_model.optimize()
                fluxes_contributions, uncertain_contributions = [], []
                for ex in possible_contributions:
                    if ex.id in sol.fluxes.keys():
                        if sol.fluxes[ex.id] >= abstol:
                            fluxes_contributions.append(ex)
                            possible_contributions.remove(ex)

                if sol.status != 'optimal' or not fluxes_contributions:
                    break
                # log confirmed contributions
                for ex_rxn in fluxes_contributions:
                    for met in ex_rxn.metabolites:
                        scores[model.id].append(met.id)

            # double-check the remaining possible contributions for excretion
            for ex_rxn in possible_contributions:
                com_model.objective = Objective(ex_rxn.flux_expression)
                sol = com_model.optimize()
                if sol.status != 'optimal' or sol.objective_value < abstol:
                    for met in ex_rxn.metabolites:
                        if met.id in scores[model.id]:
                            print("removing", met.id)
                            scores[model.id].remove(met.id)
        return scores

    @staticmethod
    def sc(cobra_models:Iterable=None, com_model=None, min_growth=0.1, n_solutions=100, abstol=1e-6, compatibilize=True):
        """Calculate the frequency of interspecies dependency in a community"""
        cobra_models, community = MSSmetana._load_models(cobra_models, com_model, compatibilize)
        com_model = community.copy()
        for rxn in com_model.reactions:
            rxn.lower_bound = 0 if 'bio' in rxn.id else rxn.lower_bound

        # c_{rxn.id}_lb: rxn < 1000*y_{species_id}
        # c_{rxn.id}_ub: rxn > -1000*y_{species_id}
        variables = {}
        constraints = []
        # TODO this can be converted to an MSCommunity object by looping through each index
        for org_model in cobra_models:
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
        for model in cobra_models:
            com_model = com_model.copy()
            other_members = [other for other in cobra_models if other.id != model.id]
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
    def smetana(cobra_models: Iterable, com_model=None, min_growth=0.1, n_solutions=100, abstol=1e-6,
                prior_values=None, compatibilize=False, sc_coupling=False):
        """Quantifies the extent of syntrophy as the sum of all exchanges in a given nutritional environment"""
        cobra_models, community = MSSmetana._load_models(cobra_models, com_model, compatibilize)
        com_model = community.copy()
        sc = None
        if not prior_values:
            mu = MSSmetana.mu(cobra_models, n_solutions, abstol, compatibilize)
            mp = MSSmetana.mp(cobra_models, com_model, abstol, compatibilize)
            if sc_coupling:
                sc = MSSmetana.sc(cobra_models, com_model, min_growth, n_solutions, abstol, compatibilize)
        elif len(prior_values) == 3:
            sc, mu, mp = prior_values
        else:
            mu, mp = prior_values

        smetana_scores = {}
        for pairs in combinations(cobra_models, 2):
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
                            smetana_scores[model1.id][model2.id] += prod([mu[model1.id].get(met.id,0), sc_score, mp_score])
        return smetana_scores