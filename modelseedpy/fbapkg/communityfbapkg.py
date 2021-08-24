class CommunityFBAPkg():
    def __init__(self, modelInfo, mediaInfo, kbase):
        '''GENERAL IMPORT '''
        self.parameters = {}
        self.variables = {}
        
        import cplex

        # import the model and media
        from modelseedpy.fbapkg import kbasemediapkg
        self.model = kbase.get_from_ws(modelInfo[0],modelInfo[1])
        media = kbase.get_from_ws(mediaInfo[0],mediaInfo[1])
        kmp = kbasemediapkg.KBaseMediaPkg(self.model)
        kmp.build_package(media)
        self.model.solver = 'optlang-cplex'

        # unambiguously defining the model objective as biomass growth
        biomass_objective = self.model.problem.Objective(
            1 * self.model.reactions.bio1.flux_expression,
            direction='max')
        self.model.objective = biomass_objective
    

    def constrain(self, element_uptake_limit = None, kinetic_coeff = None, abundances = None, msdb_path_for_fullthermo = None, lp_file = False):
        ''' APPLY CONSTRAINTS '''
        # applying uptake constraints
        element_contraint_name = ''
        if element_uptake_limit is not None:
            from modelseedpy.fbapkg import elementuptakepkg

            element_contraint_name = 'eup'
            eup = elementuptakepkg.ElementUptakePkg(self.model)
            eup.build_package(element_uptake_limit)

        # applying kinetic constraints
        kinetic_contraint_name = ''
        if kinetic_coeff is not None:
            from modelseedpy.fbapkg import commkineticpkg

            kinetic_contraint_name = 'ckp'
            ckp = commkineticpkg.CommKineticPkg(self.model)
            ckp.build_package(kinetic_coeff,abundances)

        # applying FullThermo constraints
        thermo_contraint_name = ''
        if msdb_path_for_fullthermo is not None:
            from modelseedpy.fbapkg import fullthermopkg

            thermo_contraint_name = 'ftp'
            ftp = fullthermopkg.FullThermoPkg(self.model)
            ftp.build_package({'modelseed_path':msdb_path_for_fullthermo})


        # conditionally print the LP file of the model
        if lp_file:
            from os.path import exists
            import re

            count_iteration = 0
            file_name = '_'.join([modelInfo[0], thermo_contraint_name, kinetic_contraint_name, element_contraint_name, str(count_iteration)])
            file_name += '.lp'
            while exists(file_name):
                count_iteration += 1
                file_name = re.sub('(_\d)', f'_{count_iteration}', file_name)

            with open(file_name, 'w') as out:
                out.write(str(self.model.solver))
                out.close()

    def calculate(self):
        ''' CALCULATE FLUXES ''' 
        # execute the model
        solution = self.model.optimize()
        '''pfba_solution = cobra.flux_analysis.pfba(model)'''

        # calculate the metabolic exchanges
        metabolite_uptake = {}
        self.variables['compartment_numbers'] = []
        for rxn in self.model.reactions:
            if (rxn.id[-2] == 'c' or rxn.id[-2] == 'p') and rxn.id[-1] != '0':
                compartment_number = rxn.id[-1]
                self.variables['compartment_numbers'].append(compartment_number)

                for metabolite in rxn.metabolites:
                    if metabolite.compartment == "e0":
                        flux = solution.fluxes[rxn.id]
                        if flux != 0:
                            rate_law = rxn.metabolites[metabolite]*flux
                            metabolite_uptake[(metabolite.id, compartment_number)] = rate_law

        # create empty matrices      # production[donor_compartment_number][receiver_compartment_number]
        from numpy import zeros

        number_of_compartments = len(set(self.variables['compartment_numbers']))
        self.variables['production'] = zeros((number_of_compartments, number_of_compartments)) 
        self.variables['consumption'] = zeros((number_of_compartments, number_of_compartments))

        # calculate cross feeding of extracellular metabolites   
        cross_all = []
        for rxn in self.model.reactions:
            for metabolite in rxn.metabolites:
                if metabolite.compartment == "e0":
                    # determine each directional flux rate 
                    rate_out = {compartment_number: rate for (metabolite_id, compartment_number), rate in metabolite_uptake.items() if metabolite_id == metabolite.id and rate > 0}
                    rate_in = {compartment_number: abs(rate) for (metabolite_id, compartment_number), rate in metabolite_uptake.items() if metabolite_id == metabolite.id and rate < 0}

                    # determine total directional flux rate 
                    total_in = sum(rate_in.values())
                    total_out = sum(rate_out.values())
                    max_total_rate = max(total_in, total_out)

                    # determine net flux 
                    net_flux = total_in - total_out
                    if net_flux > 0:
                        rate_out[None] = net_flux
                    if net_flux < 0:
                        rate_in[None] = abs(net_flux)


                    # establish the metabolites that partake in cross feeding 
                    for donor, rate_1 in rate_out.items():
                        if donor is not None:
                            donor_index = int(donor) - 1           
                            for receiver, rate_2 in rate_in.items():
                                if receiver is not None:
                                    receiver_index = int(receiver) - 1

                                    # assign calculated feeding rates to the production and consumption matrices
                                    rate = rate_1 * rate_2 / max_total_rate
                                    self.variables['production'][donor_index][receiver_index] += rate
                                    self.variables['consumption'][receiver_index][donor_index] += rate

        return self.variables['production'], self.variables['consumption']
    

    def visualize(self, graph = True, table = True):
        ''' VISUALIZE FLUXES ''' 
        # graph the community network
        if graph:
            from itertools import combinations
            import networkx
            
            graph = networkx.Graph()
            for num in self.variables['compartment_numbers']:
                graph.add_node(num)
            for com in combinations(self.variables['compartment_numbers'], 2):
                species_1 = int(com[0])-1
                species_2 = int(com[1])-1

                interaction_net_flux = round(self.variables['production'][species_1][species_2] - self.variables['consumption'][species_1][species_2], 4)
                if species_1 < species_2:
                    graph.add_edge(com[0],com[1],flux = interaction_net_flux)
                elif species_1 > species_2:
                    graph.add_edge(com[0],com[1],flux = -interaction_net_flux)

            pos = networkx.circular_layout(graph)
            networkx.draw_networkx(graph,pos)
            labels = networkx.get_edge_attributes(graph,'flux')
            networkx.draw_networkx_edge_labels(graph,pos,edge_labels=labels)
        
        
        # view the cross feeding matrices
        if table:
            from pandas import DataFrame as df

            species = [num+1 for num in range(len(self.variables['production']))]

            print('\nProduction matrix:')
            prod_df = df(self.variables['production'])
            prod_df.index = prod_df.columns = species
            prod_df.index.name = 'Donor'
            print(prod_df)

            print('\n\nConsumption matrix:')
            cons_df = df(self.variables['consumption'])
            cons_df.index = cons_df.columns = species
            cons_df.index.name = 'Receiver'
            print(cons_df)