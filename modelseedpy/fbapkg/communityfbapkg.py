def CommunityFBAPkg(modelInfo, mediaInfo, kbase, element_uptake_limit = None, kinetic_coeff = None, abundances = None, msdb_path_for_fullthermo = None):
    # import general modules
    from itertools import combinations
    from numpy import unique
    import cplex
    import re
    
    # import the model and media
    from modelseedpy.fbapkg import kbasemediapkg
    model = kbase.get_from_ws(modelInfo[0],modelInfo[1])
    media = kbase.get_from_ws(mediaInfo[0],mediaInfo[1])
    kmp = kbasemediapkg.KBaseMediaPkg(model)
    kmp.build_package(media)
    model.solver = 'optlang-cplex'
    
    # applying uptake constraints
    if element_uptake_limit is not None:
        from modelseedpy.fbapkg import elementuptakepkg
        eup = elementuptakepkg.ElementUptakePkg(model)
        eup.build_package(element_uptake_limit)
    
    # applying kinetic constraints
    if kinetic_coeff is not None:
        from modelseedpy.fbapkg import commkineticpkg
        ckp = commkineticpkg.CommKineticPkg(model)
        ckp.build_package(kinetic_coeff,abundances)
    
    # applying FullThermo constraints
    if msdb_path_for_fullthermo is not None:
        from modelseedpy.fbapkg import fullthermopkg
        ftp = fullthermopkg.FullThermoPkg(model)
        ftp.build_package({'modelseed_path':msdb_path_for_fullthermo})

    # unambiguously defining the model objective as biomass growth
    biomass_objective = model.problem.Objective(
        1 * model.reactions.bio1.flux_expression,
        direction='max')
    model.objective = biomass_objective

    # excute FBA
    solution = model.optimize()
    '''pfba_solution = cobra.flux_analysis.pfba(model)'''
    
    # calculate the metabolic exchanges
    metabolite_uptake = {}
    compartment_numbers = []
    for rxn in model.reactions:
        if (rxn.id[-2] == 'c' or rxn.id[-2] == 'p') and rxn.id[-1] != '0':
            compartment_number = rxn.id[-1]
            compartment_numbers.append(compartment_number)

            for metabolite in rxn.metabolites:
                if metabolite.compartment == "e0":
                    rate_law = 0
                    flux = solution.fluxes[rxn.id]
                    if flux != 0:
                        rate_law += rxn.metabolites[metabolite]*flux
                        metabolite_uptake[(metabolite.id,compartment_number)] = rate_law

    compartment_numbers = unique(compartment_numbers)
    number_of_compartments = len(compartment_numbers)

    # calculate cross feeding of extracellular metabolites      # production[donor_compartment_number][receiver_compartment_number]
    production = [[0] * number_of_compartments for loop in range(number_of_compartments)] 
    consumption = [[0] * number_of_compartments for loop in range(number_of_compartments)]
    boundary_fluxes = []
    cross_all = []
    for rxn in model.reactions:
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
                if total_in > total_out:
                    rate_out[None] = total_in - total_out
                if total_in < total_out:
                    rate_in[None] = total_out - total_in

                # establish the metabolites that partake in cross feeding 
                for donor, rate_1 in rate_out.items():
                    if donor is None:
                        receiver = list(rate_in.keys())[0]
                        boundary_fluxes.append(receiver)
                    else:
                        donor_index = int(donor) - 1
                        
                        for receiver, rate_2 in rate_in.items():
                            if receiver is None:
                                boundary_fluxes.append(donor)
                            else:
                                receiver_index = int(receiver) - 1
                        
                                # assign calculated feeding rates to the production and consumption matrices
                                rate = rate_1 * rate_2 / max_total_rate
                                production[donor_index][receiver_index] += rate
                                consumption[receiver_index][donor_index] += rate

                                
    # graph the resultant community interaction       
    from matplotlib import pyplot
    import networkx
    
    graph = networkx.Graph()
    for num in compartment_numbers:
        graph.add_node(num)
    for com in combinations(compartment_numbers, 2):
        species_1 = int(com[0])-1
        species_2 = int(com[1])-1
        
        interaction_flux = max(production[species_1][species_2], consumption[species_2][species_1])
        graph.add_edge(com[0],com[1],weight=interaction_flux)

    pos = networkx.circular_layout(graph)
    networkx.draw_networkx(graph,pos)
    labels = networkx.get_edge_attributes(graph,'weight')
    networkx.draw_networkx_edge_labels(graph,pos,edge_labels=labels)
    pyplot.show()
    return production, consumption, graph