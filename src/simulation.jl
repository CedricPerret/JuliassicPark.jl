#***********************************************
#*** User-facing functions (simulation control)
#***********************************************

using Distributed

"""
    evol_model(parameters_input, fitness_function, repro_function; 
               sweep = Dict{Symbol, Vector}(), 
               additional_parameters = Dict(), 
               migration_function = nothing, 
               genotype_to_phenotype_mapping = identity)

Main entry point to run an evolutionary simulation or parameter sweep with user-defined components.
 
# Arguments
- `parameters_input::Dict`: Dictionary specifying all simulation parameters. Must include at least `:z_ini` (initial trait generator or values). Other required parameters depend on trait type (e.g. boundaries for continuous or discrete traits).
- `fitness_function::Function`: User-defined function computing fitness from traits or populations. Must:
    - Return fitness values as the first output.
    - Accept parameters as keyword arguments (`; args...`).
    - Optionally return additional outputs for logging.
    - Accept `should_it_print::Bool` to limit costly computations.
- `repro_function::Function`: Function defining reproduction. May be in-place (`!`) or return-based. See `list_reproduction_methods()`.

# Keyword Arguments
- `sweep::Dict{Symbol, <:AbstractVector}`: Optional parameter sweep. Each key is a parameter symbol; values are tested in combination.
- `additional_parameters::Dict`: Optional dictionary of derived or fixed parameters to be merged in and saved to output.
- `migration_function::Union{Function, Nothing}`: Optional function defining migration across patches or groups.
- `genotype_to_phenotype_mapping::Function`: Optional mapping function, mainly for sexual reproduction. Defaults to `identity`.

# Returns
- `DataFrame` or `Vector{DataFrame}`: Simulation result(s) containing traits, fitness values, and optional extra outputs, depending on `:write_file`, `:split_simul`, and `:split_sweep` settings.
    - If `:write_file = true`, results are written to `.csv` files instead of returned.
    - If `:n_simul == 0`, the simulation function is returned without execution.

# Notes
- Supports both single runs and full parameter sweeps via `sweep`.
- The behaviour of saving vs. returning data is controlled by:
    - `:write_file` — write to disk or return
    - `:split_simul` — one file per simulation
    - `:split_sweep` — one file per parameter set
- Parallelisation:
    - If `:split_sweep` and `:split_simul` are both `true`, each sweep × simulation is parallelised and saved independently.
    - If no writing is requested and `:n_simul` is large, threads are used to speed up simulation replicates.
- Trait types (Boolean, Integer, Float) are inferred from `:z_ini` and related fields (`:boundaries`, `:sigma_m`, etc.).
- Additional statistics returned by `fitness_function` are automatically logged, with names controlled by `:other_output_names`.
"""

## Wrapper (prepare parameters and put default, then call each part)
## here prepare only thing that will not be varied in a sweep.
function evol_model(parameters_input, fitness_function, repro_function; sweep=Dict{Symbol, Vector}(), additional_parameters= Dict{Symbol, Function}(), migration_function = nothing, genotype_to_phenotype_mapping = identity)
    
    ## To not modify the parameters given
    if isempty(parameters_input)
        parameters = parse_commandline()
    else
        parameters = deepcopy(parameters_input)
    end
    parameters = parameters isa NamedTuple ? Dict{Any,Any}(pairs(parameters)) : parameters
    parameters = merge(get_default_parameters(), parameters)

    @assert haskey(parameters, :z_ini) "Missing required parameter: `:z_ini`. Please provide the initial values or generators for the trait(s)."

    ## We infer if the group sizes will change based on the reproduction function.
    if  contains(give_me_my_name(repro_function), "explicit") 
        #-> Group sizes are not constant.
        parameters[:n_cst] = false
    else
        parameters[:n_cst] = true
    end
    ## Parameter that will never be printed.
    append!(parameters[:parameters_to_omit],PARAMETERS_TO_ALWAYS_OMIT)
    parameters[:parameters_to_omit]= Symbol.(parameters[:parameters_to_omit])

    ## Generate the model to output (as shown in replicator, it takes parameters and its ID which is i_simul as input)
    model = get_template_model(parameters, fitness_function, repro_function; additional_parameters= additional_parameters, migration_function = migration_function, genotype_to_phenotype_mapping = genotype_to_phenotype_mapping)
    
    #--- Run
    ## If no sweep, it gives back a vector containing as single element the parameters
    run_parameter_sweep_distributed(model, sweep, parameters)

end


"""
    get_template_model(parameters_input, fitness_function, repro_function; additional_parameters=Dict(), migration_function=nothing, genotype_to_phenotype_mapping=identity)

Returns a function that runs a single evolutionary simulation with a fixed set of parameters.

This wrapper assembles all components—initialisation, fitness evaluation, mutation, reproduction, migration, and output—into a callable function `model(parameters, i_simul)` that executes a full simulation run.

# Arguments
- `parameters_input::Dict`: Dictionary of input parameters (copied internally).
- `fitness_function::Function`: User-supplied function computing fitness from phenotypes.
- `repro_function::Function`: Function specifying how reproduction is performed (in-place or return-based).
- `additional_parameters::Dict`: Derived parameters to be computed at initialisation.
- `migration_function::Function`: Optional migration step applied after reproduction.
- `genotype_to_phenotype_mapping`: Function mapping genotypes to phenotypes (default: `identity`).

# Returns
A function `model(parameters, i_simul)` that:
- Runs the simulation for the specified number of generations.
- Uses `i_simul` to set the random seed (for reproducibility).
- Returns a `DataFrame` containing the results.

# Notes
- If `genotype_to_phenotype_mapping` is not provided and `repro_function` is sexual, default mappings are used:
  - For single-locus: mean genotype (average).
  - For multilocus: additive mapping with allelic effect `:delta` (must be provided).
- Only relevant keyword arguments are passed to mutation, reproduction, and migration to ensure flexibility and modularity.
"""

function get_template_model(parameters_input, fitness_function, repro_function; additional_parameters= Dict{Symbol, Function}(), migration_function = nothing, genotype_to_phenotype_mapping = identity)
    model = function(parameters_input, i_simul)
        #*** Initialisation
        Random.seed!(i_simul)
        parameters = deepcopy(parameters_input)
        #--- Initialise second parameters which need to be derived from the given parameters (which can be directly printed)
        parameters,cst_output_name,cst_output = compute_derived_parameters!(parameters,additional_parameters;additional_parameters_to_omit=parameters[:additional_parameters_to_omit])

        n_traits = 1
        if parameters[:z_ini] isa AbstractDataFrame
            n_traits = length(filter(c -> startswith(string(c), "z"), names(parameters[:z_ini])))
            push!(parameters[:parameters_to_omit], :z_ini)
        ## Standardise z_ini
        else
            if !isa(parameters[:z_ini], Tuple)
                #-> A single trait but z_ini is not a tuple. We ensure `z_ini` is always a tuple, even if the user provides a scalar (e.g. 0.5).
                parameters[:z_ini] = tuple(parameters[:z_ini])
            end
            n_traits = length(parameters[:z_ini])
        end
        ## If it continues previous simulations, the first generation is already contained in the previous dataframe, so we do "skip it"
        (parameters[:z_ini] isa AbstractDataFrame && parameters[:n_print] == 1) && (parameters[:n_print] += 1)

        #--- Normalise the trait relevant parameters (depending of number of traits, need to unwrap or copy) + reorganise them in the right format
        mut_kwargs_names = [:sigma_m, :boundaries, :mutation_type,:bias_m] 
        _normalise_trait_parameters!(parameters,[:mu_m; mut_kwargs_names], n_traits)
        ## Move from sigma_m = [sigma_m1,sigma_m2], bias_m = [bias_m1,bias_m2] to param_trait1 = [sigma_m = sigma_m1, bias_m=bias_m1]
        ## which can be given as a whole to mutation
        mut_kwargs = _prepare_kwargs_multiple_traits(parameters, mut_kwargs_names, n_traits)

        #--- Initialise the population 
        population = initialise_population(parameters[:z_ini], parameters[:n_ini], parameters[:n_patch]; boundaries = parameters[:boundaries], simplify = parameters[:simplify],n_loci = parameters[:n_loci])

        #--- This is in case the user gives a single generator, not encapsulated in a vector because the expected input is a vector of generators for each trait. 
        if contains(give_me_my_name(repro_function), "sexual") && genotype_to_phenotype_mapping == identity
            @assert all(T <: Float64 for T in _type_of_traits(parameters[:z_ini])) "Sexual reproduction with diploids is only implemented for continuous traits (Float64)"
            if parameters[:n_loci] == 1
                println("With sexual reproduction, the genotype to phenotype mapping needs to be provided. In the case single locus, we assume average mapping with no dominance")
                genotype_to_phenotype_mapping = average_mapping
            else
                println("With sexual reproduction, the genotype to phenotype mapping needs to be provided. In the case multilocus, we assume additive effects with two possible discrete alleles")
                if !haskey(parameters, :delta)
                    error("Please provide the value of a single allelic effect")
                end
                genotype_to_phenotype_mapping = x -> additive_mapping(x,parameters[:delta])
            end
        end
        calculate_phenotype = genotype_to_phenotype_mapping === identity ? identity : (pop -> lift_map(genotype_to_phenotype_mapping, pop))

        #--- Initialise type of reproduction
        reproduction_mode = :pure
        new_population = deepcopy(population)
        if contains(give_me_my_name(repro_function), "!")
            if contains(give_me_my_name(repro_function), "WF")
                reproduction_mode = :inplace_double
            else
                reproduction_mode = :inplace_single
            end
        end
        in_place_fitness_function = contains(give_me_my_name(fitness_function), "!") ? true : false

        #--- Preprocess fitness function
        population_phenotype = calculate_phenotype(population)
        ## Standardise the output of the fitness function. See preprocess_fitness_function for details.
        instanced_fitness_function = nothing; output = Any[];
        if !in_place_fitness_function
            correction = _infer_fitness_function_correction(population,fitness_function, parameters,calculate_phenotype)
            instanced_fitness_function = preprocess_fitness_function(population_phenotype, fitness_function, parameters,correction)
            output = [[population_phenotype]; instanced_fitness_function(population_phenotype; parameters...)]
        else
            correction = -1
            ## In-place fitness functions don’t return fitness values, so we cannot infer the correction.
            ## Instead, check the function signature to ensure it matches the population structure.
            @assert _infer_fitness_function_correction_by_signature(population, fitness_function) == 0  """
            In-place fitness functions must have the whole population as trait input  
            - If your whole population is a vector of individuals, the fitness function must work on a vector. 
            - If your whole population is a vector of vectors (metapopulation), the fitness function must work on a vector of vectors.
            """
            instanced_fitness_function = preprocess_fitness_function(population_phenotype, fitness_function, parameters,correction)
            fitness = vv(0.,population_phenotype)
            output = [[population_phenotype,fitness]; instanced_fitness_function(population_phenotype, fitness; parameters...)]
        end

        fit = (output[2] isa AbstractVector && first(output[2]) isa AbstractVector) ? Iterators.flatten(output[2]) : output[2]
        @assert all(x -> x isa AbstractFloat, fit) "Fitness values must be floating-point"
        
        ## If fitness function returns a named tuple and no other output name specified, we get the names of the extras output.
        other_output_names = !isempty(parameters[:other_output_names]) ? parameters[:other_output_names] : extract_output_names(population_phenotype, fitness_function, parameters,correction)
        ## Truncate if too many names were provided
        if length(other_output_names) > length(output)-2
            extras_ignored = join(other_output_names[(length(output) - 1):end], ", ", " and ")
            resize!(other_output_names,length(output)-2)
            println("Too many extra output names were provided; the extras $extras_ignored are ignored")
        end
        
        df_res, saver = init_data_output(
            parameters[:de], [["z", "fitness"]; other_output_names],
            output, parameters[:n_gen], parameters[:n_print], parameters[:j_print],
            i_simul, parameters[:n_patch], parameters[:n_ini], parameters[:n_cst];
            output_cst_names=cst_output_name, output_cst=cst_output, n_traits
        )

        ## To avoid passing unnecessary or unused parameters (which could cause errors or reduce clarity), we explicitly filter only the relevant keyword arguments from the main `parameters` dictionary.
        repro_kwargs_names = [:n_replacement,:transition_proba, :n_pop_by_class, :n_patch, :mig_rate,:group_fitness_fun,:n_loci]
        repro_kwargs = _filter_kwargs(parameters,repro_kwargs_names)

        mig_kwargs_names = [:mig_rate]
        mig_kwargs = _filter_kwargs(parameters,mig_kwargs_names)

        #*** Run simulations
        for i_gen in 1:parameters[:n_gen]

            population_phenotype = calculate_phenotype(population)
            output[1] = population_phenotype
            #--- Calculate fitness
            if !in_place_fitness_function
                #@ To avoid reallocating memory
                #@ Note that it means that the output keep the old extras when not calculating them but it is fine since we print/use extras only when calculating them.
                res = instanced_fitness_function(population_phenotype; parameters..., should_it_print=should_it_print(i_gen, parameters[:n_print], parameters[:j_print]))
                if should_it_print(i_gen, parameters[:n_print], parameters[:j_print])
                    #-> Update all the output
                    output[2:(1+length(res))] = res
                else
                    #-> Update only fitness (the only element required for reproduction)
                    output[2] = res[1]
                end
            else
            #-> In-place fitness function
                res = instanced_fitness_function(population_phenotype, output[2]; parameters..., should_it_print=should_it_print(i_gen, parameters[:n_print], parameters[:j_print]))
                #@ Note that it means that the output keep the old extras when not calculating them but it is fine since we print/use extras only when calculating them.
                output[3:end] = res
            end
            #--- Save
            saver(df_res, i_gen, output)
            #--- Reproduce
            if reproduction_mode == :pure
                #->repro_function gives back a new population which needs to be assigned (slower)
                population = repro_function(population, output[2], parameters[:str_selection], parameters[:mu_m], mut_kwargs; repro_kwargs...)
            elseif reproduction_mode == :inplace_single
                #->repro_function is in-place (faster)
                repro_function(population, output[2], parameters[:str_selection], parameters[:mu_m], mut_kwargs; repro_kwargs...)
            elseif reproduction_mode == :inplace_double
                #->repro_function is in-place (faster)
                repro_function(population, new_population,output[2], parameters[:str_selection], parameters[:mu_m], mut_kwargs; repro_kwargs...)
                population, new_population = new_population, population 
            end
            #--- Check if all patches became empty 
            if my_isempty(population)
                error("Population collapsed")
            end
            #--- Migrate
            if migration_function !== nothing
                if contains(give_me_my_name(migration_function), "!")
                    #-> In-place
                    migration_function(population, new_population; mig_kwargs...)
                    population, new_population = new_population, population 
                else
                    #-> Allocating
                    population = migration_function(population; mig_kwargs...)
                end
            end
        end
        if parameters[:z_ini] isa AbstractDataFrame
            df_first = parameters[:z_ini]
            df_res.gen .+= maximum(df_first.gen)
            ## We duplicated the last/first generation so we remove it here.
            append!(df_first, df_res; promote=true, cols=:union)        
            ## Homogeneise the id_simul    
            df_first.i_simul .= df_first.i_simul[1]
            df_res = df_first
        end
        return df_res
    end
    return model
end


"""
    run_parameter_sweep_distributed(fun, sweep, parameters)

Runs a distributed parameter sweep over multiple parameter combinations and simulation replicates.

Handles different parallelisation strategies based on user-specified flags such as `:write_file`, `:split_sweep`, and `:split_simul`.

# Arguments
- `fun`: A function of the form `(parameters, seed) -> DataFrame`, returning results from a single simulation.
- `sweep::Dict{Symbol, Vector}`: Dictionary of parameter values to sweep over.
- `parameters::Dict`: Base parameters for all simulations. Must include fields like `:n_simul`, `:write_file`, and `:name_model`.

# Parallelisation Modes
- If `:write_file=false`: all results are accumulated and returned as a `DataFrame` or list of `DataFrame`s.
- If `:write_file=true` and `:split_sweep=true`: results for each parameter combination are saved to separate files.
- If `:split_simul=true`: each replicate is saved individually, with simulation ID in the filename.
- Uses multi-threading when possible (unless resolution is individual-level or memory usage is high).

# Returns
- If `:write_file=false`: concatenated `DataFrame` (or list of `DataFrame`s if `:split_sweep=true`).
- If `:write_file=true`: writes results to disk and returns `nothing`.

# Errors
- Throws an error if `:split_simul=true` but `:split_sweep=false`, as this may lead to file name conflicts.

# Notes
- Output filenames are generated using `build_filepath(...)`, optionally summarising swept parameter values.
- Uses progress bars and thread-local buffers for efficient aggregation under multi-threading.
"""
function run_parameter_sweep_distributed(fun, sweep, parameters)
    list_parameters_set, sweep_df = get_parameters_from_sweep(parameters, sweep)
    ##function to calculate a random seed
    n = length(list_parameters_set)
    total = n*parameters[:n_simul]
    p = Progress(total, 1)
    update!(p, 0)  # <-- Forces the bar to show immediately
    if parameters[:split_simul] && !parameters[:split_sweep] && !isempty(sweep)
        error("split_simul=true requires split_sweep=true to avoid file conflicts or ambiguity.")
    end

    if !parameters[:write_file] || (parameters[:write_file] && !parameters[:split_simul] && !parameters[:split_sweep])
        # -> No parallel write, accumulate everything
        #--- Generate data
        list_res = Vector{DataFrame}(undef, n)
        if Threads.nthreads() == 1 || parameters[:n_simul] < 10 || parameters[:de] =='i'
            #-> Small enough or take too much memory to paralelise over threads.
            for i in 1:n
                res = DataFrame()
                for i_simul in 1:parameters[:n_simul]
                    id_simul =  get_simulation_seed(parameters, i_simul)
                    sim_res = fun(list_parameters_set[i],id_simul)
                    append!(res, sim_res)
                    next!(p)
                end
                if n > 1
                    list_res[i] = hcat(repeat(sweep_df[i:i, :], inner = nrow(res)),res )
                else
                    list_res[i] = res
                end
            end
        else
            #-> Use threads to parallelise over simulations
            for i in 1:n
                simulation_results = Vector{DataFrame}(undef, parameters[:n_simul])

                Threads.@threads for i_simul in eachindex(simulation_results)
                    id_simul = get_simulation_seed(parameters, i_simul)
                    simulation_results[i_simul] = fun(list_parameters_set[i], id_simul)
                    next!(p)
                end

                # Concatenate results from all threads
                i_res = vcat(simulation_results...)
                if n > 1
                    list_res[i] = hcat(repeat(sweep_df[i:i, :], inner = nrow(i_res)), i_res)
                else
                    list_res[i] = i_res
                end
            end
        end

        ## Then either write the whole or give it back
        if parameters[:write_file]
            #-> Concatenate and save
            CSV.write(_build_filepath(parameters[:name_model], parameters, parameters[:parameters_to_omit], ".csv"; swept=sweep), vcat(list_res...))
            return nothing
        elseif parameters[:split_sweep]
            return list_res
        else
            return vcat(list_res...)
        end
    end

    if parameters[:write_file]
        if (parameters[:split_sweep] && parameters[:split_simul]) || (parameters[:split_simul] && isempty(sweep))
            # -> Parallelise both: each sweep point, each replicate, independently saved
            jobs = collect(Iterators.product(1:n, 1:parameters[:n_simul]))
            if parameters[:distributed]
                @sync @distributed for k in eachindex(jobs)
                        i, i_simul = jobs[k]
                        id_simul = get_simulation_seed(parameters, i_simul)
                        res = fun(list_parameters_set[i], id_simul)
                        CSV.write(_build_filepath(parameters[:name_model], list_parameters_set[i], parameters[:parameters_to_omit], ".csv", id_simul), res)
                end
            else
                Threads.@threads for k in eachindex(jobs)
                        i, i_simul = jobs[k]
                        id_simul = get_simulation_seed(parameters, i_simul)
                        res = fun(list_parameters_set[i], id_simul)
                        CSV.write(_build_filepath(parameters[:name_model], list_parameters_set[i], parameters[:parameters_to_omit], ".csv", id_simul), res)
                        next!(p)
                end
            end
        elseif parameters[:split_sweep]
            # -> Parallelise only over parameter sets
            if parameters[:distributed]
                @sync @distributed for i in 1:n
                    filepath = _build_filepath(parameters[:name_model],list_parameters_set[i],parameters[:parameters_to_omit],".csv")
                    for i_simul in 1:parameters[:n_simul]
                        id_simul = get_simulation_seed(parameters, i_simul)
                        res = fun(list_parameters_set[i], id_simul)
                        CSV.write(filepath, res;
                                append = i_simul > 1,
                                writeheader = i_simul == 1)
                    end
                end
            else
                Threads.@threads for i in 1:n
                    filepath = _build_filepath(parameters[:name_model],list_parameters_set[i],parameters[:parameters_to_omit],".csv")
                    for i_simul in 1:parameters[:n_simul]
                        id_simul = get_simulation_seed(parameters, i_simul)
                        res = fun(list_parameters_set[i], id_simul)
                        CSV.write(filepath, res;
                                append = i_simul > 1,
                                writeheader = i_simul == 1)
                        next!(p)
                    end
                end
            end
        end


    return nothing


    end
end


function parallel_if(distributed::Bool, range, body)
    if distributed
        @sync @distributed for k in range
            body(k)
        end
    else
        Threads.@threads for k in range
            body(k)
        end
    end
end