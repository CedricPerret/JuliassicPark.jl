#-----------------------------------------------
#*** Generate population
#-----------------------------------------------

"""
    initialise_population(z_ini, n_ini, n_patch; boundaries = nothing, simplify = true, n_loci = 0)

Initialise a population.

Arguments
- z_ini: Source for trait values. Accepts:
    - Number (copied to all individuals)
    - AbstractVector (sampled; if length == 1 it is treated as a single individual to copy)
    - Distributions.Distribution (sampled; truncated if boundaries are given)
    - Tuple of per-trait sources, one entry per trait
    - AbstractDataFrame or any Tables.jl table with columns :gen, :patch, and :z or :z1, :z2, ...
      In this case the last generation is used. The table is authoritative and n_ini, boundaries, and n_loci are ignored.
- n_ini::Int: Number of individuals per patch when sampling or copying.
- n_patch::Int: Number of patches.
- boundaries: Per-trait bounds. For one trait pass [min, max] or (min, max). For multiple traits pass one [min, max] per trait.
- simplify::Bool: If true and n_patch == 1, return a single vector of individuals. Faster for downstream functions, but outputs have less uniform dimensions.
- n_loci::Int: Number of loci for diploids. If > 0, each individual stores an n_loci × 2 genotype per trait.

Returns
- population: A vector of patches, each a vector of individuals. If simplify == true and n_patch == 1, returns a single vector.
  Individuals are Float64 for one trait, or tuples for multiple traits. With n_loci > 0, entries are genotype arrays.

Examples

# One trait, scalar input
initialise_population(1.0, 5, 2; boundaries = [0.0, 2.0])

# Distribution input with bounds
initialise_population(Normal(0, 1), 5, 2; boundaries = [-1.0, 1.0])

# Draw randomly from a vector
initialise_population([1.0, 2.0, 3.0], 5, 2)

# Two traits, independent sources
initialise_population((Normal(0, 1), [0.5, 1.0, 1.5]), 5, 2; boundaries = [(-2.0, 2.0), (0.0, 2.0)])

# Diploid with two loci
initialise_population(Normal(0, 1), 5, 2; boundaries = [-1.0, 1.0], n_loci = 0)

#Generator per patch
initialise_population([Normal(1.,0.1),Normal(0.2,0.001)], 5, 2; boundaries = [0.0, 2.0])
# If you want to provide a single scalar per patch, encapsulate them
initialise_population([[1.0],[0.5]], 5, 2; boundaries = [0.0, 2.0])

"""
function initialise_population(z_ini::AbstractDataFrame, n_ini::Int, n_patch::Int; boundaries=nothing, simplify=true, n_loci=0)
    # This is better than guessing from the vector dimensions as if n_cst = false, we can't infer when the input needs to be sampled or taken as it is.
    population = get_pop_at_last_gen(z_ini)
    if simplify && n_patch == 1 && length(population) == 1
        population = population[1]
    end
    return population
end


function initialise_population(z_ini, n_ini, n_patch::Int; boundaries=nothing, simplify = true, n_loci = 0)
    n_trait = z_ini isa Tuple ? length(z_ini) : 1
    #--- We wrap the boundaries the time of the initialisation (cannot wrap the parameter before because with a single trait, population simplifies to scalar rather than tuple)
    ## We want to get one generator per patch per trait.
    z_ini_complete = z_ini isa Tuple ? collect(Any, z_ini) : Any[z_ini]
    #--- If only one boundary set is given, we either replicate or wrap it.
    ## - If boundaries are missing, fill with `nothing`.
    ## - If a single pair of limits is given, replicate it for all traits.
    if boundaries === nothing
        boundaries = fill(nothing, n_trait)
    elseif !isa(boundaries[1], Vector) && !isa(boundaries[1], Tuple)
        boundaries = [boundaries for _ in 1:n_trait]
    end

    ## Here, it would make sense and is more elegant to consider that asexual reproduction is one loci for haploid (direct mapping)
    ## However, generalizing this rule makes it ambiguous to differentiate asexual from single-locus cases.
    ## As a temporary fix, we handle this directly here.
    n_loci_temp = n_loci == 0 ? 1 : copy(n_loci)
    ploidy =  n_loci == 0 ? 1 : 2
 
    ## We aim to have one generator per trait per patch.
    ## Note that it is not possible to give a single Real per patch,
    ## as this conflicts with the logic used to sample from a random vector.
    ## This can still be done by providing a vector of one element for each patch.
    for i_trait in 1:n_trait
        if !(z_ini_complete[i_trait] isa AbstractVector) || eltype(z_ini_complete[i_trait]) <: Real
             #-> Case 1: single generator or a flat vector of values
            z_ini_complete[i_trait] = [z_ini_complete[i_trait] for _ in 1:n_patch]
        else
            #-> Case 2: one generator per patch, we just check if there are the right size
            @assert length(z_ini_complete[i_trait]) == n_patch "Length of generators given for trait $(i_trait) ($(length(z_ini_complete[i_trait]))) does not match number of patches ($(n_patch)). Provide either one generator for the entire population, or one per patch."
        end
    end

    #--- Create populations: for each trait and patch, generate n_ini individuals.
    populations = [[
        [_generate_trait(z_ini_complete[i_trait][i_patch];
            boundaries=boundaries[i_trait], n_loci=n_loci_temp, ploidy=ploidy) for i_ind in 1:n_ini] for i_patch in 1:n_patch]
        for i_trait in 1:n_trait]
    
    #--- Restructure the population obtained
    if n_trait > 1
        population = [Tuple.(invert(getindex.(populations, i))) for i in 1:n_patch]
    else
        #-> No tuple. Each individual is a scalar.
        population = populations[1]
    end

    ## Simplify if single patch.It makes it faster for other function. 
    (simplify && n_patch == 1) && (population = only(population))

    return population
end

function _generate_trait(z_ini::Real; boundaries, n_loci, ploidy)
    boundaries !== nothing && @assert (z_ini <= boundaries[2] && z_ini >= boundaries[1]) "Initial values out of bounds."
    unwrap(fill(z_ini, n_loci, ploidy))
end

function _generate_trait(z_ini::Distributions.Distribution; boundaries, n_loci, ploidy)
    ## If it is nothing, then we have unbounded distribution
    boundaries === nothing ? unwrap(rand(z_ini, n_loci, ploidy)) : unwrap(rand(truncated(z_ini, boundaries[1], boundaries[2]), n_loci, ploidy))
end

function _generate_trait(z_ini::AbstractVector; boundaries, n_loci, ploidy)
    boundaries !== nothing && @assert (all(z_ini .<= boundaries[2]) && all(z_ini .>= boundaries[1])) "Initial values out of bounds."
    unwrap(rand(z_ini, n_loci,ploidy))
end



# function get_pop_at_gen(res::DataFrame; gen=:last, patch_col::Symbol=:patch, trait_cols=:auto)
#     # 1) Pick the generation of interest
#     selected_gen = gen === :last ? maximum(res.gen) : gen
#     # 2) Keep only rows for that generation and sort by patch for deterministic output order
#     df_gen = sort(res[res.gen .== selected_gen, :], patch_col)
#     # 3) Determine which columns are traits
#     if trait_cols === :auto
#         trait_cols = filter(c -> startswith(string(c), "z"), names(df_gen))
#     end
#     # Sort traits by numeric suffix so z1,z2,...,z10 are in the right order
#     trait_cols = sort(trait_cols, by = c -> begin
#         m = match(r"\d+$", string(c))
#         m === nothing ? 1 : parse(Int, m.match)
#     end)
#     # 4) Group rows by patch
#     grouped_by_patch = groupby(df_gen, patch_col)
#     # 5) Build population: per patch, collect individuals
#     if length(trait_cols) == 1
#         # Single trait: individuals are scalars
#         zcol = trait_cols[1]
#         return [collect(group_df[!, zcol]) for group_df in grouped_by_patch]
#     else
#         # Multiple traits: individuals are tuples (z1,z2,...)
#         n_traits = length(trait_cols)
#         return [begin
#                     # Pull columns once for cache friendliness
#                     cols = [group_df[!, trait_cols[t]] for t in 1:n_traits]
#                     [ntuple(t -> cols[t][row_i], n_traits) for row_i in 1:nrow(group_df)]
#                 end for group_df in grouped_by_patch]
#     end
# end

function get_pop_at_gen(res::DataFrame; gen=:last, patch_col::Symbol=:patch, trait_cols=:auto)
    # 1) Pick the generation of interest
    selected_gen = gen === :last ? maximum(res.gen) : gen
    # 2) Keep only rows for that generation and sort by patch for deterministic output order
    df_gen = sort(res[res.gen .== selected_gen, :], patch_col)
    # 3) Determine which columns are traits
    if trait_cols === :auto
        trait_cols = filter(c -> startswith(string(c), "z"), names(df_gen))
    end
    # Sort traits by numeric suffix so z1,z2,...,z10 are in the right order
    trait_cols = sort(trait_cols, by = c -> begin
        m = match(r"\d+$", string(c))
        m === nothing ? 1 : parse(Int, m.match)
    end)
    # 4) Group rows by patch and build a lookup dictionary
    grouped_by_patch = Dict(g[1, patch_col] => g for g in groupby(df_gen, patch_col))
    all_patches = sort(unique(res[!, patch_col]))
    # 5) Build population: per patch, collect individuals
    if length(trait_cols) == 1
        zcol = trait_cols[1]
        # Determine element type (if any nonempty patch)
        eltype_T = let first_nonempty = findfirst(p -> haskey(grouped_by_patch, p) && nrow(grouped_by_patch[p]) > 0, all_patches)
            first_nonempty === nothing ? Any :
            eltype(grouped_by_patch[all_patches[first_nonempty]][!, zcol])
        end
        return [haskey(grouped_by_patch, p) ?
                    collect(grouped_by_patch[p][!, zcol]) :
                    eltype_T[] for p in all_patches]
    else
        n_traits = length(trait_cols)
        # Infer element type from first nonempty patch
        eltype_T = let first_nonempty = findfirst(p -> haskey(grouped_by_patch, p) && nrow(grouped_by_patch[p]) > 0, all_patches)
            first_nonempty === nothing ? Any :
            Tuple{map(t -> eltype(grouped_by_patch[all_patches[first_nonempty]][!, trait_cols[t]]), 1:n_traits)...}
        end
        return [haskey(grouped_by_patch, p) ? begin
                    g = grouped_by_patch[p]
                    cols = [g[!, trait_cols[t]] for t in 1:n_traits]
                    [ntuple(t -> cols[t][r], n_traits) for r in 1:nrow(g)]
                end : eltype_T[] for p in all_patches]
    end
end

#--- Convenience wrapper for the last generation
get_pop_at_last_gen(res::DataFrame) = get_pop_at_gen(res; gen=:last, patch_col=:patch, trait_cols=:auto)



#-----------------------------------------------
#*** Prepare fitness function
#-----------------------------------------------

"""
    preprocess_fitness_function(population, fitness_function, parameters)

Preprocess a fitness function to handle different levels of input: individual, group of individuals, or metapopulation.

This function adjusts the provided `fitness_function` to ensure it can handle different cases and always returns a consistent output format. 
It adjusts (i) how it is run e.g. fitness function applies to individual but population is a vector so broadcast + (ii) reorganise the output.
The goal is to obtain an output that is a vector of size [n_variable] containing a vector of size [n_patch] of vectors of size [n_pop] (in a metapopulation) or a vector of size [n_pop].
For instance, when the fitness function gives fitness w + 1 output o, applies to an individual and the population is a vector, we want to move from [w_1,o_1],[w_2,o_2] to [w_1,w_2], [o_1, o_2]
For instance, when the fitness function gives fitness w + 1 output o, applies to an individual and the population is a vector of vector, we want to move from [[w_11,o_11],[w_21,o_21],[[w_12,o_12],[w_22,o_22]] to [[w_11,w_21,w_12,w_22],  [o_11,o_21,o_12,o_22]]

Fitness function output can be an element, or a tuple. This function uses `collect` if the output is a tuple.

# Arguments
- `population::Vector` or `population::Vector{Vector}`: The population trait.
- `fitness_function::Function`: The fitness function to preprocess.
- `parameters`: Additional parameters to pass to the fitness function.

# Returns
- `Function`: A preprocessed fitness function that handles various input cases and returns a consistent output format.

# Examples

```
using DataFrames

julia> function gaussian_fitness_function(individual; optimal, sigma)
           fitness = exp(-((individual - optimal)^2) / (2 * sigma^2))
           secondary = individual + optimal
           return fitness, secondary
       end

julia> population = [[0.5, 0.2, 0.1], [1.5, 0.8, 0.95]]

julia> instanced_fitness_function = preprocess_fitness_function(population, gaussian_fitness_function, [:optimal => 1, :sigma => 2], identity)
(Function)

# Using the fitness function directly
julia> [gaussian_fitness_function.(group; optimal=1, sigma=2) for group in population]
2-element Vector{Vector{Tuple{Float64, Float64}}}:
 [(0.9692332344763441, 1.5), (0.9231163463866358, 1.2), (0.903707077873196, 1.1)]
 [(0.9692332344763441, 2.5), (0.9950124791926823, 1.8), (0.9996875488230391, 1.95)]

# Using the preprocessed fitness function
julia> instanced_fitness_function(population; optimal=1, sigma=2)
2-element Vector{Vector{Tuple{Float64, Float64}}}:
 [[0.9692332344763441, 0.9231163463866358, 0.903707077873196], [0.9692332344763441, 0.9950124791926823, 0.9996875488230391]]
 [[1.5, 1.2, 1.1], [2.5, 1.8, 1.95]]
```
"""
function preprocess_fitness_function(population, fitness_function,parameters, correction)
    ## For now, we changed the tuple output of fitness function to vector as it is easier to work with later. 
    #--- Define the instanced fitness function
    instanced_fitness_function = nothing
    if correction == 2
        ##If single output,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
        ##If multiple output,  invert so that we have a vector by output rather by individual [(o1_ind1,o2_ind1),(o1_ind2,o2_ind2)] => ([o1_ind1, o2_ind1], [o2_ind1, o2_ind2])
        ## Same logic at higher level.
        instanced_fitness_function = function(population; parameters...)
            #Past one
            #invert([_my_invert(collect.(ensure_tuple.(fitness_function.(group; parameters...)))) for group in population])
            invert([collect(invert(ensure_tuple.(fitness_function.(group; parameters...)))) for group in population])
        end
    elseif correction == 1
        ##If single output,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
        ##If multiple output,  invert so that we have a vector by output rather by individual [(o1_ind1,o2_ind1),(o1_ind2,o2_ind2)] => ([o1_ind1, o2_ind1], [o2_ind1, o2_ind2])
        instanced_fitness_function = function(population; parameters...)
            #_my_invert(collect.(ensure_tuple.(fitness_function.(population; parameters...))))
            collect(invert(ensure_tuple.(fitness_function.(population; parameters...))))
            #collect(_my_invert(ensure_tuple.(fitness_function.(population; parameters...))))
        end
    elseif correction == 0
        ##If single output,  vectorize if to put it as the only element of the vector output [o1_ind1,o1_ind2] => [[o1_ind1,o1_ind2]]
        instanced_fitness_function = function(population; parameters...)
            collect(ensure_tuple(fitness_function(population; parameters...)))
        end
    elseif correction == -1
        #-> In-place so correction is 0 but fitness function take fitness.
        #instanced_fitness_function  = fitness_function
        if fitness_function(population,vv(0.,population);parameters...) === nothing
            #-> Replace by empty vector so we can use the same function to feed to output after
            instanced_fitness_function = function(population, fitness; parameters...)
                (fitness_function(population, fitness; parameters...))
                return Any[]
            end
        else
            #-> There are extras. We need to collect in case they are tuple and vectorise in case there is one.
            instanced_fitness_function = function(population, fitness; parameters...)
                collect(ensure_tuple(fitness_function(population, fitness; parameters...)))
            end
        end
    else
        error("Could not infer the correct input type for fitness_function.")
    end

    return instanced_fitness_function

end

"""
    extract_output_names(population, fitness_function, parameters, correction)

Extracts the names of additional outputs from a fitness function if it returns a `NamedTuple`. Only works if the function returns a `NamedTuple`; otherwise returns an empty vector.
"""
function extract_output_names(population, fitness_function, parameters,correction)
    sample_input = correction == 2 ? population[1][1] :
                    correction == 1 ? population[1] :
                    population
    position_extras_output = 0
    if correction != -1
        output = fitness_function(sample_input; parameters...)
        #-> Skip fitness
        position_extras_output = 2
    elseif correction == -1
        #-> In-place 
        output = fitness_function(sample_input, vv(0.,sample_input); parameters...)
        position_extras_output = 1
    end
    if output isa NamedTuple
        return [string(k) for k in keys(output)[position_extras_output:end]]
    else
        return String[]
    end
end


## Make a single sample run of fitness function. Use to preprocess the fitness function (for instance the size of the output)
function peek_fitness_output(population, fitness_function, parameters, correction)
    sample_input = correction == 2 ? population[1][1] :
                    correction == 1 ? population[1] :
                    population
    if correction != -1
        output = fitness_function(sample_input; parameters...)
    elseif correction == -1
        #-> In-place 
        output = fitness_function(sample_input, similar(sample_input); parameters...)
    end
    return output
end

function _infer_fitness_function_correction(population,fitness_function::Function, parameters,genotype_to_phenotype_mapping)
    ##Do not correct for genotype_to_phenotype mapping to get type in signature 
    correction = _infer_fitness_function_correction_by_signature(population,fitness_function)
    if length(methods(fitness_function).ms) > 1
        @warn "Multiple methods defined for the fitness function; falling back to manual inference."
        correction = _infer_fitness_function_correction_by_trial_error(genotype_to_phenotype_mapping(population), fitness_function, parameters)
    elseif correction < 0
        @warn "Could not infer input level from method signature; falling back to manual inference."
        correction = _infer_fitness_function_correction_by_trial_error(genotype_to_phenotype_mapping(population), fitness_function, parameters)
    end
    return(correction)
end



function _infer_fitness_function_correction_by_signature(population,fitness_function::Function)
    ## See https://discourse.julialang.org/t/obtaining-parameter-types-for-a-function/62169/8
    arg_type = fieldtypes((methods(fitness_function).ms[1].sig))[2]
    if arg_type <: AbstractVector{<:AbstractVector}
        if typeof(population) <: Vector{<:Vector}
            return 0
        elseif typeof(population) <: Vector
            return error("Fitness function requires at least a metapopulation")
        end
    elseif arg_type <: AbstractVector
        if typeof(population) <: Vector{<:Vector}
            return 1
        elseif typeof(population) <: Vector
            return 0
        end
    elseif arg_type <: Number || arg_type <: Tuple || arg_type <: Matrix
        if typeof(population) <: Vector{<:Vector}
            return 2
        elseif typeof(population) <: Vector
            return 1
        end
    end
    return -1
end


function _infer_fitness_function_correction_by_trial_error(population, fitness_function, parameters)
    need_no_change = _test_if_function_works(() -> fitness_function(population; parameters...))
    need_to_be_distributed = _test_if_function_works(() -> fitness_function.(population; parameters...))
    #--- Only test individual level if population appears to be a metapopulation
    need_to_be_double_distributed = false
    
    if isa(population, Vector{<:Vector})
        need_to_be_double_distributed = _test_if_function_works(() -> fitness_function(population[1][1]; parameters...))
    end
    results = [need_no_change, need_to_be_distributed, need_to_be_double_distributed]

    if !any(results)
            error("""
            There is an error in the fitness function. To see the details of the error, 
                please explicitly annotate the input type in your fitness function:
                - ::Number for individuals
                - ::Vector for groups
                - ::Vector{<:Vector} for metapopulations
            """)
        elseif count(identity, results) > 1
            error("""
            Ambiguity: The fitness function works for multiple resolution of population.
            ➤ Please explicitly annotate the input type in your fitness function:
                - ::Number for individuals
                - ::Vector for groups
                - ::Vector{<:Vector} for metapopulations
            """)
        else
            return findfirst(results) - 1
        end
end

function _test_if_function_works(f::Function)
    try
        f()
        return true
    catch e
        return false
    end
end



## Test
# function fit_pop(population::Vector{Vector};kwargs...)
#     [group .+ 2 for group in population]
# end
# function fit_group(population::Vector;kwargs...)
#     population .+ 2
# end
# function fit_ind(population::Float64;kwargs...)
#     population + 2
# end
# function fit2(population;kwargs...)
#         population + 2
# end
# function fit_error(population;kwargs...)
#         sqrt(-2)
# end
# population=rand(100)
# infer_fitness_function_correction(population,fit_pop,parameters_example)
# infer_fitness_function_correction(population,fit_group,parameters_example)
# infer_fitness_function_correction(population,fit_ind,parameters_example)
# infer_fitness_function_correction(population,fit2,parameters_example)
# infer_fitness_function_correction(population,fit_error,parameters_example)

# population=[rand(100) for _ in 1:2]
# infer_fitness_function_correction(population,fit_pop,parameters_example)
# infer_fitness_function_correction(population,fit_group,parameters_example)
# infer_fitness_function_correction(population,fit_ind,parameters_example)
# infer_fitness_function_correction(population,fit2,parameters_example)


#-----------------------------------------------
#*** Seed
#-----------------------------------------------

function generate_random_seed()
    return rand(10^9:10^10)
end


function get_simulation_seed(parameters::Dict, i_simul::Int)
    if haskey(parameters, :seed)
        #-> We use the seed provided by the user
        if parameters[:n_simul] == 1
            return parameters[:seed]
        else
            @assert isa(parameters[:seed], AbstractVector) "If n_simul > 1, parameter :seed must be a vector of seeds"
            @assert length(parameters[:seed]) == parameters[:n_simul] "Expected one seed per simulation"
            return parameters[:seed][i_simul]
        end
    else
        #-> We generate it
        return generate_random_seed()
    end
end

#-----------------------------------------------
#*** Sweep
#-----------------------------------------------

"""
    get_parameters_from_sweep(base::Dict, sweep::Dict{Symbol,Vector};
                    mode::Symbol = :grid,
                    filter::Function = x -> true)

Generates a list of parameter dictionaries by sweeping over combinations of values.

# Arguments
- `base`: base parameter dictionary to start from.
- `sweep`: dictionary of parameters to sweep over. Each value must be a `Vector`.
- `filter`: an optional function that receives a `Dict` and returns `true` if the combination should be included.

# Notes
- In `parameters`, the key `:sweep_grid` controls the combination strategy:
  - `:grid` (default) generates the full Cartesian product of all values.
  - `:zip` aligns values by position (like `zip()` in Julia).
# Returns
A vector of `Dict{Symbol,Any}`, one per valid parameter combination.
"""
function get_parameters_from_sweep(base::Dict, sweep::Dict{Symbol, <:AbstractVector};
                         filter::Function = x -> true)
    if isempty(sweep)
        return [base], DataFrame()
    end
    sweep_keys = collect(keys(sweep)); sweep_values = collect(values(sweep)); param_dicts = []; sweep_records = [];

    combos = base[:sweep_grid] ?
        Iterators.product(sweep_values...) :
        zip(sweep_values...)

    for combo in combos
        override = Dict(zip(sweep_keys, combo))
        params = merge(deepcopy(base), override)
        if filter(params)
            push!(param_dicts, params)
            push!(sweep_records, override)
        end
    end
    # Automatically print sweep summary
    # if !isempty(sweep_records)
    #     sweep_df = DataFrame(sweep_records)
    #     colnames = names(sweep_df)
    #     varying_cols = colnames[.!my_allequal.(eachcol(sweep_df))]
    #     println("\nParameter sweep over:")
    #     println(sweep_df)
    # else
    #     println("No valid parameter combinations after filtering.")
    # end
    return param_dicts, DataFrame(sweep_records)
end

##Does not take in account that some parameters can be a vector. Need to wrap these parameters in Ref. It is not good practice to have such oparameters as they are hard to print in csv
## Need to ref these values
##Also split the strings
function dt_parameter_sweep(parameters)
    ## To avoid Iterators.product to loop over each character of a string
    parameters_values = [v isa String ? [v] : v for v in values(parameters)]
    rename!(DataFrame(Iterators.product(parameters_values...)),[keys(parameters)...])
end

function dt_parameter_sweep(parameters_names,parameters_values)
    ## To avoid Iterators.product to loop over each character of a string
    parameters_values = [v isa String ? [v] : v for v in parameters_values]
    rename!(DataFrame(Iterators.product(parameters_values...)),[parameters_names...])
end

#-----------------------------------------------
#*** Filter and restructure kwargs
#-----------------------------------------------

"""
    _filter_kwargs(kwargs::Dict, accepted::Vector{Symbol})

Returns a `NamedTuple` containing only the keys in `accepted` that are present in `kwargs`.

Useful for selectively passing keyword arguments to functions.
"""
function _filter_kwargs(kwargs::Dict, accepted::Vector{Symbol})
    return (; (k => kwargs[k] for k in accepted if haskey(kwargs, k))...)
end
"""
    _prepare_kwargs_multiple_traits(parameters::Dict{Symbol,Any},
                                    kwarg_names::Vector{Symbol},
                                    n_traits::Int)

Extract and structure mutation-related keyword arguments for one or several traits.

- If `n_traits == 1`, return a flat `NamedTuple` of the selected parameters.  
- If `n_traits > 1`, build a per-trait `NamedTuple` and wrap them under keys
  `:mutation_parameter_1`, `:mutation_parameter_2`, …  
  Keys whose value is `nothing` for a given trait are omitted.

# Arguments
- `parameters::Dict{Symbol,Any}`: Dictionary of model parameters.
- `kwarg_names::Vector{Symbol}`: Parameter names to extract.
- `n_traits::Int`: Number of traits.

# Returns
- `NamedTuple`:  
  * Single trait → flat `NamedTuple`.  
  * Multiple traits → `NamedTuple` of the form  
    `(mutation_parameter_1 = (...,), mutation_parameter_2 = (...,), …)`.

# Example
```julia
parameters = Dict(
    :mu_m    => [0.001, 0.01],
    :sigma_m => [0.1, 0.2],
    :bias_m  => [0.0, nothing]
)
kwarg_names = [:sigma_m, :bias_m]

_prepare_kwargs_multiple_traits(parameters, kwarg_names, 2)
# Returns:
# (
#   mutation_parameter_1 = (sigma_m = 0.1, bias_m = 0.0),
#   mutation_parameter_2 = (sigma_m = 0.2,)
# )

"""
function _prepare_kwargs_multiple_traits(parameters, kwarg_names::Vector{Symbol}, n_traits)
    kwargs = _filter_kwargs(parameters, kwarg_names)
    if n_traits == 1
        return kwargs
    else
        per_trait = [ (; (k => v[i] for (k,v) in zip(keys(kwargs), values(kwargs)) if v[i] !== nothing)...) 
                      for i in 1:n_traits ]
        return (; (Symbol("mutation_parameter_$(i)") => per_trait[i] for i in 1:n_traits)...)
    end
end


"""
    _normalise_trait_parameters!(par::Dict{Symbol,Any}, names::Vector{Symbol}, n_traits::Int) -> Dict

Coerces shapes of trait-related parameters:
- If n_traits == 1: unwrap length-1 containers to scalars (except `:boundaries`, which may be a pair).
- If n_traits > 1: replicate scalars or length-1 containers across traits; assert correct lengths.
Special case: a single `[min, max]` (or `(min, max)`) for `:boundaries` is replicated per trait.
"""
function _normalise_trait_parameters!(parameters,kwargs_names, n_traits)
    for kwarg_i in kwargs_names
        haskey(parameters, kwarg_i) || continue
        parameters[kwarg_i] === nothing && continue
        if n_traits == 1
            if parameters[kwarg_i] isa AbstractVector
                if length(parameters[kwarg_i]) == 1
                    #-> To guard against boundaries being a vector
                    parameters[kwarg_i] = parameters[kwarg_i][1]
                end
            end
        elseif n_traits > 1
            if kwarg_i == :boundaries
                #-> A single vector of boundary has been provided
                if !isa(parameters[kwarg_i][1],Vector) && !isa(parameters[kwarg_i][1],Tuple)
                    parameters[kwarg_i]= fill(parameters[kwarg_i], n_traits)
                end
            elseif !(parameters[kwarg_i] isa AbstractVector) && !(parameters[kwarg_i] isa Tuple)
                #-> A scalar has been provided
                parameters[kwarg_i] = fill(parameters[kwarg_i], n_traits)
            else
                @assert length(parameters[kwarg_i]) == n_traits "$kwarg_i must be scalar or a vector of length $n_traits"
            end
        end
    end
end


#-----------------------------------------------
#*** Default parameters
#-----------------------------------------------

# --- Store mutable global defaults
const _DEFAULT_PARAMETERS = Ref{Union{Nothing, Dict{Symbol,Any}}}(nothing)

# --- Base (immutable) defaults
function _base_default_parameters()
    Dict(
        :n_gen => 1000,
        :n_patch => 1,
        :n_ini => 100,
        :n_loci => 0,
        :mu_m => 0.0001,
        :str_selection => 1.0,
        :n_print => 1,
        :j_print => 1,
        :de => 'i',
        :other_output_names => String[],
        :write_file => false,
        :name_model => "model_",
        :parameters_to_omit => String[],
        :additional_parameters_to_omit => Symbol[],
        :split_simul => false,
        :sweep_grid => true,
        :split_sweep => false,
        :distributed => false,
        :n_simul => 1,
        :simplify => true,
        :boundaries => nothing
    )
end

# --- Public accessor: returns current default parameters
function get_default_parameters()
    _DEFAULT_PARAMETERS[] === nothing ? _base_default_parameters() : deepcopy(_DEFAULT_PARAMETERS[])
end

function print_parameters(parameters;io::IO=stdout)
    parameters = parameters isa NamedTuple ? Dict{Any,Any}(pairs(parameters)) : parameters
    parameters = merge(get_default_parameters(), parameters)
    printed = Set{Symbol}()

    println(io, "Parameters\n" * "="^35)

    for (category, keys) in PARAMETER_CATEGORIES
        entries = [(k, parameters[k]) for k in keys if haskey(parameters, k)]
        if !isempty(entries)
            println(io, "\n[$category]")
            for (k, v) in entries
                desc = get(PARAMETER_DESCRIPTIONS, k, "")
                println(io, rpad(string(k), 25), "= ", rpad(string(v), 10), "   ", desc)
                push!(printed, k)
            end
        end
    end

    # Print any uncategorized parameters
    uncategorized = setdiff(keys(parameters), printed)
    if !isempty(uncategorized)
        println(io, "\n[Other parameters]")
        for k in sort(collect(uncategorized))
            println(io, rpad(string(k), 25), "= ", parameters[k])
        end
    end
end

function print_default_parameters(io::IO=stdout)
    print_parameters(get_default_parameters();io)
end



const PARAMETER_CATEGORIES = [
    "Demographic settings"    => [:n_gen, :n_ini, :n_patch, :n_loci],
    "Mutation settings"       => [:mu_m, :sigma_m, :bias_m, :boundaries, :mutation_type],
    "Reproduction settings"   => [:str_selection],
    "Output options"          => [:n_print, :j_print, :de, :other_output_names],
    "File writing"            => [:write_file, :name_model, :parameters_to_omit, :additional_parameters_to_omit],
    "Simulation control"      => [:n_simul, :split_simul, :sweep_grid, :split_sweep, :distributed, :simplify],
]

const PARAMETER_DESCRIPTIONS = Dict(
    :n_gen        => "Number of generations",
    :n_ini        => "Initial individuals per patch",
    :n_patch      => "Number of patches (groups)",
    :n_loci       => "Number of loci (for diploid traits)",
    :mu_m         => "Mutation rate per trait",
    :sigma_m      => "Mutation effect (standard deviation)",
    :bias_m       => "Bias in mutation (directional)",
    :boundaries   => "Trait bounds (min, max)",
    :mutation_type => "Mutation kernel type (optional)",
    :str_selection => "Strength of selection",
    :n_print      => "First generation to record output",
    :j_print      => "Interval between outputs",
    :de           => "Data resolution: 'g', 'p', or 'i'",
    :other_output_names => "Names for additional variables",
    :write_file   => "Whether to write results to disk",
    :name_model   => "Prefix for output filename",
    :parameters_to_omit => "Parameters excluded from filename",
    :additional_parameters_to_omit => "Additional derived parameters to exclude from output",
    :n_simul      => "Number of independent simulations",
    :split_simul  => "Whether to save each simulation replicate to a separate file. Requires :split_sweep = true. Also controls whether simulation replicates can be parallelised independently.",
    :sweep_grid => "Whether to use a full Cartesian product (`true`, default) or zip mode (`false`)",
    :split_sweep  => "Whether to save each parameter set to a separate file. Also controls whether parameter sets can be parallelised independently.",
    :distributed => "Whether to run simulations on distributed workets. (requires @everywhere for functions and imports)",
    :simplify     => "Flatten population if there is a single patch"
)

const PARAMETERS_TO_ALWAYS_OMIT = ["other_output_names",
"boundaries",
"n_cst",
"distributed",
"parameters_to_omit",
"name_model",
"additional_parameters_to_omit",
"write_file",
"split_simul",
"sweep_grid",
"split_sweep",
"simplify"]


# --- Public setter: globally override some defaults
function set_default_parameters!(overrides; strict=true)
    overrides = overrides isa NamedTuple ? Dict{Any,Any}(pairs(overrides)) : overrides
    current = _DEFAULT_PARAMETERS[] === nothing ? _base_default_parameters() : _DEFAULT_PARAMETERS[]
    _DEFAULT_PARAMETERS[] = merge(current, overrides)
end

# --- Reset to original defaults
function reset_default_parameters!()
    _DEFAULT_PARAMETERS[] = nothing
end

#-----------------------------------------------
#*** Additional parameters computed at the start of the simulation
#-----------------------------------------------

"""
    compute_derived_parameters!(parameters, additional_parameters; additional_parameters_to_omit)

Compute and insert derived parameters into the `parameters` dictionary in-place.  
Each derived parameter is defined by a user-supplied function that takes keyword arguments from `parameters`.

Typically used to generate constant values at the start of a simulation, such as carrying capacities, networks, or trait-specific constants.

# Arguments
- `parameters::Dict{Symbol, Any}`: Dictionary of model parameters. Modified in-place.
- `additional_parameters::Dict{Symbol, Function}`: Maps the name of each derived parameter to a function that returns its value. Functions must accept only keyword arguments, and those must be present in `parameters`.
- `additional_parameters_to_omit`: List of parameter names (as `Symbol` or `String`) to exclude from output saving.

# Returns
A tuple of:
- `parameters`: Updated dictionary including derived parameters.
- `cst_output_names::Vector{String}`: Names of derived parameters to be saved.
- `cst_output_values::Vector{Any}`: Corresponding values of those parameters.

# Requirements
- Each derived value must match the resolution of model output:
  - Single value (generation-level),
  - Vector of length `n_patch` (patch-level),
  - Vector of vectors with shape `(n_patch, n_ini)` (individual-level).
- If resolution is individual-level, `parameters[:n_cst]` must be `true`.

"""
function compute_derived_parameters!(parameters, additional_parameters; additional_parameters_to_omit)
    ## So that user can provide a vector of string or symbol
    additional_parameters_to_omit = Symbol.(additional_parameters_to_omit)
    cst_output_names=[]; cst_output_values=[];
    for (key, value) in additional_parameters
        if haskey(parameters, key)
            @warn "Derived parameter `:$key` already exists in the parameters dictionary. It will be replaced by the newly computed value."
        end
        value_parameters = value(; parameters...)
        #value_parameters = value_parameters isa AbstractVector ? value_parameters : [value_parameters]
        parameters[key] = value_parameters
        #--- Saving new parameters
        ##The new parameters are too long to print in the name of the output file but...
        ##they can be printed in the data output if specified by the users
        if first(string(key)) != '_' && key ∉ additional_parameters_to_omit
            # Validity checks
            @assert length(value_parameters) == 1 || length(value_parameters) == parameters[:n_patch] "The output of derived parameter `:$key` must have resolution compatible with the model: " * 
                "either one value per generation, per patch (length = n_patch), or per individual " *
                "(length = n_patch, with sub-vectors of length n_ini)."
            if isa(value_parameters[1], Vector)
                @assert length(value_parameters[1]) == parameters[:n_ini] "The output of derived parameter `:$key` must have resolution compatible with the model: " *
                    "either one value per generation, per patch (length = n_patch), or per individual " *
                    "(length = n_patch, with sub-vectors of length n_ini)."
                @assert parameters[:n_cst] == true "Derived parameter `:$key` is at individual-level resolution, " *
                    "but population size is not constant (`n_cst = false`)."
            end
            push!(cst_output_names, string(key))
            push!(cst_output_values, value_parameters)
        end
    end
    append!(parameters[:parameters_to_omit], Symbol.(keys(additional_parameters)))
    return parameters,cst_output_names,cst_output_values
end