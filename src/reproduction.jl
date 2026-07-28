#***********************************************
#*** Reproduction functions
#***********************************************


"""
    correct_fitness!(fitness)

Adds a small constant to all fitness values in-place to avoid issues when all individuals have zero fitness (relevant only in non-explicit reproduction functions). This preserves the logic that individuals with zero fitness are never chosen—unless all have zero.

# Note
- In-place modification is safe here because the true fitness values are already saved before reproduction is called (see `evol_model`).
"""
function correct_fitness!(fitness, str_selection::Float64)
    add_eps!(fitness)
    power!(fitness,str_selection)
    return nothing
end


#-----------------------------------------------
#*** ASEXUAL REPRODUCTION
#-----------------------------------------------

#*** Safe sampler
## - spot negative weights
## - copy mutable object (when genotype is a matrix) rather than copying the reference to avoid mutating multiple offspring at once

function handle_neg_fitness!(e)
    if isa(e, ArgumentError) && occursin("found negative weight", sprint(showerror, e))
        @error "Negative fitness value detected during reproduction."
        error("Fitness values must be strictly positive for Wright–Fisher reproduction.\n" *
            "You likely passed raw payoffs instead of valid fitness values.\n" *
            "Consider transforming payoffs into fitness using:\n" *
            "  • exponential mapping:    fitness = exp.(β .* payoff)\n" *
            "  • baseline shift:         fitness = payoff .- minimum(payoff) + ϵ\n" *
            "  • or any transformation ensuring fitness > 0.")
    else
        rethrow(e)
    end
end

function safe_sample(pop, w::AbstractWeights, n; replace::Bool = true)
    try
        return sample(pop, w, n; replace = replace)
    catch e
        handle_neg_fitness!(e)
    end
end

## When individuals are matrices, `sample` fills the new population with references to the same underlying matrices; 
## this can lead to multiple individuals pointing to the same object, so mutating one would also mutate the others.
function safe_sample(pop::AbstractVector{<:AbstractMatrix}, w::StatsBase.AbstractWeights, n; replace::Bool = true)
    new_pop = try
        new_pop = sample(pop, w, n; replace = replace) 
    catch e
        handle_neg_fitness!(e)
    end
    @inbounds for i in eachindex(new_pop)
        new_pop[i] = copy(new_pop[i])       # break aliasing so later mutation! does not touch parents
    end
    return new_pop
end



function safe_sample!(pop, w::StatsBase.AbstractWeights, dest; replace = true)
    try
        return sample!(pop, w, dest; replace = replace)
    catch e
        handle_neg_fitness!(e) 
    end
end

## When individuals are matrices, `sample!` fills the new population with references to the same underlying matrices; 
## this can lead to multiple individuals pointing to the same object, so mutating one would also mutate the others.
function safe_sample!(pop::AbstractVector{<:AbstractMatrix}, w::StatsBase.AbstractWeights, dest::AbstractVector{<:AbstractMatrix}; replace::Bool = true)
    try
        sample!(pop, w, dest; replace = replace)   # fills dest with references
    catch e
        handle_neg_fitness!(e)
    end
    @inbounds for i in eachindex(dest)
        dest[i] = copy(dest[i])       # break aliasing so later mutation! does not touch parents
    end
    return dest
end


#-----------------------------------------------
#*** Reproduction Moran Process
#-----------------------------------------------

"""
    reproduction_Moran_DB!(pop, fitness, str_selection, mu_m, mut_kwargs; n_replacement=1)

Simulates **death–birth Moran reproduction**:

1. **Random death**: `n_replacement` individuals are selected uniformly at random to be removed.
2. **Fitness-based birth**: New individuals are drawn with probabilities proportional to fitness^`str_selection` and mutated.

!!! note
    - This function modifies `pop` in place.
    - Individuals can replace themselves.

# Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `fitness::Vector{Float64}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m`: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.
- `n_replacement::Int`: Number of individuals to replace (default: 1).

# Returns
- `Vector{<:Any}`: Updated population after reproduction and mutation.

# Example
```julia
pop = [0.1, 0.2, 0.3, 0.4]
fitness = [1.0, 2.0, 3.0, 4.0]
str_selection = 1.0
mu_m = 0.1
mut_kwargs = (; sigma_m=0.05, boundaries=(0.0, 1.0))

reproduction_Moran_DB!(pop, fitness, str_selection, mu_m, mut_kwargs)
```
"""
function reproduction_Moran_DB!(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; n_replacement = 1, kwargs...) where T
    correct_fitness!(fitness,str_selection)
    pop[sample(eachindex(pop),n_replacement,replace=false)] = mutation.(safe_sample(pop,Weights(fitness),n_replacement),Ref(mu_m);mut_kwargs...)
    return nothing
end


"""
    reproduction_Moran_BD!(pop, fitness, str_selection, mu_m, mut_kwargs; n_replacement=1)

Simulates **birth–death Moran reproduction**:

1. **Fitness-based death**: `n_replacement` individuals are selected inversely proportional to fitness^`str_selection` to be removed.
2. **Random birth**: New individuals are drawn randomly and mutated.

!!! note
    - This function modifies `pop` in place.
    - Individuals can replace themselves.

# Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `fitness::Vector{Float64}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m`: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.
- `n_replacement::Int`: Number of individuals to replace (default: 1).

# Returns
- `Vector{<:Any}`: Updated population after reproduction and mutation.

pop = [0.1, 0.2, 0.3, 0.4]
fitness = [1.0, 2.0, 3.0, 4.0]
str_selection = 1.0
mu_m = 0.1
mut_kwargs = (; sigma_m=0.05, boundaries=(0.0, 1.0))

reproduction_Moran_BD!(pop, fitness, str_selection, mu_m, mut_kwargs)
"""
function reproduction_Moran_BD!(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; n_replacement = 1, kwargs...) where T
    correct_fitness!(fitness,str_selection)
    pop[safe_sample(eachindex(pop),Weights(1 ./ (fitness)),n_replacement,replace=false)] = mutation.(sample(pop,n_replacement),Ref(mu_m);mut_kwargs...)
    return nothing
end

"""
    reproduction_Moran_pairwise_learning!(pop::Vector{T}, fitness::Vector{Float64}, str_selection::Float64, mu_m::Float64, mut_kwargs; kwargs...) where T

Applies **pairwise comparison learning** (also called imitation dynamics) to a population, modifying it in place. 

In each learning step, a learner and a teacher are randomly sampled. The learner adopts the teacher’s trait
with probability given by the Fermi update rule:
```math
P(\text{adopt}) = \\frac{1}{1 + \\exp(-β(f_t - f_l))}
```
where f_t and f_l are the fitness of the teacher and learner, and β = str_selection is the selection strength.

This process is widely used in models of finite population dynamics (e.g., Traulsen et al. 2016).

## Arguments
- `pop::Vector{<:Any}`: Trait values for each individual.
- `fitness::Vector{Float64}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m``: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.
- `kwargs`: Additional keyword arguments (unused but accepted for interface compatibility).

## Example
```julia
using StatsBase

pop = [0.1, 0.2, 0.3, 0.4]
fitness = [1.0, 2.0, 3.0, 4.0]
str_selection = 1.0
mu_m = 0.1
mut_kwargs = (; sigma_m=0.05, boundaries=(0.0, 1.0))

reproduction_Moran_pairwise_learning!(pop, fitness, str_selection, mu_m, mut_kwargs)
```
"""
function reproduction_Moran_pairwise_learning!(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; kwargs...) where T
    learner_index, teacher_index = sample(1:length(pop),2,replace=false)
    if rand() < 1/(1+exp(-str_selection*(fitness[teacher_index] - fitness[learner_index])))
        pop[learner_index] = mutation(pop[teacher_index],mu_m;mut_kwargs...)
    end
    return nothing
end

#-----------------------------------------------
#*** Reproduction Wright–Fisher process
#-----------------------------------------------

"""
    reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)

Simulates **Wright–Fisher reproduction**, with optional support for metapopulations.

- If `pop` is a flat vector, reproduction occurs in a single unstructured population.
- If `pop` is a vector of vectors (i.e. a metapopulation), all individuals are pooled for reproduction
  and then redistributed randomly into patches of the same size (i.e. well-mixed).

1. Each offspring is drawn by sampling a parent with probability proportional to `fitness^str_selection`.
2. Offspring are mutated.
3. Previous population perishes.

!!! note
    - Population size remains constant.
    - Reproduction is asexual. See `reproduction_WF_sexual` for sexual reproduction.

# Arguments
- `pop::Vector{<:Any}` or `Vector{Vector{<:Any}}`: Trait values for each individual.
- `fitness::Vector{Float64}` or `Vector{Vector{Float64}}`: Fitness values corresponding to each individual in `pop`.
- `str_selection::Float64`: Scaling exponent for selection strength. Determines how sharply fitness differences affect reproduction.
- `mu_m`: mutation_probability
- `mut_kwargs`: Further arguments used for mutation.

# Returns
- Same structure as `pop`: a new population or metapopulation after reproduction.

# Examples
```julia
pop = [1.0, 2.0, 3.0, 4.0]
fitness = [0.1, 0.2, 0.3, 10.0]
str_selection = 1.0
mu_m = 0.2
mut_kwargs = (; sigma_m = 1.0, boundaries = (0.0, 5.0))
new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)

# Well-mixed metapopulation
using BenchmarkTools
pop = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
fitness = [[0.1, 0.2], [0.3, 10.0], [0.2, 0.1]]
new_pop = reproduction_WF(pop, fitness, str_selection, mu_m, mut_kwargs)
@btime reproduction_WF!(pop,new_pop,  fitness, str_selection, mu_m, mut_kwargs)
```
"""
#--- Wright–Fisher: single population (vector)
function reproduction_WF(pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m, mut_kwargs; kwargs...) where T
    #@ Equivalent to Moran death–birth reproduction with n_replacement = population size,
    #@ but avoids explicitly sampling deaths, which improves performance.  
    correct_fitness!(fitness,str_selection)
    offspring_pop = safe_sample(pop,Weights(fitness),length(pop))
    mutation!(offspring_pop,mu_m;mut_kwargs...)
    return(offspring_pop)
end

#--- Wright–Fisher: metapopulation (vector of vectors)
function reproduction_WF(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...) where T
    #@ Reduce faster than Iterators.flatten
    #---1. Pooled all individuals into a single population and reproduce
    metapopulation=reproduction_WF(reduce(vcat,pop),reduce(vcat,fitness),str_selection,mu_m,mut_kwargs;kwargs...)
    #---2. Redistribute the pooled population into patches
    #@ Iterators.partition faster than own partition. 
    population = [collect(part) for part in Iterators.partition(metapopulation, length(pop[1]))]
end

#--- Wright–Fisher: in-place and single population (vector)
function reproduction_WF!(pop::Vector{T},new_pop::Vector{T},fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...) where T
    correct_fitness!(fitness,str_selection)
    safe_sample!(pop,Weights(fitness),new_pop)
    mutation!(new_pop,mu_m;mut_kwargs...)
end

#--- Wright–Fisher: in-place and metapopulation (vector of vectors)
function reproduction_WF!(pop::Vector{Vector{T}}, new_pop::Vector{Vector{T}}, fitness::Vector{Vector{Float64}}, str_selection::Float64, mu_m, mut_kwargs; kwargs...) where {T}
    n_groups = length(pop);group_size = length(pop[1]);
    #--- Correct fitness
    correct_fitness!(fitness,str_selection)
    #--- Generate global pool of index of parents
    parents_idx = safe_sample(1:(group_size * n_groups), StatsBase.Weights(reduce(vcat,fitness)), (group_size * n_groups))
    #--- Reassign the offspring
    for j in 1:n_groups
        for i in 1:length(pop[j])
            #--- Get the parent (need to transform i and j in a global index)
            parent_id  = parents_idx[(j - 1) * group_size + i]        
            #--- Assign (need to get the value of the parent using the global index to pair of index)                
            new_pop[j][i] = pop[(parent_id - 1) ÷ group_size + 1   ][(parent_id - 1) % group_size + 1  ]
        end
    end
    mutation!(new_pop,mu_m;mut_kwargs...)
    return nothing
end


#--- Internal: Wright–Fisher reproduction independently in each patch
#! Not exported — used when reproduction occurs locally within each subpopulation
"""
    _reproduction_WF_patchwise(pop, fitness, str_selection, mu_m, mut_kwargs)

Applies **Wright–Fisher reproduction independently to each patch** in a structured population.

Each group in `pop` is treated as an isolated population:
1. Offspring are sampled within the patch with probabilities ∝ `fitness^str_selection`
2. Offspring are mutated.
(3. Groups are not mixed.)

Use this function if reproduction is local and structure is maintained across generations.
To simulate **global reproduction** across a metapopulation, use `reproduction_WF` instead.

# Arguments
- `pop::Vector{Vector{<:Any}}`: Trait values for each individual in each group.
- `fitness::Vector{Vector{Float64}}`: Fitness values for each individual in each group.
- `str_selection::Float64`: Selection strength.
- `mu_m`: Mutation parameter.
- `mut_kwargs`: Keyword arguments passed to `mutation!`.

# Returns
- `Vector{Vector{<:Any}}`: New structured population with same grouping.

@internal
"""
function _reproduction_WF_patchwise(pop::Vector{Vector{T}}, fitness::Vector{Vector{Float64}},
                                    str_selection::Float64, mu_m, mut_kwargs; kwargs...) where T
    [reproduction_WF(group, fitness_group, str_selection, mu_m, mut_kwargs; kwargs...)
     for (group, fitness_group) in zip(pop, fitness)]
end

function _reproduction_WF_patchwise!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...) where T
    for i in eachindex(pop)
        reproduction_WF!(pop[i], new_pop[i], fitness[i], str_selection, mu_m, mut_kwargs; kwargs...)
    end
    return nothing
end


"""
    reproduction_WF_island_model_fixed_migrant_fraction(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate, kwargs...)

Simulates **Wright–Fisher reproduction with a fixed expected migrant fraction** in a structured population.

For each offspring position in a focal group:

1. With probability `1 - mig_rate`, the offspring is philopatric and its parent is sampled from the focal group with probability proportional to `fitness^str_selection`.
2. With probability `mig_rate`, the offspring is a migrant and its parent is sampled from the pooled population of all other groups with probability proportional to `fitness^str_selection`.
3. Mutation is applied to the resulting offspring.

This procedure fixes the expected proportion of migrant offspring at `mig_rate`. Parents from different nonlocal groups compete with one another for migrant positions, but local and nonlocal parents do not compete for the same offspring position.

!!! note
    - Total population size and group sizes remain constant.
    - `mig_rate` is the probability that an offspring position is filled by a migrant.
    - Migration is symmetric across groups.
    - Reproduction is asexual.
    - At least two groups are required when `mig_rate > 0`.

# Arguments
- `pop::Vector{Vector{T}}`: Trait values for individuals in each group.
- `fitness::Vector{Vector{Float64}}`: Corresponding fitness values.
- `str_selection::Float64`: Selection strength, used as an exponent applied to fitness.
- `mu_m`: Mutation probability.
- `mut_kwargs`: Keyword arguments passed to `mutation!`.

# Keyword Arguments
- `mig_rate::Float64`: Probability that an offspring position is filled by a migrant.
- `kwargs...`: Additional keyword arguments accepted for interface compatibility.

# Returns
- `Vector{Vector{T}}`: New structured population after reproduction, migration, and mutation, with the same grouping.

# Example
```julia
pop = [fill(i, 100) .+ 0.01 * collect(1:100) for i in 1:5]
fitness = [rand(100) for _ in 1:5]
mu_m = 0.1
str_selection = 1.0
mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 5.0))

new_pop = reproduction_WF_island_model_fixed_migrant_fraction(
    pop, fitness, str_selection, mu_m, mut_kwargs;
    mig_rate = 0.1
)
```
"""

function reproduction_WF_island_model_fixed_migrant_fraction(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    group_size = length(pop[1]) ; n_groups = length(pop);
    correct_fitness!(fitness,str_selection)
    migrants_flag=[rand(group_size) .< mig_rate for i in 1:n_groups]
    new_pop = [Vector{T}(undef, group_size) for i in 1:n_groups]
    for i in 1:n_groups
        #@ It is weirdly faster to flatten even if we have many small patches, than to draw an index and then map it on the vector of vector (see dev)
        new_pop[i][migrants_flag[i]] .= safe_sample(vcat_except(pop,i), StatsBase.Weights(vcat_except(fitness,i)), count(migrants_flag[i]))
        new_pop[i][.!(migrants_flag[i])] .= safe_sample(pop[i], StatsBase.Weights(fitness[i]), count(.!migrants_flag[i]))
        mutation!(new_pop[i],mu_m;mut_kwargs...)
    end
    return(new_pop)
end

function reproduction_WF_island_model_fixed_migrant_fraction!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    group_size = length(pop[1]) ; n_groups = length(pop);
    #--- Prep fitness
    correct_fitness!(fitness,str_selection)
    if group_size > 1
        #--- Generate philopatric population by reproducing everyone in their patches (non-philopatric will be replaced later)
        for j in eachindex(pop)
            safe_sample!(pop[j],Weights(fitness[j]),new_pop[j])
        end

        #--- Now, we draw the non-philopatric offspring
        ## Buffer
        offspring_of_migrants = falses(group_size)
        #--- For each group...
        for j in 1:n_groups
            #--- Draw who is an offspring_of_migrants
            for i in 1:group_size
                offspring_of_migrants[i] = rand() < mig_rate
            end
            #--- if none, skip the rest (no need to calculate the pop and fitness without focal group)
            if !any(offspring_of_migrants)
                #-> No immigrants in this group
                continue
            end
            #--- Otherwise build the pop and fitness without focal group
            idx = findall(offspring_of_migrants)  
            other_group   = vcat_except(pop, j)
            other_fitness = StatsBase.Weights(vcat_except(fitness, j))
            #--- Draw parents of migrants
            parents_of_migrants = safe_sample(other_group, other_fitness, length(idx))
            #--- Assign parents
            new_pop[j][idx] = parents_of_migrants
        end

    else
        #-> For groups of 1, we use another method to avoid generating each time the other_group and other_fitness
        #--- Generate philopatric population by copying everyone (non-philopatric will be replaced later)
        for j in 1:length(pop)
            new_pop[j][1] = copy(pop[j][1])
        end
        offspring_of_migrants = rand(length(pop)) .< mig_rate
        #--- We draw from the global pool and redraw if we get the same individual as parent (impossible since migration is true)
        #--- Draw offspring_of_migrants
        idx = findall(offspring_of_migrants)  
        ## Skip if no migrants
        if !isempty(idx)
            ## Generate weight fitness only once (need to flatten)
            w_fitness = StatsBase.Weights(getindex.(fitness, 1))
            #--- Draw parents of migrants
            #@ faster to generate all migrants then correct than do it one by one
            parents_of_migrants = safe_sample(1:length(pop), w_fitness, length(idx))
            ## Redraw only those that landed on their own patch, until clean
            while true
                bad = parents_of_migrants .== idx
                any(bad) || break
                parents_of_migrants[bad] = safe_sample(1:length(pop), w_fitness, count(bad))
            end
                #--- Assign parents
            for t in eachindex(idx)
                new_pop[idx[t]][1] = pop[parents_of_migrants[t]][1]
            end
        end
    end
    #--- Mutate
    mutation!(new_pop,mu_m;mut_kwargs...)
    return nothing
end

"""
    reproduction_WF_island_model_hard_selection!(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate, kwargs...)

Simulates **Wright–Fisher reproduction with hard selection and island-model migration** in a structured population.

All individuals compete for each breeding spot. The weight of each potential parent is given by its fitness, multiplied by a factor determined by its source group:
1. Individuals from the target group receive weight `(1 - mig_rate) * fitness^str_selection`.
2. Individuals from each other group receive weight `(mig_rate / (n_groups - 1)) * fitness^str_selection`.
3. The offspring of the target group are sampled from this single weighted parent pool.
4. Mutation is applied to the offspring.

This implements a **hard-selection regime** because local and nonlocal parents compete in the same lottery. Consequently, groups with greater total fitness can contribute more offspring to the next generation.

For efficiency, the population is flattened once. For each target group, a corresponding flat vector of parent weights is generated and used to sample a complete group of offspring.

!!! note
    - Total population size and group sizes remain constant.
    - All groups are assumed to have the same size.
    - Migration is symmetric across groups.
    - `mig_rate` weights the contribution of nonlocal parents; it is not necessarily the expected realised fraction of migrants.
    - Reproduction is asexual.

# Arguments
- `pop::Vector{Vector{T}}`: Trait values for individuals in each group.
- `fitness::Vector{Vector{Float64}}`: Corresponding fitness values.
- `str_selection::Float64`: Selection strength, used as an exponent applied to fitness.
- `mu_m`: Mutation probability.
- `mut_kwargs`: Keyword arguments passed to `mutation!`.

# Keyword Arguments
- `mig_rate::Float64`: Weight assigned to nonlocal parental contributions relative to local contributions.
- `kwargs...`: Additional keyword arguments accepted for interface compatibility.

# Returns
- `Vector{Vector{T}}`: New structured population after reproduction and mutation, with the same grouping.

# Example
```julia
pop = [fill(i, 100) .+ 0.01 * collect(1:100) for i in 1:5]
new_pop = [zeros(Float64,100) for i in 1:5]
fitness = [rand(100) for _ in 1:5]
mu_m = 0.1
str_selection = 1.0
mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 5.0))

reproduction_WF_island_model_hard_selection!(
    pop, new_pop, fitness, str_selection, mu_m, mut_kwargs;
    mig_rate = 0.1
)
```
"""
function reproduction_WF_island_model_hard_selection!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    group_size = length(pop[1]) ; n_groups = length(pop);
    #--- Special case: one group only
    if n_groups == 1
        reproduction_WF!(pop[1],new_pop[1],fitness[1],str_selection,mu_m,mut_kwargs;kwargs...)
        return nothing
    end

    #--- Prep fitness
    correct_fitness!(fitness,str_selection)

    #--- Flatten the population into a single global parent pool
    flat_pop = reduce(vcat,pop)

    #--- Buffer for the weighted fitness of all possible parents
    #@ The order corresponds to flat_pop: group 1, group 2, ..., group n_groups
    flat_fitness_with_mig = Vector{Float64}(undef, group_size * n_groups)

    #--- For each target group...
    for i in 1:n_groups
        #--- Build the parent weights for the target group
        for j in 1:n_groups
            #@ This form was empirically faster than assigning fitness and multiplying in place
            if i == j
                flat_fitness_with_mig[((j-1)*group_size + 1):(j*group_size)] .= fitness[j] .* (1 - mig_rate)
            else
                flat_fitness_with_mig[((j-1)*group_size + 1):(j*group_size)] .= fitness[j] .* (mig_rate / (n_groups - 1))
            end
        end
        #--- Draw all offspring from the single weighted parent pool
        safe_sample!(flat_pop, StatsBase.Weights(flat_fitness_with_mig), new_pop[i])
    end
    #--- Mutate
    mutation!(new_pop,mu_m;mut_kwargs...)
    return nothing
end

"""
    reproduction_WF_island_model_soft_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate, kwargs...)

Simulates Wright–Fisher reproduction with **soft selection** and island-model migration in a structured population.

Under soft selection, competition is **local**: for each offspring position in a focal group, the parent is chosen
1. from the same group with probability `1 - mig_rate`, with sampling probability proportional to fitness within that group; or
2. from a randomly chosen other group with probability `mig_rate`, with the parent then sampled proportionally to fitness within that source group.

This is a soft-selection scheme because selection is always evaluated **within groups**, and each group contributes the same number of offspring slots to the next generation.

For efficiency, this is implemented in a different way:
1. Each group first produces a full new generation independently using Wright–Fisher reproduction within the group.
2. Then, a subset of offspring positions (determined by `mig_rate`) is replaced by immigrants drawn from other groups.

This is usually faster as mig_rate tends to be low. 

!!! note
    - Maintains constant group sizes.
    - Migration is uniform and unidirectional (target → source group sampled independently).
    - Reproduction is asexual.

# Arguments
- `pop::Vector{Vector{T}}`: Trait values for each individual in each group.
- `fitness::Vector{Vector{Float64}}`: Corresponding fitness values.
- `str_selection::Float64`: Selection strength (exponent applied to fitness).
- `mu_m`: Mutation probability.
- `mut_kwargs`: Keyword arguments passed to `mutation`.

# Keyword Arguments
- `mig_rate::Float64`: Probability that an individual is replaced by a migrant.
- `kwargs...`: Additional arguments passed to `mutation`.

# Returns
- `Vector{Vector{T}}`: New structured population after local reproduction and migration.

# Example
```julia
pop = [fill(i, 100) .+ 0.01 * collect(1:100) for i in 1:5]
fitness = [rand(100) for _ in 1:5]
mu_m = 0.1
str_selection = 1.0
mut_kwargs = (; sigma_m = 0.5, boundaries = (0.0, 5.0))

new_pop = reproduction_WF_island_model_soft_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate = 0.1)

@btime new_pop = reproduction_WF_island_model_soft_selection(pop, fitness, str_selection, mu_m, mut_kwargs; mig_rate = 0.1)
```

"""
function reproduction_WF_island_model_soft_selection(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    n_groups = length(pop); groups_size = length.(pop)
    #--- Local reproduction
    ## This also corrects fitness and check for negative fitness
    new_pop=_reproduction_WF_patchwise(pop, fitness, str_selection, mu_m,mut_kwargs)
    #--- Number of migrants per group
    n_migrant_by_group =sum.([rand(groups_size[i]) .< mig_rate for i in 1:n_groups])
    #--- For each migrant, identify from which patch they come from.
    patch_of_parent=[random_int_except(1,n_groups,i,n_migrant_by_group[i]) for i in 1:n_groups]
    #--- Select parent of migrant based on fitness
    for i in 1:n_groups
        if n_migrant_by_group[i] == 0
            continue
        end
        new_pop[i][sample(eachindex(pop[i]),n_migrant_by_group[i],replace=false)] = mutation.(vcat(safe_sample.(pop[patch_of_parent[i]],StatsBase.Weights.(fitness[patch_of_parent[i]]),1)...),Ref(mu_m);mut_kwargs...)
    end
    return(new_pop)
end

function reproduction_WF_island_model_soft_selection!(pop::Vector{Vector{T}},new_pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m, mut_kwargs; mig_rate, kwargs...) where T
    n_groups = length(pop); groups_size = length.(pop)
    #--- Local reproduction
    ## This also corrects fitness and check for negative fitness
    _reproduction_WF_patchwise!(pop, new_pop, fitness, str_selection, mu_m,mut_kwargs)
    #--- Number of migrants per group
    n_migrant_by_group =sum.([rand(groups_size[i]) .< mig_rate for i in 1:n_groups])
    #--- For each migrant, identify from which patch they come from.
    patch_of_parent=[random_int_except(1,n_groups,i,n_migrant_by_group[i]) for i in 1:n_groups]
    #--- Select parent of migrant based on fitness
    for i in 1:n_groups
        if n_migrant_by_group[i] == 0
            continue
        end
        new_pop[i][sample(eachindex(pop[i]),n_migrant_by_group[i],replace=false)] = mutation.(vcat(safe_sample.(pop[patch_of_parent[i]],StatsBase.Weights.(fitness[patch_of_parent[i]]),1)...),Ref(mu_m);mut_kwargs...)
    end
    return(nothing)
end

"""
    reproduction_WF_copy_group_trait!(population, group_level_trait, fitness, str_selection, mu_m, mut_kwargs; group_fitness_fun, kwargs...)

Simulates **Wright–Fisher reproduction**, where individuals adopt the group-level trait (e.g. an institutional rule) of a successful group.

1. Each group’s fitness is computed using `group_fitness_fun`, applied to the individual fitness values.
2. For each group, a new set of individuals is created by sampling group-level traits according to group fitness (scaled by `str_selection`).
3. Mutation is applied in-place to each individual trait after copying.

This model is designed for scenarios where individuals *copy the trait in place* within a group, such as in cultural evolution models where group norms or strategies are transmitted collectively.

### Arguments
- `population::Vector{Vector{T}}`: The current structured population, grouped into patches.
- `group_level_trait::Vector{T}`: A single trait value associated with each group.
- `fitness::Vector{Vector{Float64}}`: Individual fitness values within each group.
- `str_selection::Float64`: Selection strength at the group level.
- `mu_m`: Mutation probability.
- `mut_kwargs`: Additional keyword arguments for `mutation!`.
- `group_fitness_fun`: Function mapping each group’s individual fitness values to a group-level fitness.
- `kwargs...`: Extra keyword arguments (ignored; for compatibility).

### Returns
- `Vector{Vector{T}}`: New structured population. Each individual copies the trait of a sampled group, followed by mutation.

### Example
```julia
group_fitness_fun = x -> mean(x)
pop = [[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]]
fitness = [[1000.0, 1.0, 1.0], [1.0, 1.0, 1.0]]
group_level_trait = mean.(pop)
mut_kwargs = (sigma_m = 0.1, boundaries = (0.0, 1.0))
new_pop = reproduction_WF_copy_group_trait(pop, group_level_trait, fitness, 1.0, 0.1, mut_kwargs; group_fitness_fun = group_fitness_fun)
```
"""
function reproduction_WF_copy_group_trait!(population::Vector{Vector{T}},new_population::Vector{Vector{T}},group_level_trait, fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs;group_fitness_fun,kwargs...) where T
    #--- Calculate group fitness and selection weights
    group_fitness = group_fitness_fun.(fitness)   
    correct_fitness!(group_fitness,str_selection)
    w_group_fitness = StatsBase.Weights(group_fitness)
    groups_sizes = length.(population)
    #--- Sample offspring per group into preallocated buffers
    @inbounds for j in 1:length(population)
        safe_sample!(group_level_trait, w_group_fitness, new_population[j])
    end
    #--- Mutate offspring in place
    mutation!(new_population, mu_m; mut_kwargs...)
    return nothing
end


#-----------------------------------------------
#*** Reproduction: explicit (individual-based, offspring explictly generated)
#-----------------------------------------------

"""
    reproduction_explicit_poisson(pop, fitness, str_selection, mu_m, mut_kwargs; kwargs...)

Simulates reproduction by generating an explicit number of offspring for each individual, drawn from a Poisson distribution centered on the `fitness^str_selection` of each individual.

!!! note
    - Reproduction is asexual.

### Arguments
- `pop::Vector{<:Any}`: Trait values for each individual in the population.
- `fitness::Vector{Float64}`: Fitness values for each individual.
- `str_selection::Float64`: Selection strength; exponent applied to fitness before reproduction.
- `mu_m::Float64`: Mutation probability.
- `mut_kwargs`: Additional keyword arguments passed to `mutation`.

### Keyword Arguments
- `kwargs...`: Ignored; included for interface compatibility.

### Returns
- `Vector{<:Any}`: New population vector containing all offspring after reproduction and mutation.

### Example
```julia
pop = [1.0, 2.0, 3.0, 4.0]
fitness = [0.1, 0.2, 0.3, 10.0]
str_selection = 1.0
mu_m = 0.2
mut_kwargs = (sigma_m=0.1, boundaries=(0.0, 1.0))

new_pop = reproduction_explicit_poisson(pop, fitness, str_selection, mu_m, mut_kwargs)
```
"""
#--- single population (vector)
function reproduction_explicit_poisson(pop::Vector{T}, fitness::Vector{Float64}, str_selection::Float64, mu_m, mut_kwargs; kwargs...) where T
    offspring = T[]
    for i in eachindex(pop)
        try
            n = rand(Poisson(fitness[i]^str_selection))
            for _ in 1:n
                push!(offspring, mutation(pop[i], mu_m; mut_kwargs...))
            end
        catch e
            if isa(e, DomainError) && occursin("Poisson: the condition λ >= zero(λ) is not satisfied", sprint(showerror, e))
                @error "Negative fitness value detected during reproduction." 
                error("Fitness values must be strictly positive for explicit reproduction.\n" *
                    "You likely passed raw payoffs instead of valid fitness values.\n" *
                    "Consider transforming payoffs into fitness using:\n" *
                    "  • exponential mapping:    fitness = exp.(β .* payoff)\n" *
                    "  • baseline shift:         fitness = payoff .- minimum(payoff) + ϵ\n" *
                    "  • or any transformation ensuring fitness > 0.")
            else
                rethrow(e)
            end
        end
    end
    return offspring
end


#--- Wrapper for metapopulation (vector of vectors)
function reproduction_explicit_poisson(pop::Vector{Vector{T}},fitness::Vector{Vector{Float64}},str_selection::Float64, mu_m, mut_kwargs; kwargs... ) where T
    #@ No significant differences if replacing by a loop
    [reproduction_explicit_poisson(group,fitness_group,str_selection,mu_m, mut_kwargs; kwargs... ) for (group,fitness_group) in zip(pop,fitness)]
end

#-----------------------------------------------
#*** SEXUAL REPRODUCTION
#-----------------------------------------------

#-----------------------------------------------
#*** Reproduction Wright–Fisher process
#-----------------------------------------------

"""
    reproduction_WF_sexual(pop, fitness, str_selection, mu_m, sigma_m, lower_bound, upper_bound)

Implements Wright–Fisher reproduction with sexual recombination across multiple loci.

Each individual is represented as a vector of diploid loci, where each locus is a `Matrix{Float64}` representing the allele pair. Offspring are generated by sampling two parents from the population, with probabilities proportional to their fitness raised to `str_selection`. Each offspring inherits one allele per locus from each parent, followed by potential mutation.

!!! note
    - The same parent may be paired with multiple others (no monogamy).
    - Each element of `pop` is a vector of matrices, each matrix representing a diploid locus.
    - The function assumes all individuals have the same number of loci.

### Arguments
- `pop::Vector{Vector{Matrix{Float64}}}`: The population, where each individual is a vector of diploid loci represented as `Matrix{Float64}`.
- `fitness::Vector{Vector{Float64}}`: Fitness values corresponding to individuals in `pop`. Each subvector contains the fitness of individuals in one subpopulation.
- `str_selection::Float64`: Strength of selection; fitness is raised to this power before sampling.
- `mu_m::Float64`: Mutation probability per allele.
- `sigma_m::Float64`: Standard deviation of mutations.
- `lower_bound::Float64`: Lower bound for allele values.
- `upper_bound::Float64`: Upper bound for allele values.

### Returns
- `Vector{Vector{Matrix{Float64}}}`: A new generation with potentially mutated offspring, structured identically to `pop`.

### Example
```julia
pop = [[rand(2, 2), rand(2, 2)], [rand(2, 2), rand(2, 2)]]
fitness = [[0.1, 0.2], [0.3, 10.0]]
str_selection = 1.0
mu_m = 0.2
mut_kwargs = (sigma_m=0.1, boundaries=(0.0, 1.0))
offspring = reproduction_WF_sexual(pop, fitness, 1.0, mu_m, mut_kwargs)
```
"""
#--- Wright–Fisher: single population (vector)
function reproduction_WF_sexual(pop,fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    correct_fitness!(fitness,str_selection)
    w_fitness = Weights(fitness)
    offspring = Vector{eltype(pop)}(undef, length(pop))
    for i in 1:length(fitness)
        #--- Sample two parents randomly (without replacement)
        parent1, parent2 = safe_sample(pop, w_fitness, 2, replace=false)
        offspring[i]=mutation(recombination(parent1,parent2),mu_m;mut_kwargs...)
    end
    return(offspring)
end

#@ We use in-place mutation for genotype represented as matrix to gain performance
function reproduction_WF_sexual(pop::Vector{Matrix},fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    correct_fitness!(fitness,str_selection)
    w_fitness = Weights(fitness)
    offspring = Vector{eltype(pop)}(undef, length(pop))
    for i in 1:length(fitness)
        #--- Sample two parents randomly (without replacement)
        parent1, parent2 = safe_sample(pop, w_fitness, 2, replace=false)
        offspring[i]=recombination(parent1,parent2)
        mutation!(offspring[i],mu_m;mut_kwargs...)
    end
    return(offspring)
end

#--- Wright–Fisher: metapopulation (vector of vectors)
function reproduction_WF_sexual(pop,fitness::Vector{Vector{Float64}},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    metapopulation=reproduction_WF_sexual(reduce(vcat,pop),reduce(vcat,fitness),str_selection,mu_m,mut_kwargs;kwargs...)
    population = [collect(part) for part in Iterators.partition(metapopulation, length(pop[1]))]
end

#--- Wright–Fisher: in-place and single population (vector)
function reproduction_WF_sexual!(pop, new_pop,fitness::Vector{Float64},str_selection::Float64,mu_m,mut_kwargs; kwargs...)
    correct_fitness!(fitness,str_selection)
    w_fitness = Weights(fitness)
    for i in 1:length(fitness)
        #--- Sample two parents randomly (without replacement)
        parent1, parent2 = safe_sample(pop, w_fitness, 2, replace=false)
        recombination!(new_pop[i],parent1,parent2)
    end
    #@ This is faster than doing in-place mutation for each genotype
    mutation!(new_pop,mu_m;mut_kwargs...)
    return(nothing)
end

#--- Wright–Fisher: in-place and metapopulation (vector of vectors)
function reproduction_WF_sexual!(pop, new_pop, fitness::Vector{Vector{Float64}}, str_selection::Float64, mu_m, mut_kwargs; kwargs...)
    group_size = length(pop[1]); n_groups   = length(pop)
    correct_fitness!(fitness, str_selection) 
    w_global = StatsBase.Weights(reduce(vcat, fitness)) 

    #--- For each child in each group, sample two global parents and recombine in place
    @inbounds for j in 1:n_groups
        for i in 1:group_size
            #--- sample two distinct parents from global pool
            parent1_id, parent2_id = safe_sample(1:(group_size * n_groups), w_global, 2; replace=false)
            #--- transform global index into group/within-group indices
            g1 = (parent1_id - 1) ÷ group_size + 1; i1 = (parent1_id - 1) % group_size + 1
            g2 = (parent2_id - 1) ÷ group_size + 1; i2 = (parent2_id - 1) % group_size + 1
            #--- Recombine into preallocated child
            recombination!(new_pop[j][i], pop[g1][i1], pop[g2][i2])
        end
    end
    mutation!(new_pop, mu_m; mut_kwargs...)
    return nothing
end


#--- Function to simulate free recombination
@inline recombination(parent1,parent2) = hcat(sample.(eachrow(parent1)),sample.(eachrow(parent2)))
#--- Wrapper if multiple traits
@inline recombination(parent1::Tuple,parent2::Tuple) = recombination.(parent1,parent2)
#--- In-place version
#@ Faster if we don't specify matrix of Float
@inline function recombination!(child::AbstractMatrix, parent1::AbstractMatrix, parent2::AbstractMatrix)
    for r in axes(child, 1)
        child[r, 1] = parent1[r, rand(axes(parent1, 2))]; child[r, 2] = parent2[r, rand(axes(parent2, 2))];
    end
    return nothing
end
#--- In-place version + tuple
@inline function recombination!(child::Tuple{Vararg{<:AbstractMatrix,N}},parent1::Tuple{Vararg{<:AbstractMatrix,N}},parent2::Tuple{Vararg{<:AbstractMatrix,N}}) where {N}
    for i in 1:N
        recombination!(child[i], parent1[i], parent2[i])
    end
    return nothing
end

#*** Programmatic registry of reproduction functions
# Pretty-printer driven by the registry
#*** Programmatic registry of reproduction functions
# Single source of truth for reproduction methods
# applies_to ∈ (:pop, :metapop, :both)
_REPRO_FUNS = [
    (name = :reproduction_Moran_DB!,                f = reproduction_Moran_DB!,                desc = "Death–birth Moran (in-place)",                       needs = Symbol[],                    applies_to = :pop),        # single population only
    (name = :reproduction_Moran_BD!,                f = reproduction_Moran_BD!,                desc = "Birth–death Moran (in-place)",                        needs = Symbol[],                    applies_to = :pop),        # single population only
    (name = :reproduction_Moran_pairwise_learning!, f = reproduction_Moran_pairwise_learning!, desc = "Pairwise comparison imitation (in-place)",            needs = Symbol[],                    applies_to = :pop),        # single population only

    (name = :reproduction_WF,                       f = reproduction_WF,                       desc = "Wright–Fisher reproduction",                          needs = Symbol[],                    applies_to = :both),       # has Vector and Vector{Vector} methods
    (name = :reproduction_WF!,                      f = reproduction_WF!,                      desc = "Wright–Fisher reproduction (in-place)",               needs = Symbol[],                    applies_to = :both),       # has Vector and Vector{Vector} methods

    (name = :reproduction_WF_island_model_fixed_migrant_fraction,  f = reproduction_WF_island_model_fixed_migrant_fraction,  desc = "Island model, fixed fraction, global competition", needs = [:mig_rate],              applies_to = :metapop),    # metapop only
    (name = :reproduction_WF_island_model_fixed_migrant_fraction!, f = reproduction_WF_island_model_fixed_migrant_fraction!, desc = "Island model, fixed fraction, global competition (in-place)", needs = [:mig_rate], applies_to = :metapop),    # metapop only
    (name = :reproduction_WF_island_model_hard_selection!, f = reproduction_WF_island_model_hard_selection!, desc = "Island model, hard selection, global competition (in-place)", needs = [:mig_rate], applies_to = :metapop),    # metapop only
    (name = :reproduction_WF_island_model_soft_selection,  f = reproduction_WF_island_model_soft_selection,  desc = "Island model, soft selection, local competition",  needs = [:mig_rate],              applies_to = :metapop),    # metapop only
    (name = :reproduction_WF_island_model_soft_selection!, f = reproduction_WF_island_model_soft_selection!, desc = "Island model, soft selection, local competition (in-place)",  needs = [:mig_rate], applies_to = :metapop),    # metapop only

    (name = :reproduction_explicit_poisson,         f = reproduction_explicit_poisson,         desc = "Explicit number of offspring (Poisson)",              needs = Symbol[],                    applies_to = :both),       # vector and wrapper for metapop

    (name = :reproduction_WF_sexual,                f = reproduction_WF_sexual,                desc = "Wright–Fisher, sexual recombination",                 needs = [:n_loci],                   applies_to = :both),       # single and metapop overloads
    (name = :reproduction_WF_sexual!,                f = reproduction_WF_sexual!,                desc = "Wright–Fisher, sexual recombination (in-place)",                 needs = [:n_loci],                   applies_to = :both),       # single and metapop overloads

    # Requires extra positional args not wired by evol_model. Keep listed for completeness.
    (name = :reproduction_WF_copy_group_trait!,      f = reproduction_WF_copy_group_trait!,      desc = "Copy group-level trait using group fitness",          needs = [:group_level_trait, :group_fitness_fun], applies_to = :metapop),
]

list_reproduction_functions() = _REPRO_FUNS

# Pretty-printer driven by the registry
function list_reproduction_methods(io::IO=stdout)
    println(io, "Available reproduction methods\n")
    w = maximum(length(string(x.name)) for x in _REPRO_FUNS)
    for x in _REPRO_FUNS
        scope = x.applies_to === :both ? "pop + metapop" :
                x.applies_to === :pop  ? "pop" :
                                         "metapop"
        n = string(x.name)
        pad = repeat(" ", w - length(n) + 2)
        println(io, "  ", n, pad, "- ", x.desc, "  [", scope, "]")
        if !isempty(x.needs)
            println(io, "     needs: ", join(x.needs, ", "))
        end
    end
end
