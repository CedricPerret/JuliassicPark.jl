#***********************************************
#*** Mutation function
#***********************************************

#--- The logic behind the dev
## To make the interface flexible and efficient, we use multiple dispatch for core logic and add `kwargs...` to allow passing extra arguments when needed.
## This allows for julia to specialise each function while being flexible

"""
    mutation(trait_value, mu; kwargs...)

Mutates a trait value `trait_value` with a mutation probability `mu`. The type of `trait_value` determines the kind of mutation applied:
- For `Float64`, mutation applies to a continuous trait and the mutated value is drawn according to the selected mutation type within the specified boundaries.
- For `Int`, mutation applies to a discrete trait and the mutated value is a random integer within the boundaries range.
- For `Bool`, mutation applies to a boolean trait and the mutated value is obtained by flipping the value.

### Arguments
- `trait_value`: Trait value (can be `Float64`, `Int`, or `Bool`).
- `mu::Float64`: Mutation probability.
- `kwargs...`: Additional keyword arguments depending on `trait_value`'s type.

### Keyword arguments (depending on trait type)

**For continuous (`Float64`) traits:**
- `mutation_type::Symbol`: Type of distribution to use (`:normal`, `:gumbel`, etc.).
- `sigma_m::Float64`: Standard deviation or scale.
- `bias_m::Float64` (optional): Directional bias (used by some mutation types).
- `boundaries::Tuple{Float64, Float64}`: Lower and upper bounds.

**For integer (`Int`) traits:**
- `boundaries::Tuple{Int, Int}`: Range for possible new values (excluding `trait_value`).

**For boolean (`Bool`) traits:**
- No additional arguments required.

### Returns
- The mutated or unchanged trait value, matching the input type.

### Examples
```julia
mutation(0.5, 0.1; mutation_type=:normal, sigma_m=0.1, boundaries=(0.0, 1.0))  # Float64
mutation(2, 0.2; boundaries=(0, 4))                                            # Int
mutation(true, 0.5)                                                            # Bool
mutation(0.5, 0.99; mutation_type=:normal_clamped, sigma_m=2, boundaries=(0.0, 1.0)) #Clamped
```
"""
#--- Fallback
function _trait_type_error(x, mu)
    error("No mutation defined for trait of type $(typeof(x))")
end

mutation(x, mu_m; kwargs...) = _trait_type_error(x, mu_m)

#--- Boolean traits
function mutation(trait_value::Bool, mu_m::Float64; kwargs...)
    rand() < mu_m ? !trait_value : trait_value
end

#--- Discrete traits
function mutation(trait_value::Int, mu_m::Float64; boundaries, kwargs...)
    if rand() < mu_m
        return(random_int_except(boundaries[1],boundaries[2],trait_value))
    else
        return(trait_value)
    end
end

#--- Continuous traits
## Wrapper that dispatches to the type-specific mutation function via `Val{T}`. Each method is specialised by Julia's compiler.
@inline function mutation(x::Float64, mu_m::Float64; mutation_type=:normal, kwargs...)
    return mutation(x, mu_m, Val(mutation_type); kwargs...)
end

# Mutates a Float64 trait using the specified mutation type `T` (e.g. :normal, :gumbel).
# Specialising on `Val{T}` lets Julia generate fast, type-specific code for each mutation type.
function mutation(trait_value::Float64, mu_m::Float64, ::Val{T}; kwargs...) where T
    if rand() < mu_m
        d = _get_mutation_distribution(trait_value, Val(T); kwargs...)
        return rand(d)
    else
        return trait_value
    end
end


#--- Multiple trait 
#@ Always mutating is faster: branching creates type instability for tuples; Julia specializes well on fixed-size tuples.
function mutation(trait_tuple::Tuple, mu_tuple; kwargs_per_trait...)
    ntuple(i -> mutation(trait_tuple[i], mu_tuple[i]; kwargs_per_trait[i]...), length(mu_tuple))
end

#--- Sexual
function mutation(genotype::AbstractMatrix, mu_m; kwargs...)
    return mutation.(genotype, mu_m;kwargs...)
end


"""
    _get_mutation_distribution(x::Float64, ::Val{:type}; sigma_m, boundaries, kwargs...)

Returns a mutation distribution centered around the trait value `x`, used to generate bounded mutations.

Supported mutation types:
- `:normal`: Truncated `Normal(x, sigma_m)` distribution.
- `:normal_clamped`: Clamped `Normal(x, sigma_m)` distribution.
- `:gumbel`: Truncated `Gumbel(x + bias_m, sigma_m)` distribution.

Arguments:
- `x::Float64`: Current trait value.
- `::Val{:type}`: Mutation type (`Val(:normal)`, `Val(:gumbel)`, etc.).
- `sigma_m::Float64`: Standard deviation or scale parameter of the distribution.
- `boundaries::Tuple{Float64, Float64}`: Lower and upper bounds.
- `bias_m::Float64` (optional): Shift used in the `:gumbel` distribution.

Returns:
- A distribution that can be sampled with `rand`.

Throws:
- An error if the mutation type is not implemented.

Examples:
```julia
_get_mutation_distribution(0.5, Val(:normal); sigma_m=0.1, boundaries=(0.0, 1.0))
_get_mutation_distribution(0.5, Val(:gumbel); sigma_m=0.1, bias_m=0.2, boundaries=(0.0, 1.0))
"""
function _get_mutation_distribution(x::Float64, ::Val{:normal}; sigma_m, boundaries, kwargs...)
    Truncated(Normal(x, sigma_m), boundaries...)
end

function _get_mutation_distribution(x::Float64, ::Val{:normal_clamped}; sigma_m, boundaries, kwargs...)
    censored(Normal(x, sigma_m), boundaries...)
end

#---Gumbel distribution for biased mutation (see Henrich paper 2017 Tasmanian ...)
function _get_mutation_distribution(x::Float64, ::Val{:gumbel}; sigma_m, bias_m=0.0, boundaries, kwargs...)
    Truncated(Gumbel(x + bias_m, sigma_m), boundaries...)
end

function _get_mutation_distribution(x::Float64, ::Val{T}; kwargs...) where T
    error("Unsupported mutation type: $T. To see implemented mutation type, check doc of _get_mutation_distribution")
end


"""
    mutation!(population::Vector{Float64}, mu_m; mutation_type=:normal, boundaries=(0.0, 1.0), kwargs...)

Efficient in-place mutation of a vector of continuous traits.

This specialised method avoids repeated dispatch and keyword parsing by resolving the mutation type and parameters once before the loop. It mutates only a subset of traits selected with probability `mu`, drawing values according to the selected mutation type.

Arguments:
- `population::Vector{Float64}`: Vector of continuous trait values.
- `mu_m`: Mutation probability.
- `mutation_type::Symbol`: Type of distribution (`:normal`, `:gumbel`, etc.).
- `boundaries::Tuple{Float64, Float64}`: Lower and upper bounds for the mutated trait values.
- `kwargs...`: Additional parameters for the mutation distribution (e.g. `sigma_m`, `bias_m`).

Returns:
- Mutates `population` in place. No return value.

Benchmark tip:
```julia
pop = rand(1000)
@btime mutation!(pop, 0.1; sigma_m=0.1, boundaries=(0.0, 1.0))
pop2 = rand(1000)
@btime for i in eachindex(pop2)
   pop2[i] = mutation(pop2[i], 0.1; sigma_m=0.1, boundaries=(0.0,1.0))
end

Examples:
pop = fill(0.5, 1000)
mutation!(pop, 0.3; sigma_m=0.1, boundaries=(0.0, 1.0))
mutation!(pop, 0.01; mutation_type=:gumbel, sigma_m=0.1, bias_m=0.2, boundaries=(0.0, 1.0))
```
"""
function mutation!(population::Vector{Float64},mu_m; mutation_type::Symbol = :normal, boundaries = (0.0, 1.0), kwargs...)
    T = Val(mutation_type)
    @inbounds for i in eachindex(population)
        if rand() < mu_m
            d = _get_mutation_distribution(population[i], T; boundaries=boundaries, kwargs...)
            population[i] = rand(d)
        end
    end
end

function mutation!(population::Vector{Vector{Float64}},mu_m; mutation_type::Symbol = :normal, boundaries = (0.0, 1.0), kwargs...)
    T = Val(mutation_type)
    @inbounds for i in eachindex(population)
        for j in eachindex(population[i])
            if rand() < mu_m
                distribution = _get_mutation_distribution(population[i][j], T; boundaries=boundaries, kwargs...)
                population[i][j] = rand(distribution)
            end
        end
    end
end


function mutation!(population::Vector{T},mu_m; kwargs...) where T
    #@ Apparently, broadcasting does not allow conditional skipping so better to use loop
    @inbounds for i in eachindex(population)
        population[i] = mutation(population[i], mu_m; kwargs...)
    end
end

function mutation!(population::Vector{Vector{T}},mu_m; kwargs...) where T
    @inbounds for j in eachindex(population)
        mutation!(population[j],mu_m; kwargs...)
    end
end

"""
    mutation!(genotype::Matrix{Float64}, mu_m; mutation_type=:normal, kwargs...)

In-place mutation of a diploid genotype represented as a `Matrix{Float64}`.

Each allele (matrix element) has an independent probability `mu` of mutating. The mutated
value is drawn from a distribution specified by `mutation_type` and `kwargs...`.

# Arguments
- `genotype`: A 2D matrix representing diploid alleles (rows: loci, columns: alleles).
- `mu_m`: Per-allele mutation probability.
- `mutation_type`: Type of mutation distribution (`:normal`, `:gumbel`, etc.).
- `kwargs...`: Additional mutation parameters (e.g., `sigma_m`, `boundaries`, `bias_m`).

# Returns
- The genotype is mutated in place.
"""
#--- Non-haploid with allelic values (multiple loci)
function mutation!(genotype::Matrix{Float64},mu_m; mutation_type=:normal, kwargs...)
    T = Val(mutation_type)
    for i in axes(genotype,1), j in axes(genotype,2)
        if rand() < mu_m
            distribution = _get_mutation_distribution(genotype[i,j], T; kwargs...)
            genotype[i,j] = rand(distribution)
        end
    end
end

#--- Non-haploid with +/- effect (multiple loci)
function mutation!(genotype::BitMatrix,mu_m; mutation_type=:normal, kwargs...)
    T = Val(mutation_type)
    for i in axes(genotype,1), j in axes(genotype,2)
        if rand() < mu_m
            genotype[i,j] = !(genotype[i,j])
        end
    end
end

function mutation!(genotype::Vector{Matrix},mu_m; mutation_type=:normal, kwargs...)
    @inbounds for i in eachindex(genotype)
        mutation!(genotype[i], mu_m;mutation_type, kwargs...)
    end
end

function mutation!(genotype::Vector{Vector{Matrix}},mu_m; mutation_type=:normal, kwargs...)
    @inbounds for i in eachindex(genotype)
        for j in eachindex(genotype[i])
            mutation!(genotype[i][j], mu_m;mutation_type, kwargs...)
        end
    end
end





