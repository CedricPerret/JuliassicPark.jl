#***********************************************
#*** Migrations function
#***********************************************
#! In-place migration function returns false to skip migration.

"""
    random_migration(metapop; mig_rate)
    random_migration!(metapop, new_metapop; mig_rate)

Redistribute individuals across patches by independent migration.

Each individual has probability `mig_rate` of moving to a randomly chosen *other* patch.

# Arguments
- `metapop::Vector{Vector{T}}`: A list of patches, each containing individuals.
- `new_metapop, metapop`: Preallocated destination and source metapopulations for the in-place version.
- `mig_rate::Real`: Probability of migration for each individual.
- `kwargs...`: Ignored; included for interface compatibility.

# Behaviour
- `random_migration` is the allocating version: it constructs and returns a new metapopulation.
- `random_migration!` is the in-place version and the recommended entry point: **it automatically chooses the fastest method**.
    - If the expected number of migrants is low, it uses the fully in-place kernel (`_random_migration!`), avoiding allocations.
    - If many migrants are expected, it calls the allocating version internally and copies into `new_metapop`.

# Returns
- `random_migration`: a new metapopulation.
- `random_migration!`: `nothing`; the result is written into `new_metapop`.

# Notes
- Migration never returns an individual to its own patch (`random_int_except` is used).
- If `mig_rate == 0` or only one patch exists, return flag to not do migration
"""
function random_migration!(metapop::Vector{Vector{T}}, new_metapop ; mig_rate::Float64, kwargs...) where T
    n_patches = length(metapop)
    #-> If no migration possible, skip
    (n_patches == 1 || mig_rate == 0.) && return false
    
    ## Cutoff chosen by trial and error.
    if sum(length.(metapop)) * mig_rate < 2000
        #-> Few migrants: use the non-allocating version
        _random_migration!(metapop, new_metapop ; mig_rate)
    else
        #-> Many migrants: use the allocating version and copy results in-place
        copyto!(new_metapop, random_migration(metapop; mig_rate))
    end
    return nothing
end

function random_migration(metapop::Vector{Vector{T}}; mig_rate::Float64, kwargs...) where T
    n_patches = length(metapop)
    if n_patches == 1 || mig_rate == 0.
        return(copy(metapop))
    end
    #--- Identify the migrants
    migrants_index = [rand(length(group)) .< mig_rate  for group in metapop]
    #--- Get the traits of the migrants
    migrants = [group[idxs] for (group, idxs) in zip(metapop, migrants_index)]
    #--- Remove the migrants to start generating the new metapopulation
    new_metapop = [group[.!idxs] for (group, idxs) in zip(metapop, migrants_index)]
    #--- Make the migrants migrate to new random patch
    for i in eachindex(metapop)
        #-> For each patch
        for j in eachindex(migrants[i])
            #-> and each migrant within this patch
            ## Make them go to another random patch
            push!(new_metapop[random_int_except(1,n_patches,i)],migrants[i][j])
        end
    end
    return(new_metapop)
end

function _random_migration!(metapop::Vector{Vector{T}}, new_metapop ; mig_rate::Float64, kwargs...) where T
    n_patches = length(metapop)
    #--- Start by copying old population into new_metapop
    for group_i in 1:n_patches
        empty!(new_metapop[group_i])
        append!(new_metapop[group_i], metapop[group_i])
    end
    #--- Remove migrants in-place
    for group_i in 1:n_patches
        # We iterate from the end so that deleteat! does not shift upcoming indices.
        # deleteat!(A, i) only shifts elements *after* index i. Any elements *before* i keep their indices.
        for ind_i in length(metapop[group_i]):-1:1
            if rand() < mig_rate
                push!(new_metapop[random_int_except(1,n_patches,group_i)], metapop[group_i][ind_i])
                deleteat!(new_metapop[group_i], ind_i)
            end
        end
    end
end


#*** Programmatic registry of migration functions
# Single source of truth for migration methods
_MIGRATION_FUNS = [
    (name = :random_migration!,
     f = random_migration!,
     desc = "Independent migration (auto-select in-place or allocating version)",
     needs = [:mig_rate]),

    (name = :random_migration,
     f = random_migration,
     desc = "Independent migration, allocating version",
     needs = [:mig_rate]),
]


list_migration_functions() = _MIGRATION_FUNS

# Pretty-printer driven by the registry
function list_migration_methods(io::IO=stdout)
    println(io, "Available migration methods\n")
    if isempty(_MIGRATION_FUNS)
        println(io, "  (none registered)")
        return
    end
    funs = sort(_MIGRATION_FUNS; by = x -> String(x.name))
    w = maximum(length(string(x.name)) for x in funs)
    for x in funs
        n = string(x.name)
        pad = repeat(" ", w - length(n) + 2)
        println(io, "  ", n, pad, "- ", x.desc)
    end
end


