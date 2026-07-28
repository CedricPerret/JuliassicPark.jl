module JuliassicPark

using SplitApplyCombine: invert
using IterTools: partition
#To avoid conflict on mean with StatsBase
using Distributions
using StatsBase
using DataFrames
using Random
using CSV
using ProgressMeter

#--- Top-level user interface. Contains evol_model.
include("simulation.jl")
#--- I/O
include("input.jl")
include("output.jl")
#--- Functions used to simulate different part of the life cycle
include("mutation.jl")
include("reproduction.jl")
include("migration.jl")
#--- Genotype to phenotype mapping useful for sexual reproduction
include("mapping.jl")
#--- utility functions
include("utils.jl")


export evol_model,
mutation, mutation!,
list_reproduction_methods, list_reproduction_functions, list_migration_methods, list_migration_functions,
reproduction_WF, reproduction_WF!, reproduction_Moran_DB!, reproduction_Moran_BD!, reproduction_Moran_pairwise_learning!,
reproduction_WF_copy_group_trait!, 
reproduction_WF_island_model_fixed_migrant_fraction,reproduction_WF_island_model_fixed_migrant_fraction!,
reproduction_WF_island_model_hard_selection, creproduction_WF_island_model_hard_selection!, reproduction_WF_island_model_soft_selection, reproduction_WF_island_model_soft_selection!,
reproduction_explicit_poisson,
reproduction_WF_sexual, reproduction_WF_sexual!,
average_mapping, additive_mapping, recombination,
random_migration, random_migration!,
initialise_population,
get_default_parameters, print_parameters, print_default_parameters, set_default_parameters!,  reset_default_parameters!,
random_int_except, remove_index, sample_except, create_interval, create_interval_from_size, fill_array_with_missing,  invert_3D,
curve_plateau,curve_plateau_sublinear, curve_sigmoid, curve_sigmoid_decreasing,
vv_empty,vv,vv_rand,vv_big,
coef_linear_regression, normalised, random_pairing, random_grouping, argextreme, power!,
get_parameters_from_sweep, lift_map,add_eps!, get_ind,
@extras

end