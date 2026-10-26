#**********************************************************
#***             JuliassicPark.jl Example Usage
#**********************************************************

## Load JuliassicPark
using JuliassicPark

## Optional packages for examples:
# - Plots.jl and DataFramesMeta.jl: to visualise results
# - Distributions.jl: to initialise traits from probability distributions
# - BenchmarkTools.jl: to test performance

using Plots
using DataFramesMeta
using Distributions
using BenchmarkTools

# If you don’t have them already, install with:
# using Pkg
# Pkg.add(["Plots", "DataFramesMeta", "Distributions", "BenchmarkTools"])

#-----------------------------------------------------------
#*** 1. Minimal Example — Random Fitness
#    Demonstrates the minimum required structure to run evolutionary model
#-----------------------------------------------------------

function dummy_fitness_function(z::Number; kwargs...)
    return rand() + z
end

parameters_example = (
    z_ini = Normal(0.1, 0.1),
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0]
)

res = evol_model(parameters_example, dummy_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 2. Gaussian Fitness Function
#    Demonstrates setting parameters and including an additional output to be saved
#-----------------------------------------------------------

function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return fitness, distance_to_optimal
end

parameters_example = Dict(
    :z_ini => Normal(0.1, 0.1),
    :n_gen => 1500,
    :patch_size => 10000,
    :n_patch => 1,
    :str_selection => 1.0,
    :mu_m => 0.005,
    :sigma_m => 0.1,
    :boundaries => [0.0, 1.0],
    :extras_names => ["distance_to_optimal"],
    :write_file => false,
    :parameters_to_omit => ["n_loci", "log_every"],
    :log_by => 'g',
    :optimal => 0.6,
    :sigma => 1.0
)


res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
@with res plot(:gen, :global_mean_distance_to_optimal)

# You can modify parameters directly to run different simulations
parameters_example[:n_gen] = 300
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
@with res plot(:gen, :global_mean_z, ylims = [0, 1])

## Use Named tuple output to avoid specifying the name of other output.
function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return (;fitness, distance_to_optimal)
end

parameters_example[:extras_names] = []
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 2.1 Conditional output with @extras
#    Only compute extra output when necessary for performance
#-----------------------------------------------------------

function gaussian_fitness_function(z::Number; optimal, sigma, should_it_print = true, kwargs...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    @extras begin
        distance_to_optimal = (z - optimal)^2
        sleep(0.001)  # Simulate a costly computation
    end
    return fitness, distance_to_optimal
end

parameters_example[:n_gen] = 200 
parameters_example[:patch_size] = 5 
parameters_example[:log_every] = 1
@time evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

parameters_example[:log_every] = 100
@time evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 2.2 Output Configurations
#-----------------------------------------------------------

#--- (a) No output file
parameters_example[:write_file] = false
parameters_example[:split_rep] = false
parameters_example[:split_sweep] = false
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#--- (b) Write to single file
parameters_example[:write_file] = true
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#--- (c) Write to single file + Thread-parallelised simulations
parameters_example[:split_rep] = false
parameters_example[:split_sweep] = false
parameters_example[:n_rep] = 15
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#--- (d) Split by simulation
parameters_example[:split_rep] = true
parameters_example[:n_rep] = 2
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)

#-----------------------------------------------------------
#*** 2.3 Parameter Sweep
#-----------------------------------------------------------

parameter_sweep = Dict(:sigma => [1.0, 2.0], :mu_m => [0.05, 0.1])
parameters_example[:n_rep] = 2
parameters_example[:log_every] = 1999
parameters_example[:log_by] = 'g'

#--- (a) All results in one file
parameters_example[:split_sweep] = false
parameters_example[:split_rep] = false
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)

#--- (b) One file per parameter set
parameters_example[:split_sweep] = true
parameters_example[:split_rep] = false
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)

#--- (c) One file per replicate and per parameter set
parameters_example[:split_sweep] = true
parameters_example[:split_rep] = true
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF; sweep = parameter_sweep)


#-----------------------------------------------
#*** 3. Threshold public good game
# Demonstrates (i) boolean trait, (ii) fitness function applying to group
#-----------------------------------------------

#--- Model description
# - Cooperators pay a cost `c`
# - Groups receive a public good if enough individuals cooperate.

function threshold_public_good_game(z::Vector{Bool};a,c,threshold,kwargs...)
    n_player = length(z)
    payoff = ((a * n_player) * (sum(z) > threshold*n_player)) .- (fill(c,length(z)) .* z)
    fitness = c .+ payoff
    return fitness
end

parameters_example = (z_ini = [true,false], mu_m = 0.001,
n_gen = 100, patch_size = 50,n_patch=100, log_by = 'g', n_rep = 100, 
a=2,c=1., threshold = 0.65)

res=evol_model(parameters_example,threshold_public_good_game,reproduction_WF)

#--- Results: Two equilibria. Some simulations converge toward enough cooperators, others to full defection
@with res plot(:gen,:global_mean_z,group=:i_simul,ylims=[0,1],legends=false)

#-----------------------------------------------
#*** 4. Dyadic game
# Demonstrates (i) structured population, (ii) integer trait, (iii) derived parameters
#-----------------------------------------------

#--- Model description
# - Individuals are paired randomly and play a dyadic game (e.g. Prisoner's Dilemma)

using Random

# Construct the payoff matrix from R, S, T, P
function define_payoff_matrix(; R, S, T, P, kwargs...)
    [P T; S R]
end


function get_payoff_from_payoff_matrix(pair_of_players, _payoff_matrix)
    i, j = Int(pair_of_players[1]) + 1, Int(pair_of_players[2]) + 1
    return [_payoff_matrix[i, j], _payoff_matrix[j, i]]
end

# Fitness function: pair individuals, compute fitness using payoff matrix
function fitness_function_dyadic_game(population::Vector; _payoff_matrix, kwargs...)
    index_population_by_pairs = random_pairing(collect(1:length(population)))
    fitness = zeros(length(population))
    for i in 1:length(index_population_by_pairs)
        fitness[index_population_by_pairs[i]] .= get_payoff_from_payoff_matrix(
            population[index_population_by_pairs[i]], _payoff_matrix
        )
    end
    return fitness
end

#--- Simulation: unstructured population (single patch)
# Integer traits (0 = defector, 1 = cooperator)

reset_default_parameters!()

parameters_example = (
    z_ini = 1,
    mu_m = 0.005,
    boundaries = [0, 1],
    patch_size = 1000,
    log_by = 'g',
    n_gen = 2000,
    R = 4.75,
    S = 0,
    T = 5,
    P = 0,
    _payoff_matrix = Derived(define_payoff_matrix)
)

#-> Results: Cooperation collapses in well-mixed population
res = evol_model(parameters_example, fitness_function_dyadic_game, reproduction_WF)

@with res plot(:gen, :global_mean_z, group = :i_simul, ylims = [0, 1], legend = false)

#-> Partial cooperation is maintained in presence of spatial structure
parameters_example = (z_ini = 1, mu_m = 0.005, boundaries= [0,1],patch_size=10, log_by ='g',n_gen=2000, 
n_patch = 500, mig_rate= 0.05,
R = 4.75, S = 0, T = 5, P = 0,_payoff_matrix = Derived(define_payoff_matrix))

res=evol_model(parameters_example,fitness_function_dyadic_game,reproduction_WF_island_model_hard_selection!)

@with res plot(:gen,:global_mean_z,group=:i_simul,ylims=[0,1],legends=false)

#-----------------------------------------------
#*** 5. Conflict game
# Demonstrates (i) fitness function taking metapopulation as input
#-----------------------------------------------

#--- Model description
# - Individuals invest in costly traits that benefit their group
# - Group success depends on the aggregate investment (nonlinearly)

function conflict_game(z::Vector{Vector{Float64}}; k, benefit, cost, kwargs...)
    output_conflict = (sum.(z) .^ k) ./ sum(sum.(z) .^ k)
    fitness = [
        ones(length(z[i])) .+ fill(benefit * output_conflict[i], length(z[i])) .- (z[i] .* cost)
        for i in 1:length(z)
    ]
    return fitness, output_conflict
end

reset_default_parameters!()

parameters_example = (
    z_ini = 0.2,
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0],
    n_patch = 10,
    patch_size = 50,
    n_gen = 5000,
    log_by = 'p',
    extras_names = ["Benefit_group"],
    k = 5,
    cost = 0.1,
    benefit = 5.0
)

res = evol_model(parameters_example, conflict_game, reproduction_WF)

#--- Results: Group-level conflict dynamics can lead to cycles of increase and decrease in hostility
@with res plot(:gen, :mean_z, group = :patch, ylims = [0, 1], legend = false)


#-----------------------------------------------
#*** 6. Two traits (Cobb–Douglas fitness landscape)
# Demonstrates (i) multiple traits
#-----------------------------------------------

#--- Model description
# - Each trait is under stabilizing selection toward a separate optimum

function cobb_douglas(z::Tuple; o_x, o_y, sigma, alpha, args...)
    fitness = (exp(-((z[1] - o_x)^2) / sigma^2) ^ alpha) *
              (exp(-((z[2] - o_y)^2) / sigma^2) ^ (1 - alpha))
    return fitness
end

parameters_example = (
    z_ini = (0.0, 0.0),                       
    mu_m = [0.005, 0.005],                  
    sigma_m = [0.1, 0.1],                    
    boundaries = [[-5.0, 5.0], [-5.0, 5.0]], 
    n_gen = 5000,
    n_rep = 20,
    log_by = 'g',
    o_x = -2.0,
    o_y = 2.0,
    sigma = 3.0,
    alpha = 0.5
)

res = evol_model(parameters_example, cobb_douglas, reproduction_WF)

#--- Results: Traits diverge toward separate optima
@with res plot(:gen, :global_mean_z1, group = :i_simul, ylims = [-2.5, 2.5], legend = false)
@with res plot!(:gen, :global_mean_z2, group = :i_simul, ylims = [-2.5, 2.5], legend = false)

#-----------------------------------------------
#*** 7. Ecological model of growth
# Demonstrates (i) parameters computed at the beginning of the simulation, (ii) varying population size, (iii) migration
#-----------------------------------------------

#--- Model description
# - Each group has its own carrying capacity `K` drawn from a Normal distribution

function ecological_fitness(z::Vector{<:Vector}; r, K, kwargs...)
    fitness = [fill(r / (1 + length(z[i]) / K[i]), length(z[i])) for i in 1:length(z)]
    group_size = length.(z)
    return (;fitness, group_size)
end

function draw_K(; mean_K, sd_K, n_patch, kwargs...)
    rand(Normal(mean_K, sd_K), n_patch)
end

reset_default_parameters!()

parameters_example = Dict(
    :z_ini => true,
    :patch_size => 10, :n_patch => 20,
    :n_gen => 100,
    :mean_K => 100, :sd_K => 30,
    :log_by => 'p',
    :r => 2,
    :mig_rate => 0.05,
    :K => Derived(draw_K)
)


#--- Results: Groups grow to different sizes based on drawn K values
res = evol_model(parameters_example, ecological_fitness, reproduction_explicit_poisson, migration_function= random_migration!)
@with res plot(:gen, :group_size, group = :patch, legend = false)

#--- Results: High migration rate homogeneises the group sizes
parameters_example[:mig_rate] = 0.8
res = evol_model(parameters_example, ecological_fitness, reproduction_explicit_poisson, migration_function= random_migration!)
@with res plot(:gen, :group_size, group = :patch, legend = false)

#-----------------------------------------------
#*** 8. Disruptive selection and polymorphism (sexual reproduction)
# Demonstrates (i) diploid traits
#-----------------------------------------------

#--- Model description
# - Compare asexual reproduction with diploid sexual reproduction
# - Fitness favors extreme phenotypes (disruptive selection)

function disruptive_fitness(z::Vector;  strength, kwargs...)
    return exp.( ((z .- mean(z)).^2) ./ strength)
end

reset_default_parameters!()

parameters_example = Dict(
    :z_ini => 0.0,
    :sigma_m => 0.01,
    :mu_m => 0.005,
    :boundaries=>[-1.,1.],
    :n_gen => 2000,
    :patch_size => 500,
    :n_loci => 0,
    :allele_eff => 0.5,
    :log_by => 'i',
    :strength =>0.005,
    :log_every=>20
)

#-> Results: Disruptive selection maintains polymorphism
res = evol_model(parameters_example, disruptive_fitness, reproduction_WF)
@with res scatter(:gen, :z, legend = false)


#-> Sexual reproduction introduces intermediary phenotypes (third branch in the center.)
parameters_example[:n_loci] = 1;
res_sexual = evol_model(parameters_example, disruptive_fitness, reproduction_WF_sexual)
@with res_sexual scatter(:gen, :z,  legend = false)


#-----------------------------------------------
#*** 9. Performance tips: in-place fitness function
# Demonstrates how to use in-place fitness function
#-----------------------------------------------

function gaussian_fitness_function!(z::Vector{Float64}, fitness; optimal, sigma, args...)
    for i in 1:length(z)
        fitness[i] = exp(-(z[i] - optimal)^2 / sigma^2)
    end
    distance_to_optimal = (z .- optimal).^2
    return distance_to_optimal
end

parameters_example = Dict(
    :z_ini => Normal(0.1, 0.1),
    :n_gen => 3000,
    :patch_size => 10000,
    :n_patch => 1,
    :str_selection => 1.0,
    :mu_m => 0.005,
    :sigma_m => 0.1,
    :boundaries => [0.0, 1.0],
    :extras_names => ["distance_to_optimal"],
    :write_file => false,
    :parameters_to_omit => ["n_loci", "log_every"],
    :log_by => 'g',
    :optimal => 0.5,
    :sigma => 1.0
)

res = evol_model(parameters_example, gaussian_fitness_function!, reproduction_WF);

function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return fitness, distance_to_optimal
end

function gaussian_fitness_function_pop(z::Vector{Float64}; optimal, sigma, args...)
    fitness = exp.(.- (z .- optimal).^2 / sigma^2)
    distance_to_optimal = (z .- optimal).^2
    return fitness, distance_to_optimal
end


#@ Performance comparison. 
@btime res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF);
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
#-> faster
@btime res = evol_model(parameters_example, gaussian_fitness_function_pop, reproduction_WF);
@with res plot(:gen, :global_mean_z, ylims = [0, 1])
#-> less allocation
@btime res = evol_model(parameters_example, gaussian_fitness_function!, reproduction_WF);
@with res plot(:gen, :global_mean_z, ylims = [0, 1])

