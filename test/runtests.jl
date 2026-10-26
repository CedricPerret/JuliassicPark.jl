using JuliassicPark
using Test

"""
run_smoke_tests(parameters, fitness_function;
                methods = list_reproduction_functions(),
                details = ['i','p','g'],
                structures = (:pop, :metapop),
                io = stdout)

For each reproduction method and data resolution `log_by`, run evol_model once and report:
  [OK]    method log_by structure
  [ERROR] method log_by structure: error message

Skips methods that require extra positional args not wired by evol_model
(e.g. reproduction_WF_copy_group_trait).

This is a smoke test: success = no exception thrown.
"""
function run_smoke_tests(parameters, fitness_function;
                         methods = list_reproduction_functions(),
                         details = ['i','p','g'],
                         structures = (:pop, :metapop),
                         io = stdout)

    errors = String[] 

    for m in methods
        # Skip methods that need extra positional args evol_model does not pass
        if any(need in (:group_level_trait, :group_fitness_fun) for need in m.needs)
            continue
        end
        # Skip sexual when we test non-float
        # if contains(string(m.name),"sexual") &&
        #    ((any(typeof.(parameters[:z_ini]) .== Bool)) || any(typeof.(parameters[:z_ini]) .== Int))
        #     continue
        # end
        # Skip sexual with non-float traits, except multilocus with Boolean alleles (default additive mapping)
        if contains(string(m.name), "sexual") &&
        (any(typeof.(parameters[:z_ini]) .== Int) ||
            (any(typeof.(parameters[:z_ini]) .== Bool) && get(parameters, :n_loci, 0) <= 1))
            continue
        end

        supported = m.applies_to === :both ? (:pop, :metapop) : (m.applies_to,)
        for structure in (s for s in structures if s in supported)
            for log_by in details
                # fresh params, keep the run very small
                params = deepcopy(parameters)
                params[:log_by]      = log_by
                params[:n_gen]   = 5
                params[:patch_size]   = 100
                params[:log_from] = 1
                params[:log_every] = 2

                # minimal extras
                if :mig_rate in m.needs && !haskey(params, :mig_rate)
                    params[:mig_rate] = 0.1
                end
                if :n_loci in m.needs && get(params, :n_loci, 0) == 0
                    params[:n_loci] = 1
                end

                if structure === :pop
                    params[:n_patch]   = 1
                    params[:flatten_single_patch]  = true
                else
                    params[:n_patch]   = 2
                    params[:flatten_single_patch]  = false
                end

                try
                    evol_model(params, fitness_function, m.f)
                    println(io, "[OK]    ", lpad(string(m.name), 45),
                            "  log_by=", log_by, "  structure=", structure)
                    flush(io)
                catch err
                    msg = sprint(showerror, err)
                    println(io, "[ERROR] ", lpad(string(m.name), 45),
                            "  log_by=", log_by, "  structure=", structure, "  : ", msg)
                    flush(io)
                    push!(errors,
                          "[ERROR] $(string(m.name)) log_by=$(log_by) structure=$(structure): $(msg)")
                end
            end
        end
    end

    return errors 
end



#*** Single trait 
function gaussian_fitness_function(z::Number; optimal, sigma, args...)
    # We add 1 to avoid that the population collapses with explicit pop dynamics
    fitness = 1. + exp(-(z - optimal)^2 / sigma^2)
    distance_to_optimal = (z - optimal)^2
    return (;fitness, distance_to_optimal)
end

parameters_example = Dict(
    :z_ini => 0.1,
    :n_gen => 5,
    :patch_size => 100,
    :n_patch => 2,
    :str_selection => 1.0,
    :mu_m => 0.005,
    :sigma_m => Derived((; kwargs...) -> rand([0.1, 0.2])),
    :boundaries => [0.0, 1.0],
    :write_file => false,
    :log_by => 'g',
    :optimal => 0.6,
    :sigma => 1.0,
    :log_every => 1
)



#res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF)
@testset "smoke tests – single trait" begin
    errors = run_smoke_tests(parameters_example, gaussian_fitness_function)
    @test isempty(errors)
end

function parameter_constraints_test(p)
    @assert p[:sigma] > 0 "Standard deviation must be positive"
end

#As first version. To improve later on.
@testset "basic invariants – simple WF run" begin
    local params = deepcopy(parameters_example)
    params[:n_gen] = 5
    params[:patch_size] = 100
    params[:n_patch] = 2
    params[:log_by] = 'i'

    res = evol_model(params, gaussian_fitness_function, reproduction_WF; parameter_constraints = parameter_constraints_test)

    @test !isempty(res)
    @test "fitness" ∈ names(res)
    @test "z" ∈ names(res) 

    # Invalid parameters
    params[:sigma] = 0.0
    @test_throws AssertionError evol_model(params, gaussian_fitness_function, reproduction_WF; parameter_constraints = parameter_constraints_test)
end

#*** Two traits float
function gaussian_fitness_function_two_traits(z::Tuple; optimal, sigma, args...)
    # We add 1 to avoid that the population collapses with explicit pop dynamics
    fitness = 1. + exp(-(z[1] - optimal)^2 / sigma^2)
    distance_to_optimal = (z[1] - optimal)^2
    return (;fitness, distance_to_optimal)
end

parameters_example_two = Dict(
    :z_ini => (0.1,0.1),
    :n_gen => 5,
    :patch_size => 100,
    :n_patch => 2,
    :str_selection => 1.0,
    :mu_m => 0.005,
    :sigma_m => (0.1, 0.1),
    :boundaries => ([0.0, 1.0],[0.0, 1.0]),
    :write_file => false,
    :log_by => 'g',
    :optimal => 0.6,
    :sigma => 1.0,
    :log_every => 1
)

## Debugging
#res = evol_model(parameters_example_two, gaussian_fitness_function_two_traits, reproduction_WF)
@testset "smoke tests – two traits" begin
    errors = run_smoke_tests(parameters_example_two, gaussian_fitness_function_two_traits)
    @test isempty(errors)
end


#*** Three traits of different types
function test_fitness_function(z::Tuple; optimal, sigma, c, args...)
    # We add 1 to avoid that the population collapses with explicit pop dynamics
    fitness = 1. + exp(-(z[1] - optimal)^2 / sigma^2)  - c * z[2] + rand() * float(z[3])
    distance_to_optimal = (z[1] - optimal)^2
    return (;fitness, distance_to_optimal)
end

parameters_example_three = Dict(
    :z_ini => (0.1, true, 5),
    :n_gen => 3,
    :patch_size => 100,
    :n_patch => 1,
    :n_loci => 0,
    :str_selection => 1.0,
    :mu_m => [0.005,0.005,0.005],
    :sigma_m => (0.1, nothing, nothing),
    :boundaries => ([0.0, 1.0], nothing, [0, 5]),
    :write_file => false,
    :log_by => 'i',
    :optimal => 0.6,
    :sigma => 1.0,
    :c => 0.5
)

## Debugging
@testset "smoke tests – three mixed traits" begin
    errors = run_smoke_tests(parameters_example_three, test_fitness_function)
    @test isempty(errors)
end

#*** Sexual, multiple loci (default additive mapping, Boolean alleles)
function fitness_multilocus(z::Number; optimal, sigma, args...)
    return 1. + exp(-(z - optimal)^2 / sigma^2)
end

parameters_multilocus = Dict(
    :z_ini => [true, false],       # each allele drawn at random
    :n_loci => 5,
    :allele_eff => 0.1,
    :n_gen => 5,
    :patch_size => 100,
    :n_patch => 2,
    :str_selection => 1.0,
    :mu_m => 0.01,
    :write_file => false,
    :optimal => 0.2,
    :sigma => 1.0
)

@testset "smoke tests – sexual, multiple loci" begin
    sexual_methods = filter(m -> contains(string(m.name), "sexual"), list_reproduction_functions())
    errors = run_smoke_tests(parameters_multilocus, fitness_multilocus; methods = sexual_methods)
    @test isempty(errors)
end