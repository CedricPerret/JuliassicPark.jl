# JuliassicPark.jl 
![CI](https://github.com/CedricPerret/JuliassicPark/actions/workflows/CI.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21647950.svg)](https://doi.org/10.5281/zenodo.21647950)

**JuliassicPark.jl** is a lightweight and flexible Julia package for simulating evolutionary models with customizable fitness functions.

It is built for researchers and modelers who need both flexibility and convenience. JuliassicPark lets you define complex model logic and rich fitness functions, while handling the rest automatically: simulation loops, parameter management, data output, and parallel execution. 

The goal is simple: spend less time on boilerplate, and more time exploring ideas.

---

## Installation

Install it from General registry:

```julia
using Pkg
Pkg.add("JuliassicPark")
```

---
## Quick start

The main function is `evol_model`, which runs a complete evolutionary simulation. The minimum required is to:

1. Write a **fitness function** that maps trait values to fitness.
2. Define **the evolving traits** by setting keys in a dictionary or named tuple:
   - `z_ini`: initial value
   - `mu_m`: mutation rate
   - `sigma_m`: size of mutation steps (for continuous traits)
   - `boundaries`: trait range (for discrete and continuous traits)
3. Provide values for the **parameters used by your fitness function** by adding them to the same dictionary or named tuple.
4. Choose a **reproduction method** that builds the next generation.

That is enough to run. For example,

```julia
using JuliassicPark

function gaussian_fitness_function(z::Number; optimal, sigma, kwargs...)
    fitness = exp(-(z - optimal)^2 / sigma^2)
    return fitness
end

parameters_example = (
    z_ini = 0.1,
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0],
    optimal = 0.5,
    sigma = 0.3)

res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF!)
```

That is the whole model. You get back a tidy `DataFrame`, ready to plot or analyse:

| i_simul | gen | patch | ind |   z   | fitness |
|---------|-----|-------|-----|-------|---------|
| 42      | 1   | 1     | 1   | 0.100 | 0.169   |
| 42      | 1   | 1     | 2   | 0.100 | 0.169   |
| 42      | 2   | 1     | 1   | 0.104 | 0.175   |
| 42      | 2   | 1     | 2   | 0.100 | 0.169   |
| ⋮       | ⋮   | ⋮     | ⋮   | ⋮     | ⋮       |

**From there, a lot is a couple of parameters away.** Add it to the parameters and rerun.

- Want 10 replicates? `n_rep = 10`.

- A population structured in 20 patches of 10? `patch_size = 10`, `n_patch = 20`.

- Mean trait per patch instead of one row per individual? `log_by = 'p'`. Per generation? `log_by = 'g'`.

-  Five thousand generations, but saved every hundredth? `n_gen = 5000`, `log_every = 100`.

- Results on disk rather than in memory? `write_file = true`.

- Varying a parameter over a whole range? 

```julia
res = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF!;
                 sweep = Dict(:sigma => [0.1, 0.2, 0.4, 0.8]))
```

In every case the output adapts: a column for the replicate or the patch, columns averaged and renamed, a column for whatever you swept.

---

## Examples

Examples are provided in the file `examples/basic_example.jl` in the repository. I recommend reading them together with this README.

Note that:
- The example file is a plain `.jl` script but it is meant to be run step by step like a notebook (as in R).
- These examples use additional packages such as `Plots.jl` and `DataFramesMeta.jl`. They are **not required** for running JuliassicPark itself, but install them if you want to reproduce figures.

---

## Status and feedback

JuliassicPark is currently in a beta `0.x` stage. The main features are implemented and it currently handles continuous and discrete traits, single or multiple traits, phenotypes or explicit multilocus genotypes, structured populations, and a range of reproduction schemes (Wright–Fisher, Moran, explicit agent-based, sexual).

I use it for my own research, but it has not been heavily tested across a wide range of models. If you run into errors, unexpected behaviour, or have suggestions, please contact me at cedric.perret.research [at] gmail [dot] com.

---

## Tutorial

To build an evolutionary model, you define 3 things: (i) the trait(s) that evolve, (ii) what happens at each generation which determines the success of individuals (fitness) and (iii) how the next generation is generated from the current one. Here, you do that by writing a fitness function, by creating a `Dict` or `NamedTuple` with the relevant keys and values explained below (called `parameters` from here on), and by choosing a reproduction function. 

### Defining traits

Define traits by adding the key `z_ini` and its value to `parameters` and, when needed, `boundaries`.

The type of `z_ini` sets the type of the trait:
- **Boolean** (`true`, `false`) for a discrete trait with two possible values, e.g. `true = cooperator`, `false = defector`.
- **Integer** (`0`, `1`, `2`, …) for a discrete trait with multiple possible values, e.g. `0 = defector`, `1 = cooperator`, `2 = tit-for-tat`. 
- **Float** (`0.2`, `1.5`, …) for a continuous trait, e.g. `0.3 = contribution to a public good`.

The value(s) of `z_ini` set the initial value(s) of the trait:
- **Single value** — all individuals receive the same initial trait value.  
- **Vector** — trait values are drawn randomly (with replacement) from the vector.  
- **Distribution** (from *Distributions.jl*) — trait values are sampled from the distribution; if `boundaries` is given, sampling is truncated to that interval.
- **DataFrame** — a DataFrame with columns `gen`, `patch`, and `z` (or `z1`, `z2`, … for multiple traits). See [Model a shock](#model-a-shock).

The values of `boundaries` determine the possible trait values (e.g. `(1, 5)` → {1,2,3,4,5} for Integer, and any number between 1 and 5 for Float). This is not required for Boolean traits. 

For [Multiple Traits](#multiple-traits), use a tuple. Each element defines one trait. For example, `z_ini = (0.5, 2)` means that each individual carries one continuous and one discrete trait. 

#### Mutation

A key element of an evolutionary model is that traits can mutate. At minimum, provide the mutation probability during reproduction by adding the key `mu_m` and its value to `parameters`. The effect of mutation depends on the trait type:

- **Discrete traits with two values (Boolean)** — the value is flipped (`true` ↔ `false`).  
- **Discrete traits with multiple values (Integer)** — the value is replaced by another integer within the allowed range.  
- **Continuous traits (Float)** — a new value is drawn according to `mutation_type` within `boundaries`. Available options are:
  - `:normal` (default) — truncated Normal distribution centred on the current trait value.
  - `:normal_clamped` — Normal mutation with values outside `boundaries` clamped to the nearest boundary.
  - `:gumbel` — truncated Gumbel distribution; use `bias_m` to introduce directional bias.

### What happens at each generation

Define what happens at each generation by writing a fitness function that takes values of traits and gives back the resulting fitness. I made the code flexible so that you can write it as you like, either at the level of:
- **One individual**: takes one individual's traits and returns that individual's fitness. Do this when *fitness depends only on the focal individual*.
- **One group**: takes a vector of individuals' traits and returns one fitness value per individual in that group. Do this when *fitness depends on interactions within a group*.
- **The whole metapopulation**: takes a vector of groups of traits and returns one fitness value for each individual in each group. Do this when *fitness depends on interactions between groups*.

The fitness function should follow this structure:
```julia
function my_fitness_function(trait; param1, param2, kwargs...)
    fitness = ...  # compute fitness
    return fitness
end
```
- It **must have a single non-keyword argument** (that is, a single argument before `;`). This is the trait value(s): a scalar or tuple for a single individual, a vector for a group, or a vector of vectors for an entire metapopulation.  
- It **must have `kwargs...` at the end of its parameter list**.  
- It **must return fitness as the first output**. The fitness must match the structure of the input trait (scalar, vector, or vector of vectors).  

To avoid ambiguity, **specify the type** for the first argument trait (for example `Number`, `Tuple`, `Vector`, or `Vector{<:Vector}`). Otherwise, the package will try to infer it automatically, and if there is an error inside the fitness function, the error may show up in confusing places. 

To run simulations, add any arguments used by your fitness function (e.g. `optimal`, `sigma`) and their values to `parameters`.

### Generating the next generation

Define how the next generation is generated by choosing a reproduction function. The standard methods are already implemented and optimised. You can see what is available with:

```julia
list_reproduction_methods() ## To see the list and requirements
list_reproduction_functions() ## To access directly the functions
```

Some of the most commonly used are:
- `reproduction_WF!` — Wright–Fisher reproduction. Non-overlapping generations.  
- `reproduction_Moran_DB!` — Moran process, death–birth update.  Overlapping generations.
- `reproduction_explicit_poisson` — explicit offspring number, drawn from a Poisson distribution.  
- `reproduction_WF_sexual` — Wright–Fisher reproduction with sexual recombination (diploid, multilocus).  

Note that we use the term *reproduction* in a broad sense. It can also represent processes such as learning or cultural transmission. For instance ``reproduction_Moran_pairwise_learning!`` is the function classically used in models with pairwise learning e.g. Traulsen et al, (2006).


### Running simulation

You are ready to run the simulations. To change settings, like make the simulations longer or the population smaller, simply assign a new value to the relevant parameters in `parameters` (e.g. add `n_gen = 5000` to it). Their names and purposes are listed in the complete [list of parameters](#list-of-parameters), and can be inspected with:
```julia
print_default_parameters()         # Print default values
get_default_parameters()           # Access current defaults
```
To change the defaults for future simulations (see [Workflow](#Workflow)), use:
```julia
set_default_parameters!(...)       # Override defaults globally
reset_default_parameters!()        # Reset to built-in defaults
```

---

### Choosing what to record

Now that your model is up and running, the next step is choosing what you want to record from it. `evol_model` returns a `DataFrame`. Each row gives the trait(s) value(s) and the fitness by generation (`log_by = 'g'`), by patch (`log_by = 'p'`), or by individual (`log_by = 'i'`), depending on `log_by`. They come along with the corresponding generation, patch and individual IDs. Each run also has an ID, which serves as the random seed for reproducibility.

You can save extra variables by adding them to the return of the fitness function
```julia
return fitness, extra1, extra2
```
or as a named tuple:
```julia
return (; fitness, extra1, extra2)
```
Extra variables need to match one of the resolutions: one value per individual, per patch or per generation (see [compute extras only when needed](#compute-extras-only-when-needed) to make it quick).

An example output with a single trait, one extra variable `distance_to_optimal`, a population structured in two patches of size 2, two generations and individual-level resolution (`log_by = 'i'`) is:

| i_simul | gen | patch | ind |   z   | fitness | distance_to_optimal |
|---------|-----|-------|-----|-------|---------|---------------------|
| 42      | 1   | 1     | 1   | 0.25  | 0.0019  | 0.0625              |
| 42      | 1   | 1     | 2   | 0.60  | 0.3679  | 0.0100              |
| 42      | 1   | 2     | 3   | 0.23  | 0.0007  | 0.0729              |
| 42      | 1   | 2     | 4   | 0.57  | 0.6130  | 0.0049              |
| 42      | 2   | 1     | 1   | 0.21  | 0.0002  | 0.0841              |
| 42      | 2   | 1     | 2   | 0.61  | 0.2980  | 0.0121              |
| 42      | 2   | 2     | 3   | 0.23  | 0.0007  | 0.0729              |
| 42      | 2   | 2     | 4   | 0.57  | 0.6130  | 0.0049              |

Column names are taken, in order of priority:
1. from `extras_names`, if you set it;
2. otherwise from the field names, if you returned a named tuple;
3. otherwise get assigned automatically `V1`, `V2`, …

#### Resolution handling

You can choose whether results are recorded at the level of the **generation**, **patch**, or **individual** using the `log_by` parameter. If a variable is defined at a finer level than the requested output, it is averaged and renamed: `mean_<var>` for individual → patch or patch → generation, `global_mean_<var>` for individual → generation. If a variable is defined at a coarser level than the requested output, its value is repeated for the corresponding patches or individuals.

Starting from the individual-level output shown above, changing the resolution to patch-level results (`log_by = 'p'`):

| i_simul | gen | patch | mean_z | mean_fitness | mean_distance_to_optimal |
|---------|-----|-------|--------|--------------|--------------------------|
| 42      | 1   | 1     | 0.43   | 0.185        | 0.036                    |
| 42      | 1   | 2     | 0.40   | 0.307        | 0.039                    |
| 42      | 2   | 1     | 0.41   | 0.149        | 0.048                    |
| 42      | 2   | 2     | 0.40   | 0.307        | 0.039                    |

And for generation-level resolution (`log_by = 'g'`):

| i_simul | gen | global_mean_z | global_mean_fitness | global_mean_distance_to_optimal |
|---------|-----|---------------|---------------------|---------------------------------|
| 42      | 1   | 0.41          | 0.246               | 0.038                           |
| 42      | 2   | 0.41          | 0.228               | 0.044                           |

Here, `global_mean_z` is the average of **all individuals in the population**.

#### Writing on disk

If `write_file = true`, results are saved to disk in a CSV file named `name_model-param1=value1-param2=value2-...csv` (see [Loading the data in R](#Loading_the_data_in_R) to import many of them in R).

Use the field `parameters_to_omit` to exclude parameters from the filename.

To keep filenames short:
- Values > 1000 are abbreviated as `1.0k`
- Float values are rounded to 5 digits
- File parameters are printed using only the filename (e.g. `network.csv` => `network`)
- Distributions are formatted as `Name_param1_param2`, e.g. `Normal_0.0_0.5`

#### Saving multiple runs

By default, multiple runs are saved on the same file. Set `split_rep = true` to save each replicate on a different file, and `split_sweep = true` to save each parameter combination on a different file (see [Parameter sweep](#Parameter_sweep))


---

## Working with more complex models

---

### Multiple Traits

Multiple traits are represented internally as tuples. To include several traits, set `z_ini` to a tuple of initial values, one per trait. For example:

```julia
z_ini = (true, 0.2, 2)
```

#### Mutations

Mutation parameters (such as `mu_m` and `sigma_m`) can be given as a single value or as a tuple.
- If a single value is provided, it is applied to all traits.
- If different values are provided, provide a tuple of the same length as the number of traits.

For parameters that apply only to certain traits (e.g. `sigma_m`), you **must use `nothing`** for traits to which they do not apply. The model cannot infer this automatically. For example:
 ```julia
mu_m = (0.01, 0.01, 0.01)
sigma_m = (nothing, 0.1, nothing)
boundaries = (nothing, [0.0, 1.0], [0, 5])
```

#### Access in fitness function

A population with multiple traits is represented as tuples, e.g. each individual is `(z1, z2, …)`, which means that in your fitness function you access the first trait of individual i in patch j with:
- `z[1]` if `z` is an individual,  
- `z[i][1]` if `z` is a population,  
- `z[j][i][1]` if `z` is a metapopulation.  

For metapopulations, you can use the helper `invert_3D` to reorganize the structure from `[patch][individual][variable]` into `[variable][patch][individual]`.  
This makes it easy to extract one trait across all patches and individuals: the first element of the result corresponds to the first trait across the whole metapopulation.

---

### Migration

Migration between patches is implemented in two ways:

- **Integrated in reproduction**: use a reproduction function that already includes migration, such as `reproduction_WF_island_model_hard_selection!` or `reproduction_WF_island_model_soft_selection`. In this case, you only need to provide `mig_rate` as a parameter.
- **Separate migration step**: specify a migration function, which is applied after reproduction. This function receives the population and arguments from parameters.

See the list of migration functions and their requirements using:

```julia
list_migration_methods()   ## To see the list and requirements
list_migration_functions() ## To access directly the functions
```

---

### Sexual Reproduction & multiple loci

When you set `n_loci > 0`, traits are no longer stored as simple numbers. Instead, each individual carries a genotype, represented as a matrix of alleles (rows = loci, columns = 2 alleles). The fitness function still takes as input the phenotype derived from this genotype, which is generated using the `genotype_to_phenotype_mapping` function.

The following defaults mapping apply:
- For a single locus: the phenotype is the **average** of the two alleles, corresponding to a purely additive model without dominance. Alleles must be continuous (`z_ini` a Float).
- For multiple loci: the phenotype is the **sum of allelic effects**, scaled by `allele_eff`. This implements an additive multilocus model with equal effect size per allele. Alleles must be Boolean (`z_ini = true/false`).

You can override the defaults by providing your own ``genotype_to_phenotype_mapping`` function. Note that this function must be defined at the level of an individual genotype (a single matrix), and not at the population level (a vector of matrices).

Be careful to use a reproduction function containing `sexual` in the name if you have `n_loci > 0`.

---
## Advanced usage
---

### Parameter sweep

Explore multiple parameter values automatically by providing a dictionary where:
- **Keys** are parameter names (symbols)
- **Values** are vectors of values to test

For example, this
```julia
parameter_sweep = Dict(:sigma => [1.0, 2.0], :mu_m => [0.05, 0.1])
evol_model(parameters_example, gaussian_fitness_function, reproduction_WF!; sweep = parameter_sweep)
```
runs 4 simulations
- (`sigma = 1.0`, `mu_m = 0.05`)  
- (`sigma = 1.0`, `mu_m = 0.1`)  
- (`sigma = 2.0`, `mu_m = 0.05`)  
- (`sigma = 2.0`, `mu_m = 0.1`)  

By default, all possible combinations are generated automatically (Cartesian product). If you set `sweep_grid = false`, values are combined **by position** (like `zip` in Julia). Using the same example:
- (`sigma = 1.0`, `mu_m = 0.05`)  
- (`sigma = 2.0`, `mu_m = 0.1`)  

#### Sweep output

By default, all results are returned in a single DataFrame, with each varying parameter included as a separate column. If `split_sweep = true`:
- With `write_file = false`, the function returns a list of DataFrames, one per parameter set.
- With `write_file = true`, each parameter set is saved to a separate file.

 See [Parallelisation and output splitting](#parallelisation-and-output-splitting) for more details on saving and splitting.

---

### Parallelisation and output splitting

JuliassicPark.jl lets you run many simulations side by side. In practice, there are two things to decide:
1) how to split the output files,
2) whether to distribute the work across CPU cores or workers.


#### How to split output files

This is controlled by two parameters:

| `split_sweep` | `split_rep` | Behaviour |
|---|---|---|
| `false` | `false` | All simulations and parameter sets are combined into one DataFrame in memory, or one CSV if `write_file = true`. When memory use is moderate (`log_by ≠ 'i'`), runs are parallelised across threads. |
| `true`  | `false` | One file per parameter set. For each set, all replicates are run and written together. Parallelisation happens across parameter sets only. |
| `true`  | `true`  | One file per replicate and per parameter set. This allows parallelisation across parameter sets and across replicates. This is the most scalable option for long sweeps. |
| `false` | `true`  | Not allowed. All replicates would try to write to the same file. An error is raised. |

When different parameter sets are saved to the same file, the varying parameters are added as columns. When replicates are saved separately, the simulation ID is included in the filename.

#### How to distribute the work

- Threads (shared memory).  
  This is the default when results are kept in memory or when writing a single combined file. Threads are simple to use and fast for medium-sized jobs on one machine.

- Distributed workers (separate processes).  
  If you set `distributed = true`, parameter sets and replicates can be sent to multiple workers. This is useful for large sweeps and cluster jobs. Each worker has its own memory and must be ready to run your model.

#### Preparing distributed runs

If you turn on `distributed`, make sure every worker knows about your code and packages (see the [Julia manual on distributed computing](https://docs.julialang.org/en/v1/manual/distributed-computing/)).

```julia
using Distributed
addprocs(4)  # or what your machine or cluster provides

@everywhere using JuliassicPark

# Put your model code in a file so it can be loaded everywhere.
@everywhere include("my_model_code.jl")
```
This includes any user-defined functions used by the model, including functions used in `Derived`.

---

### Parameters derived at runtime

You can define parameters that are computed **at the start of the simulation**, rather than fixed in advance. This is useful when these values varies across run or depend on other parameters, for example, drawing carrying capacities from a distribution, generating a network based on a chosen network type, or ensuring that initial trait values are always far from the current optimum, whatever that optimum is.

To do this, add it as a normal parameter with as value the function that will be used to generate it, and wrapped in `Derived`. The function can use any argument in the `parameters` dictionary, except other derived parameters. It **must take its arguments only as keywords (after `;`)**, and include `kwargs...` at the end for compatibility.

```julia
function calculate_carrying_capacity(; mean_K, sigma_K, n_patch, kwargs...)
    rand(Normal(mean_K, sigma_K), n_patch)
end

parameters = (
    # ... your trait and simulation parameters
    mean_K = 5,
    sigma_K = 2,
    n_patch = 10,
    K = Derived(calculate_carrying_capacity))

evol_model(parameters, fitness_function, repro_function)
```

#### Saving derived parameters

Since they can be long, derived parameters are saved in the output table instead of the filename. To prevent a derived parameter from being saved, you can either:

- Start its name with an underscore (e.g. `_hidden_variable`), or
- Add its name to the `parameters_to_omit` list.

Any parameter that is saved must have a resolution consistent with the simulation:
- One value per generation,
- One value per patch,
- One value per individual (requires constant group and population size).

#### Sweeping over derived parameters

Sweeps cannot be done over derived parameters. Instead, it is easier to sweep over a parameter feeding into the function, for instance:

```julia
function calculate_carrying_capacity(; K_type, mean_K, sigma_K, n_patch, kwargs...)
    if K_type == :constant
        fill(float(mean_K), n_patch)
    elseif K_type == :normal
        rand(Normal(mean_K, sigma_K), n_patch)
    else
        error("Unknown K_type: $K_type")
    end
end

evol_model(parameters, fitness_function, repro_function; sweep = Dict(:K_type => [:constant, :normal]))
```

---

### Setting parameter constraints

It is often useful to write the constraints on some parameters to avoid that another user (or yourself from the future) run the wrong simulation.

To do this, write a function that takes the parameters as input and checks them with `@assert`, then pass it to `evol_model` through the `parameter_constraints` keyword:
```julia
function check_constraints(p)
    @assert p[:sigma] > 0 "sigma must be positive"
    @assert 0 < p[:eta] < 1 "eta must be between 0 and 1"
end

evol_model(parameters, fitness_function, reproduction_WF!;
           parameter_constraints = check_constraints)
```

The constraints are checked once for each parameter combination, before running the simulations. If one is not satisfied, an error is raised.

---

### Model a shock

To model a change during a simulation, run the first part, then give its result as `z_ini` to a second run with the new parameters:

```julia
res1 = evol_model(parameters_example, gaussian_fitness_function, reproduction_WF!)
res2 = evol_model((; parameters_example..., z_ini = res1, optimal = 0.8), gaussian_fitness_function, reproduction_WF!)
```

The last generation of `res1` is used as the starting population, and `res2` contains both parts, with generations continuing from the last one of `res1`.

Note that:
- The first run must be saved with one row per individual (`log_by = 'i'`, the default).
- The table must contain a single simulation.
- The number of patches and the patch size are taken from the table.
- Derived parameters are computed again. To keep the same values, compute them once and pass them as normal parameters.

---

## Tips

---

### Workflow

#### Setting up a script

A good way to organise an analysis is to set the defaults once, keep a small dictionary for the parameters you want to vary, and modify it as you go.

```julia
set_default_parameters!((
    n_gen = 5000,
    n_patch = 100,
    log_by = 'g',
    write_file = true,
    z_ini = 0.1,
    mu_m = 0.005,
    sigma_m = 0.1,
    boundaries = [0.0, 1.0]))

parameters = Dict(
    :patch_size => 1000,
    :optimal => 0.5,
    :sigma => 0.1)

res = evol_model(parameters, gaussian_fitness_function, reproduction_WF!)

parameters[:patch_size] = 500     # smaller population
res = evol_model(parameters, gaussian_fitness_function, reproduction_WF!)

parameters[:sigma] = 0.3  
res = evol_model(parameters, gaussian_fitness_function, reproduction_WF!)
```

Be careful: changes accumulate. `parameters` keeps everything you modified earlier in the session, so the third run above uses both `patch_size = 500` and `sigma = 0.3`. Rebuild the dictionary when you want a clean starting point.

#### Loading the data in R

If you prefer to analyse your results in R, a small helper script is provided in `scripts/R/import_data.R`. It provides a function to **load and combine all matching `.csv` files in the current working directory into a single `data.table` (or `data.frame`)**. The helper requires the R package `data.table`.

```r
f_import_data(name_file, listVar = c(), parse_value = TRUE, as_data_frame = FALSE)
```

* `name_file` specifies the parts that filenames must contain. Several strings can be provided to successively filter the files.
* `listVar` specifies which parameters in the filenames should be extracted and added as columns.
* `parse_value = TRUE` automatically converts parameter values to numbers or logical values when possible.
*  `as_data_frame = FALSE` controls the output format: `FALSE` (default) returns a `data.table`, while `TRUE` returns a standard `data.frame`.

For example, suppose you ran a parameter sweep with a population size of `1000`, varied `sigma`, and used `split_sweep`, so that each parameter combination was saved in a separate file. You can load and combine all corresponding files with:

```r
data <- f_import_data("patch_size=1000", "sigma")
```

This returns a single `data.table` containing all matching files, with an additional `sigma` column containing the value extracted from each filename.

Several filters and parameters can also be specified:

```r
data <- f_import_data(c("patch_size=1000", "mu_m=0.01"), c("sigma", "n_patch"))
```

Parameters whose values are vectors are automatically expanded into separate columns. For example, a filename containing `optimal=[-2.0,2.0]` will produce the columns `optimal1` and `optimal2`.

Because `-` is used to separate parameters in JuliassicPark filenames, use `_` rather than `-` within multi-word parameter names or string values.

### Performance

#### Compute extras only when needed

If you only save output occasionally (e.g. every 100 generations with `log_every = 100`, see [List of parameters](#list-of-parameters)), you can skip computing extra variables at every step. To do so:
1. Include `should_it_print = true` as a keyword argument in the fitness function.
2. Wrap the code computing the extra variables in an `@extras` block.

```julia
function my_fitness_function(ind; param1, param2, should_it_print = true, kwargs...)
    fitness = ...  # compute fitness
    @extras begin
        extra1 = ...  # only runs if should_it_print == true
        extra2 = ...
    end
    return fitness, extra1, extra2
end
```

**Warning:** 
Avoid reusing variables already defined before `@extras`. If a variable is already defined before the block, it will be overwritten with `NaN` when skipping computation. 


#### In-place fitness functions

JuliassicPark.jl also directly supports **in-place fitness functions** to reduce memory allocation. Instead of returning a newly allocated fitness array, these functions write fitness values into a preallocated array. In-place functions can be defined at either the population or metapopulation level.
An in-place fitness function:
- Takes the current population as the first argument.  
- Takes the fitness array as the second argument.  
- Assigns values directly into `fitness`.  
Example:

```julia
function gaussian_fitness_function!(z::Vector{Float64}, fitness; optimal, sigma, kwargs...)
    for i in 1:length(z)
        fitness[i] = exp(-(z[i] - optimal)^2 / sigma^2)
    end
    distance_to_optimal = (z .- optimal).^2
    return distance_to_optimal
end
```

Here the fitness function writes results into the provided `fitness` vector.  

**Warning:**
    - **Make sure** that the function name ends with `!`.
    - In-place fitness functions only work with reproduction functions that **preserve group sizes**, since the fitness array is preallocated assuming a constant number of individuals per group.

#### Reproduction function

The code runs much faster when:

- you use reproduction functions where group size is constant, since output arrays can be preallocated
- you use in-place reproduction functions (those ending with `!`).

#### Memory allocation

Reducing memory allocation is often one of the most effective ways to speed up a simulation. Fitness is evaluated many times per generation and across many generations, so even small allocations can add up substantially.

To reduce allocations:
- Prefer explicit loops with in-place updates over patterns that create temporary arrays.
- Use the `@.` macro to broadcast all operations in an expression at once. This helps avoid accidental temporary arrays from repeated broadcasting.
- Reuse preallocated buffers when possible rather than creating new arrays inside frequently called functions.
- When slicing arrays, consider `@views` to avoid unnecessary copies.

---
## Design choices
---

### Representation of traits and population

A group is represented as a vector of trait values, and a population as a vector of such groups. The length of the outer vector is the number of patches, and the length of each inner vector is the group size. If there is only one patch, the population can be a single vector of traits.

_Why?_ The alternative would be a matrix representation, which can be faster for fixed group sizes but fails when groups differ in size. Vectors of vectors are more flexible, and Julia’s methods work efficiently with this representation.

Multiple traits for an individual are stored as a tuple.

_Why?_ Tuples are immutable, lightweight, and clearly distinguish “one individual with multiple traits” from “a group of individuals.”

Genotypes are represented as matrices, with rows corresponding to loci and columns to alleles (diploid by default).

_Why?_ This makes genotypes easy to recognize, and phenotypes can be derived quickly from them using mapping functions.

### All as parameters, not agent-based

A common approach is to code each individual as an agent with attributes (traits) and methods (e.g. mutate). This works well in some contexts, but evolutionary models often do not require a full agent-based framework. Most processes reduce naturally to operations on vectors of traits, which are faster and easier to read.

For the same reason, the package does not use an object-oriented style. Julia is not designed for that paradigm, and keeping everything as simple data structures with parameter dictionaries makes the code transparent, lightweight, and closer to the mathematical models used in evolutionary theory.

---

## Full function signature

---

```julia
evol_model(parameters, fitness_function, reproduction_method; 
sweep=Dict{Symbol, Vector}(), migration_function = nothing, genotype_to_phenotype_mapping = identity, parameter_constraints = nothing)
```

### 🧾 Arguments

#### Required:

| Argument              | Type                         | Description                                                                 |
|-----------------------|------------------------------|-----------------------------------------------------------------------------|
| `parameters`          | `Dict` or `NamedTuple`       | All model settings: initial traits, population size, mutation rules, etc. See [Running simulation](#running-simulation). |
| `fitness_function`    | `Function`                   | User-defined function that returns fitness (and optionally extra variables). See [What happens at each generation](#what-happens-at-each-generation).        |
| `reproduction_method` | `Function`       | Function describing how the next generation is built depending on the fitness. Built-in name (e.g. `reproduction_WF!`) or custom function. See [Reproduction](#reproduction). |

#### Optional:
| Argument              | Type                         | Description                                                                 |
|-----------------------|------------------------------|-----------------------------------------------------------------------------|
| `sweep`          | `Dict`       |  A dictionary specifying which parameters to vary across runs. Triggers automatic parameter sweep. See [Parameter Sweep](#parameter-sweep). |
| `migration_function`    | `Function`                   | Function describing if and how migration happens after reproduction. Built-in name (e.g. `random_migration`) or custom function. |
| `genotype_to_phenotype_mapping` | `Function`       | Function describing how phenotype is calculated from genotype. Default functions are defined for sexual reproduction. |
| `parameter_constraints` | `Function` | Function that checks the parameters with `@assert` before any simulation runs. Called once per parameter combination. Defaults to `nothing`. See [Setting parameter constraints](#setting-parameter-constraints). |

---

## List of parameters

These are the parameters already in place that you can use to control the simulations.

```
n_gen                     => Number of generations
patch_size                     => Initial number of individuals per patch
n_patch                   => Number of patches (groups)
n_loci                    => Number of loci (for diploid traits)
mu_m                      => Mutation rate per trait

str_selection             => Strength of selection (scale fitness)
log_from                   => First generation to record output
log_every                   => Interval between outputs
log_by                        => Data resolution: 'g', 'p', or 'i'
extras_names        => Custom names for extra variables returned by the fitness function; overrides field names if a NamedTuple is used
write_file                => Whether to write results to disk
name_model                => Prefix for output filename
parameters_to_omit        => Parameters excluded from filename
n_rep                   => Number of independent simulations
split_rep               => Whether to save each simulation replicate to a separate file. Requires split_sweep = true. Also controls whether simulation replicates can be parallelised independently.
sweep_grid                => Whether to use a full Cartesian product (`true`, default) or zip mode (`false`)
split_sweep               => Whether to save each parameter set to a separate file. Also controls whether parameter sets can be parallelised independently.
distributed               => Whether to run simulations on distributed workers. (requires @everywhere for functions and imports)
flatten_single_patch                  => Flatten population if there is a single patch
seed                     => Optional random seed. One integer for a single simulation, or a vector of seeds for multiple simulations. Generated randomly if omitted.
```



---

## Project Structure

- `src/` — core simulation logic split into mutation, reproduction, migration, simulation engine
- `test/` — unit tests
- `examples/basic_example.jl` — demo models

## License

MIT License. See `LICENSE` file for details.