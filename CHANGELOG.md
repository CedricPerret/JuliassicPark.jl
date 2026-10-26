# Changelog

## [0.3.0] - 2026-10-26
### Breaking
- Renamed built-in parameters: `n_print` → `log_from`, `j_print` → `log_every`, `de` → `log_by`, `simplify` → `flatten_single_patch`, `n_simul` → `n_rep`, `split_simul` → `split_rep`, `n_ini` → `patch_size`, `delta` → `allele_eff`, `other_output_names` → `extras_names`. Old names still work with a deprecation warning.
- Derived parameters are now given directly in `parameters` as `name = Derived(f)`. The `additional_parameters` keyword of `evol_model` is removed.
- Removed `additional_parameters_to_omit`. Derived parameters to hide go in `parameters_to_omit` (the old name still works with a warning).
- Sweeping over `Derived` values is not allowed. Sweep over a normal parameter read by the function instead.
- `z_ini` given as a table must be a `DataFrame` (not a view) containing a single simulation.
- `random_migration` and `random_migration!` only work with explicit reproduction, since they change patch sizes. An error is raised otherwise.
- `boundaries` no longer has a default value, including in `mutation!`.

### Added
- `parameter_constraints`: user function checking the parameters before running.

### Fixed
- Mutation rate set to 0 when `sigma_m` is 0, avoiding an error with zero-variance truncated distributions.
- Filename abbreviation works for vectors of non-numbers.
- In-place migration with `mig_rate = 0` no longer reverts to the previous population. In-place migration functions now return `false` when migration is skipped. Migration is skipped with a single patch.
- Continuing from a DataFrame: `n_patch` and `patch_size` are taken from the table, and a table without `i_simul` keeps the new run's ID.
- Sexual reproduction with several loci and the default mapping now runs (requires Boolean alleles).
- `reproduction_WF_island_model_hard_selection!` is now exported.
- A swept parameter not in `parameters` is now included in the filename, so sweeps no longer overwrite each other's file.
- Function to save output for the case `n_cst = false` and `log_by = i`.
- Simulation, in-place fitness function + @extras. Do not change output when it should not print, so it does not write NaN.

## [0.2.4] - 2026-10-20
- Change the sort parameters when building filename so it works without dependencies (dictionnary cannot be sort with the base julia)

## [0.1.3] - 2026-12-03
- Update simulations so that only the first element of output is updated when should_it_print = false. It avoids compatibility problem with the default values of extra, e.g. NaN instead of [NaN].

## [0.1.2] - 2026-05-03
- Fix fitness-type assertion for non in-place fitness functions by validating values (works even when the fitness container is Vector{Any} / Any[]).
- Fix trait-count inference in init_data_output by passing n_trait explicitly (instead of inferring via _leaf_type/fieldcount).
- Update _leaf_type documentation with a warning that recursion stops at Any when containers are type-erased (e.g. Vector{Any}).
- Fix make_output_corrector to split multi-trait phenotype outputs using n_traits, robust to Vector{Any} containers.
- Fix preprocess_fitness_function output standardisation to avoid type erasure (moving collect to occur after output reorganisation/inversion so heterogeneous per-call outputs don’t force Vector{Any}).

## [0.1.1] - 2026-03-03
- Cleaned up dependencies to reduce installation and precompilation time.
- No user-facing API changes; behavior should be identical to 0.1.0.

## [0.1.0] - 2026-02-09
- Initial public release of JuliassicPark.