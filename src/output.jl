#***********************************************
#*** Output and I/O functions
#***********************************************

#-----------------------------------------------
#*** Name of the file (contains parameters)
#-----------------------------------------------

"""
    _abbreviate(x)

Converts a number, string, or distribution into a short string. This is useful to avoid unreadable long name file.

- Numbers ≥ 1000 are shown as `"1.0k"`, others are rounded.
- Filenames like `"path/file.csv"` become `"file"`.
- Distributions are encoded as `"Name_param1_param2"`.

Used for output filenames or compact parameter summaries.
"""
function _abbreviate(x)
    if x isa Real
        if x >= 1000
            return(string((x/1000))*"k")
        elseif (x%10==0)
            return(string(Int(x)))
        else
            return(string(round(x,digits=5)))
        end
    elseif x isa Vector
        return string(round.(x,digits=5))
    elseif x isa String
        if any(occursin.([".csv",".txt"],x))
            #findlast return nothing. We use something to change it to 0
            return(SubString(x,something(findlast('/',x),0)+1,findlast('.',x)-1))
        else
            return(string(x))
        end
    elseif x isa Distribution
        string(nameof(typeof(x)))*"_"*join([string(getfield(x, parameter_distrib)) for parameter_distrib in fieldnames(typeof(x))], "_")
    else
        return(string(x))
    end
end

#--- Internal function: builds base filename string from parameter values
function _build_filename(wd::String, parameters::Dict, parameters_to_omit::Vector{Symbol}, format::String; swept=Dict{Symbol,Vector}())
    parameters_copy=deepcopy(parameters)
    #--- Remove parameters we do not want to save
    for i in parameters_to_omit
        delete!(parameters_copy,i)
    end
    #--- Replace the parameter we vary by first -> last value
    for (key, sweep_vals) in swept
        if haskey(parameters_copy, key)
            parameters_copy[key] = string(sweep_vals[1], "→", sweep_vals[end])
        end
    end
    #--- Generate name of the file with -key=value for each parameter (sorted for reproducibility), with abbreviated values
    name_file = join(["-" * string(k) * "=" * _abbreviate(v) for (k, v) in sort(collect(parameters_copy); by=first)], "")
end

#--- Public function: returns full file path
function _build_filepath(wd::String, parameters::Dict, parameters_to_omit::Vector{Symbol}, format::String; swept=Dict{Symbol,Vector}())
    name_file = _build_filename(wd, parameters, parameters_to_omit, format; swept=swept)
    return wd * name_file * format
end

#--- Variant for split_simul=true: appends simulation ID
function _build_filepath(wd::String,parameters::Dict,parameters_to_omit::Array{Symbol,1},format::String,i_simul::Int64; swept=Dict{Symbol,Vector}())
    name_file = _build_filename(wd, parameters, parameters_to_omit, format; swept=swept)
    return wd * name_file * "-S="*string(i_simul)*format
end

#-----------------------------------------------
#*** Saving function
#-----------------------------------------------

"""
A macro for optionally computing and printing "extra" values inside a fitness function.

The trick is that all variables assigned in the block are first initialized to `[]`, so that they are always defined even if not computed. The expressions in the block are only executed if `should_it_print == true`.

!!! warning
    This macro will overwrite any existing variables of the same name with `[]` if they are reassigned in the block. Avoid name collisions with earlier code.
"""
macro extras(block)
    # Extract all top-level expressions (may include LineNumberNodes)
    stmts = block isa Expr && block.head == :block ? block.args : [block]
    # Keep only assignment expressions like `x = ...`
    assignments = [stmt for stmt in stmts if stmt isa Expr && stmt.head == :(=)]
    # Extract the left-hand side variable names from assignmentss
    vars = [stmt.args[1] for stmt in assignments]
    # Generate code to define each variable as `NaN` 
    # Using NaN here is not ideal, since NaN is normally reserved for invalid numerical results and it propagates through arithmetic. 
    #However, it avoids compatibility issues (such as Float64 values being promoted to Any and causing errors in the output). 
    # + this is acceptable, because these values are never used in further computations — they only serve as placeholders so that the user can always return a single, consistent tuple.
    init = [:( $(esc(var)) = NaN ) for var in vars]   
    #init = [:( $(esc(var)) = nothing ) for var in vars]
    #init = [:( $(esc(var)) = [] ) for var in vars]
    # Return full quoted expression: define all vars, then conditionally assign
    return quote
        $(init...)   # Always define all variables as `nothing`
        if $(esc(:should_it_print))
            $(esc.(stmts)...) # Only evaluate assignments if needed
        end
    end
end



"""
    _infer_variable_resolution(output::Vector; n_patch::Int = 1)

Identifies the resolution of the variable output in the evolutionary model, that is whether it is measured at the generation, group (patch), or individual level.

Used internally to determine how simulation outputs should be reshaped or aggregated before saving.

# Arguments
- `output`: A vector containing outputs from a model run. Each element may correspond to data at a different level.
- `n_patch`: Number of patches in the model (default is 1). Used to distinguish between patch-level and individual-level structures.

# Returns
A `Vector{Char}` of the same length as `output`, where each character indicates the type of structure:
- `'G'`: Generation-level
- `'g'`: Generation-level (encapsulated in a vector)
- `'p'`: Patch-level
- `'i'`: Individual-level (in a nested vector = as a structured population)
- `'I'`: Individual-level (in a single vector = as a well-mixed population)

These types account for different ways users may return outputs from a model. For example:
- `'G'` handles scalar values instead of `[x]`
- `'I'` handles flat vectors when there's a single patch or multiple traits
- Multiple traits can be encoded as tuples inside these structures

# Example
test = [1, [1, 2], [[1, 2, 3], [4, 5, 6]], [1, 2, 3, 4, 5, 6]]

_infer_variable_resolution(test, n_patch=2)
# → ['G', 'p', 'i', 'I']

_infer_variable_resolution(test)
# → ['G', 'g', 'i', 'I']
"""

function _infer_variable_resolution(output;n_patch=1)
    level_of_detail = Vector{Char}(undef, length(output))
        for (i,o) in enumerate(output)
            if !isa(o, AbstractVector)
                #->{Any}
                level_of_detail[i] = 'G'
            elseif isa(o[1],Union{AbstractVector, Tuple})
                #->{Vector{Vector}}
                level_of_detail[i] = 'i'
            else
                if length(o) == n_patch
                    #-> {Vector}[n_patch]
                    level_of_detail[i] = 'p'
                elseif length(o) == 1
                    #->{Vector}[1]
                    level_of_detail[i] = 'g'
                else
                    #->{Vector}[n_pop * n_patch]
                    level_of_detail[i] = 'I'
                end
            end
        end
    return(level_of_detail)
end

"""
    init_data_output(de, output_names, output_example, n_gen, gen_first_print, print_every,
                     i_simul, n_patch, n_pop, n_cst;
                     output_cst_names=[], output_cst=[]) -> (DataFrame, Function)

Create a `DataFrame` to store model outputs and return a saving function to fill it during a simulation.

Each element of the output is classified by its structure (generation-, patch-, or individual-level), and transformed to match the desired level of detail specified by the user in `de`.

If the detail level (`de`) is higher than the level of the variable (e.g., trying to print a scalar per individual), the output is repeated. If it is lower (e.g., a per-individual value printed per generation), it is averaged.

If `n_cst = true`, the population size is assumed constant, and the DataFrame is preallocated for faster writing. If `n_cst = false`, data is appended dynamically each time output is saved.

Trait values are extracted from the first element of `output_example`, and if they are stored as tuples, one column is created per trait (`z1`, `z2`, etc.). Other columns are named using `output_names`, or default names like `V1`, `V2`, etc.

# Arguments
- `de::Char`: Desired level of data resolution (`'g'` for generation, `'p'` for patch, `'i'` for individual).
- `output_names::Vector{String}`: Names of the outputs. If too short, default names are added.
- `output_example::Vector`: Example output used to infer structure and types.
- `n_gen::Int`: Total number of generations.
- `gen_first_print::Int`: First generation at which output is saved.
- `print_every::Int`: Number of generations between two saves.
- `i_simul`: Simulation ID.
- `n_patch::Int`: Number of patches in the population.
- `n_pop::Int`: Number of individuals per patch.
- `n_cst::Bool`: Whether the population size is constant.
- `output_cst_names::Vector{String}`: Names of constant outputs (optional).
- `output_cst::Vector`: Constant outputs for each row (optional).

# Returns
- A `DataFrame` ready to store simulation results.
- A function `save_data_to_df(df, i_gen, output)` that writes outputs at the correct position in the DataFrame.

# Notes
- Generation = one call to the reproduction function (Moran or Wright–Fisher logic).
- Constant outputs (e.g. `ID`, `K`, or parameters) are repeated automatically depending on `de` and `n_cst`.

# Output types (inferred from `output_example`)
- `'G'`: Scalar → generation-level
- `'g'`: Vector of length 1 → generation-level
- `'p'`: Vector of length `n_patch` → patch-level
- `'i'`: Vector of `n_patch` vectors of length `n_pop` → individual-level (structured)
- `'I'`: Vector of length `n_patch * n_pop` → individual-level (well-mixed)
"""
function init_data_output(de,output_names::Vector{String},output_example,n_gen::Int, gen_first_print::Int, print_every::Int, i_simul, n_patch::Int, n_pop::Int, n_cst::Bool; output_cst_names=[],output_cst=[], n_traits)
    #--- Safety checks
    @assert de in ['g', 'p', 'i'] "Invalid value for `de`: expected one of 'g', 'p', or 'i'."
    @assert n_patch ≥ 1 "Expected `n_patch` to be a positive integer."
    @assert n_pop ≥ 1 "Expected `n_pop` to be a positive integer."
    @assert n_gen ≥ 1 "Expected `n_gen` to be a positive integer."
    @assert gen_first_print ≥ 1 "Expected `gen_first_print` to be a positive integer."
    @assert print_every ≥ 1 "Expected `print_every` to be a positive integer."
    @assert gen_first_print ≤ n_gen "First print generation must be ≤ total number of generations."
    for i in 1:length(output_example)
            @assert length(output_example[i]) == 1 || length(output_example[i]) == n_patch || (length(output_example[i]) == n_pop && !(isa(output_example[i][1], Vector))) || (isa(output_example[i][1], Vector) && length(output_example[i][1]) == n_pop) "The variable to save `:$(Symbol(output_names[i]))` must have resolution compatible with the model: " *
            "either one value per generation, per patch (length = n_patch), or per individual " *
            "(length = n_patch, with sub-vectors of length n_ini)."
    end

    #--- Number of generation to print 
    gen_printed = gen_first_print:print_every:n_gen
    n_gen_printed = length(gen_printed)
    #--- Adjust output names if not enough names are provided
    if length(output_names) < length(output_example) 
        output_names = [output_names;["V"*string(i) for i in 1:(length(output_example) -  length(output_names))]]
    end
    #--- Multiple traits?
    ## If yes, the first vector needs to be decoupled into one vector per trait.
    #@show output_example
    correct_output_for_n_trait = _make_output_corrector(output_example[1],n_traits)
    if n_traits > 1
        output_names = [["z"*string(i) for i in 1:n_traits]; output_names[2:end]]
    end
    #--- Identify the resolution of the variables measured
    level_output = _infer_variable_resolution(correct_output_for_n_trait(output_example),n_patch=n_patch)
    type_output_cst = _infer_variable_resolution(output_cst,n_patch=n_patch)

    ## Each output is modified (correction_function) to fit the resolution required by the user.
    #*** Results required at the generation level
    if de == 'g'
        #--- Create the vector of functions to apply to each output 
        correction_function = replace(level_output,
        'g'=>(x->x[1]),
        'G'=>(x->x),
        #Case where we calculate for instance variance in a patch and then do the mean of that. 
        #I don't see how we will get NaN at individual level for now so I did not implemented there.
        #could be moved to n_cst = false
        'p'=>(x->mean(filter(!isnan, x))),
        'i'=>(x->mean(vcat(x...))),
        'I'=>(x->mean(x)))
        ## with corresponding names
        output_names=string.(replace(level_output,'g'=>"",'G'=>"",'p'=>"mean_",'i'=>"global_mean_",'I'=>"global_mean_"),output_names)
        #! We do not print the lower levels of details of the constant output
        cst_to_keep =  type_output_cst .∈ Ref(['g', 'G'])
        output_cst, output_cst_names, type_output_cst = getindex.((output_cst, output_cst_names, type_output_cst), Ref(cst_to_keep))
        empty_output = [i[] for i in _leaf_type.([correction_function[i](correct_output_for_n_trait(output_example)[i]) for i in eachindex(correction_function)])]


        #--- Initialise dataFrame
        if n_cst
            df_res = DataFrame(i_simul=repeat([i_simul], inner=n_gen_printed), gen=gen_printed)
            ## Add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,n_gen_printed)),
            'G'=>(x->repeat([x],n_gen_printed)))
            for i in eachindex(output_cst_names)
                #! We do not output the lower levels.
                df_res[!, output_cst_names[i]] = correction_function_cst[i](output_cst[i])
            end
            last_column = ncol(df_res)
            ## Add the column of variables
            for i in 1:length(output_names)
                T = eltype(empty_output[i])
                df_res[!, output_names[i]] = Vector{T}(undef, n_gen_printed)
            end
            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen -> floor(Int, (i_gen - gen_first_print) / print_every) + 1

        else
            ## Create empty dataframe of the right type
            # #--- The type of output is either float if not at the level of generation (since we do the mean) or given by the output_example
            # type_output = [i[] for i in  _leaf_type.(correct_output_for_n_trait(output_example))]
            # for (i, lvl) in enumerate(level_output)
            #     if lvl in ['i', 'I', 'p']
            #         type_output[i] = Float64[]
            #     end
            # end
            df_res = DataFrame([repeat([Int64[]],2);[i[] for i in  _leaf_type.(output_cst)];empty_output],
            ["i_simul","gen",output_cst_names...,output_names...]
            )
            ## Create the function to add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->x),
            'G'=>(x->[x]))
            get_base_output = (i_gen, n_pop) -> [[i_simul], [i_gen],[correction_function_cst[i](output_cst[i]) for i in eachindex(output_cst_names)]...]
        end
    #*** Results at the PATCH level
    elseif de == 'p'
        #--- Create the vector of functions to apply to each output 
        #@Iterators.partition faster than list comprehension
        correction_function = replace(level_output,
        'g'=>(x->fill(x[1],n_patch)),
        'G'=>(x->fill(x,n_patch)),
        'p'=>(x->x),
        'i'=>(x->mean.(x)),
        ## Apply only to case where patch = 1 but we keep the format of vector of vector
        'I'=>(x->[mean(x)]))
        # 'I'=>(x->mean.(collect(Iterators.partition(x,n_pop)))))
        ## with corresponding names
        output_names = string.(replace(level_output,'g'=>"",'G'=>"",'p'=>"",'i'=>"mean_",'I'=>"mean_"),output_names)
        #! We do not print the lower levels of details of the constant output
        cst_to_keep =  type_output_cst .∈ Ref(['g', 'G','p'])
        output_cst, output_cst_names, type_output_cst = getindex.((output_cst, output_cst_names, type_output_cst), Ref(cst_to_keep))
        empty_output = [i[] for i in _leaf_type.([correction_function[i](correct_output_for_n_trait(output_example)[i]) for i in eachindex(correction_function)])]
        #--- Initialise dataFrame
        if n_cst
            total_length_output = n_gen_printed * n_patch
            df_res = DataFrame(i_simul=repeat([i_simul], inner = total_length_output),
                gen=repeat(gen_printed, inner=n_patch),
                patch=repeat(1:n_patch, outer=n_gen_printed))

            ## Add the constant output. 
            # If population size is constant, we can write down all the constant output now.
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,total_length_output)),
            'G'=>(x->repeat([x],total_length_output)),
            'p'=>(x->repeat(x,outer=n_gen_printed)))
            for i in eachindex(output_cst_names)
                df_res[!, output_cst_names[i]] = correction_function_cst[i](output_cst[i])
            end
            last_column = ncol(df_res)
            ## Add the variables columns 
            for i in 1:length(output_names)
                T = eltype(empty_output[i])
                df_res[!, output_names[i]] = Vector{T}(undef, total_length_output)
            end
            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen -> (n_patch*(floor(Int,(i_gen-gen_first_print)/print_every))+1):n_patch*(1+floor(Int,(i_gen-gen_first_print)/print_every))
        else
            #-> Population size is not constant so we need to write down the constant output each time we want to save
            ##Create empty dataframe of the right type
            # #--- The type of output is either float if not at the level of generation (since we do the mean) or given by the output_example
            # type_output = [i[] for i in  _leaf_type.(correct_output_for_n_trait(output_example))]
            # for (i, lvl) in enumerate(level_output)
            #     if lvl in ['i', 'I']
            #         type_output[i] = Float64[]
            #     end
            # end
            df_res = DataFrame([repeat([Int64[]],3);[i[] for i in  _leaf_type.(output_cst)];empty_output],
            ["i_simul","gen","patch",output_cst_names...,output_names...])
            ## Create the function to add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,n_patch)),
            'G'=>(x->repeat([x],n_patch)),
            'p'=>(x->x))
    
            get_base_output = (i_gen, n_pop) -> [fill(i_simul, n_patch),
            fill(i_gen, n_patch),
            collect(1:n_patch),
            [correction_function_cst[i](output_cst[i]) for i in eachindex(output_cst_names)]...]
        end
    #*** Print variables at the level of the INDIVIDUAL 
    elseif de == 'i'
        #--- Create the vector of functions to apply to each output 
        correction_function = replace(level_output,
        'g'=>(x->fill(x[1],n_pop*n_patch)),
        'G'=>(x->fill(x,n_pop*n_patch)),
        'p'=>(x->repeat(x,inner=n_pop)),
        'i'=>(x-> reduce(vcat,x)),
        'I'=>(x-> x))
        empty_output = [i[] for i in _leaf_type.([correction_function[i](correct_output_for_n_trait(output_example)[i]) for i in eachindex(correction_function)])]
        #--- Initialise dataFrame
        if n_cst
            total_length_output = n_gen_printed * n_patch * n_pop
            df_res = DataFrame(i_simul=repeat([i_simul], inner=total_length_output),
                gen=repeat(gen_printed, inner=n_patch * n_pop),
                patch=repeat(1:n_patch, outer=n_gen_printed, inner=n_pop),
                ind=repeat(1:n_pop*n_patch, outer= n_gen_printed))

            ## Add the constant output. 
            correction_function_cst=replace(type_output_cst,
            'g'=>(x->repeat(x,total_length_output)),
            'G'=>(x->repeat([x],total_length_output)),
            'p'=>(x->repeat(x,outer=n_gen_printed, inner=n_pop)),
            'i'=>(x->repeat(reduce(vcat,x),outer=n_gen_printed)),
            'I'=>(x->repeat(x,outer=n_gen_printed)))
            corrected_output_cst = [correction_function_cst[i](output_cst[i]) for i in eachindex(correction_function_cst)]
            for i in eachindex(output_cst_names)
                df_res[!, output_cst_names[i]] = corrected_output_cst[i]
            end

            last_column = ncol(df_res)
            ## Add the variables columns 
            for i in 1:length(output_names)
                T = eltype(empty_output[i])
                df_res[!, output_names[i]] = Vector{T}(undef, total_length_output)
            end

            ## Create the function to calculate the position of the output
            calculate_position_output = i_gen ->(n_patch*n_pop*(floor(Int,(i_gen-gen_first_print)/print_every))+1):n_patch*n_pop*(1+floor(Int,(i_gen-gen_first_print)/print_every))
        else
            #-> Varying population size
            #-> This is the only case where the size of the output is not constant so we need to take group sizes as input
            #--- Create the vector of functions to apply to each output 
            correction_function = replace(level_output,
            'g'=>((x, group_sizes)->fill(x[1],sum(group_sizes))),
            'G'=>((x, group_sizes)->fill(x,sum(group_sizes))),
            'p'=>((x, group_sizes)-> collect(Iterators.flatten((fill(x[i],group_sizes[i]) for i in 1:n_patch)))),
            'i'=>((x, group_sizes)-> reduce(vcat,x)),
            'I'=>((x, group_sizes)-> x))
            ##Create empty dataframe of the right type
            df_res = DataFrame([repeat([Int64[]],4);[i[] for i in  _leaf_type.(output_cst)];empty_output],
            ["i_simul","gen","patch","ind",output_cst_names...,output_names...]
            )
            ## Create the function to add the constant output. 
            correction_function_cst= correction_function
            get_base_output = (i_gen, n_pop) -> [fill(i_simul,  sum(n_pop)),
            fill(i_gen, sum(n_pop)),
            collect(Iterators.flatten((fill(i, n_pop[i]) for i in 1:length(n_pop)))),
            collect(Iterators.flatten((collect(1:n_pop[i]) for i in 1:length(n_pop)))),
            [correction_function_cst[i](output_cst[i]) for i in eachindex(output_cst_names)]...]
        end
    end
    #*** Create saving function

    if n_cst
        save_data_to_df = function (df, i_gen, output)
            
            if should_it_print(i_gen,gen_first_print,print_every) == true
                corrected_output_for_n_trait = correct_output_for_n_trait(output)
                corrected_output = [correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)]
                #@show corrected_output_for_n_trait
                for i in eachindex(corrected_output)
                    df[calculate_position_output(i_gen), last_column+i] = corrected_output[i]
                end
            end
        end
    else
        save_data_to_df = function (df, i_gen, output)
            if should_it_print(i_gen,gen_first_print,print_every) == true
                corrected_output_for_n_trait = correct_output_for_n_trait(output)
                #@show corrected_output_for_n_trait
                if isempty(output[1])
                    return
                end
                ## Here we add a single row. For append, we need a vector while for the n_cst = true, scalar are ok.
                group_sizes = length.(output[2]) 

                if de == 'g'
                    corrected_output = map(x->[x],[correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)])
                elseif de == 'p'
                    corrected_output = [correction_function[i](corrected_output_for_n_trait[i]) for i in eachindex(correction_function)]
                elseif de == 'i'
                    #-> Size of output varies
                    corrected_output = [correction_function[i](corrected_output_for_n_trait[i], group_sizes) for i in eachindex(correction_function)]
                end
                # @show output[1]
                # @show group_sizes
                # @show get_base_output(i_gen,group_sizes)
                # @show corrected_output
                append!(df,
                rename(DataFrame([get_base_output(i_gen,group_sizes);
                corrected_output],:auto),
                names(df))
                )
            end

        end
    end
    return(df_res, save_data_to_df)
end

should_it_print(i_gen,n_print,j_print) = (i_gen-n_print) % j_print == 0 && i_gen >= n_print

function _filter_cst_for_de(type_output_cst::Vector{Char}, de::Char)
    drop = de == 'g' ? ['i','I','p'] :
           de == 'p' ? ['i','I']     :
                       Char[]
    keep = .!(type_output_cst .∈ Ref(drop))
    return type_output_cst[keep]
end


"""
    _make_output_corrector(first_output::AbstractVector,n_traits:Int)

Return a function `correct(output)::Vector` that:
- splits the first element (traits) into one vector per trait if the trait type is a Tuple,
- otherwise returns the output unchanged.

Examples
--------
# population: Vector{Tuple}
correct = _make_output_corrector([(0.1,0.2), (0.3,0.4)],2)
correct([[ (0.1,0.2), (0.3,0.4) ], :foo])  # => [ [0.1,0.3], [0.2,0.4], :foo ]

# metapopulation: Vector{Vector{Tuple}}
correct = _make_output_corrector([[(0.1,0.2)], [(0.3,0.4)]],2)
correct([[ [(0.1,0.2)], [(0.3,0.4)] ], :foo])  # => [ [[0.1],[0.3]], [[0.2],[0.4]], :foo ]
"""
function _make_output_corrector(z::AbstractVector, n_traits::Int)
    n_traits == 1 && return identity
    # value-based: treat as metapop if first element is a vector
    is_metapop = !isempty(z) && first(z) isa AbstractVector

    if is_metapop
        #-> metapopulation: Vector{Vector{Tuple}}
        return function (output)
            z = output[1]
            per_trait = [map(p -> getfield.(p, i), z) for i in 1:n_traits]
            [per_trait; output[2:end]]
        end
    else
        #-> population: Vector{Tuple}
        return function (output)
            z = output[1]
            per_trait = [getfield.(z, i) for i in 1:n_traits]
            [per_trait; output[2:end]]
        end
    end
end
