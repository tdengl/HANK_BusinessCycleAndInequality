# __precompile__(false)
# Code runs on Julia 1.6
# ------------------------------------------------------------------------------
## Package Calls
# ------------------------------------------------------------------------------
# Packages used: Plots Distributions BenchmarkTools JLD2 FileIO DataFrames ForwardDiff
# SparseArrays LinearAlgebra Random LaTeXStrings MatrixEquations Roots KrylovKit JSON
# SpecialFunctions FFTW Parameters Setfield MCMCChains StatsPlots Optim CSV 
# OrderedCollections Flatten FieldMetadata MKL

module HANKEstim

if Sys.islinux() # issues encountered when using mkl with macos + more than 1 thread
    using MKL
end
if Sys.iswindows() # issues encountered when using mkl with macos + more than 1 thread
    using MKL
end
using Plots, Distributions, BenchmarkTools, JLD2, FileIO, DataFrames, ForwardDiff
using SparseArrays, LinearAlgebra, Random, LaTeXStrings, MatrixEquations
using Roots, KrylovKit, JSON; using SpecialFunctions: erf; using FFTW: dct
using Parameters, Setfield, MCMCChains, StatsPlots, Optim, CSV, OrderedCollections
using Flatten, FieldMetadata; import Flatten: flattenable

export LinearResults, linearize_full_model, EstimResults, find_mode, load_mode, montecarlo,
        find_steadystate, prepare_linearization, @writeXSS, @make_fn, @make_fnaggr,
        @make_struct, @make_struct_aggr, compute_steadystate, SteadyResults,
        mylinearinterpolate3, Ksupply, Fsys, SGU, Fsys_agg, SGU_estim, SolveDiffEq,
        mode_finding, likeli, kalman_filter, kalman_filter_smoother, measurement_error, rwmh,
        ModelParameters, NumericalParameters, EstimationSettings,
        Tauchen, EGM_policyupdate, Kdiff, distrSummaries, @generate_equations,
        @make_deriv, @make_deriv_estim, prioreval, FSYS_MIT_ALL_SHOCKS

include("1_Model/input_aggregate_names.jl")

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
include("1_Model/Parameters.jl")
include("3_NumericalBasics/Structs.jl")
include("6_Estimation/prior.jl")

e_set = EstimationSettings()
@make_struct IndexStruct
@make_struct_aggr IndexStructAggr

include("2_includeLists/include_NumericalBasics.jl")
include("2_includeLists/include_HetAgentsFcns.jl")
include("2_includeLists/include_LinearizationFunctions.jl")
include("2_includeLists/include_Estimation.jl")


@make_fn produce_indexes
@make_fnaggr produce_indexes_aggr


@doc raw"""
    compute_steadystate()

Compute steady state.

# Returns
`struct` `SteadyResults`, containing returns of [`find_SS()`](@ref)
"""
function compute_steadystate(m_par)
  #Calculate steady state capital stock
  KSS, VmSS, VkSS, distrSS, n_par, m_par = find_steadystate(m_par)

  # Prepare steadys state information for linearization
  XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, n_par, m_par,
      CDF_SS, CDF_m, CDF_k, CDF_y, distrSS = prepare_linearization(KSS, VmSS, VkSS, distrSS, n_par, m_par)
  
      println("Number of DCTs for Vm:")
      println(length(compressionIndexes[1]))

      println("Number of DCTs for Vk:")
      println(length(compressionIndexes[2]))

      println("Number of DCTs for COP:")
      println(length(compressionIndexes[3]))
      

  # Reduce further based on importance in dynamics at initial guess 
  if n_par.further_compress
    compressionIndexes, indexes, n_par = control_reduc_loop(compressionIndexes, XSS,
                                         m_par, n_par, indexes, distrSS,shock_names)
  end

    println("Number of DCTs for reduced Vm:")
    println(length(compressionIndexes[1]))

    println("Number of DCTs for reduced Vk:")
    println(length(compressionIndexes[2]))

    println("Number of DCTs for reduced COP:")
    println(length(compressionIndexes[3]))
    

  return SteadyResults(XSS, XSSaggr, indexes, indexes_aggr, compressionIndexes, 
                       n_par, m_par, CDF_SS, CDF_m, CDF_k, CDF_y, distrSS)
end

@doc raw"""
    linearize_full_model()

Linearize the full solution (i.e. with idiosyncratic states and controls) around the steady state,
using [`SGU()`](@ref).

# Returns
`struct` `LinearResults`, containing
- `A::Array{Float64,2}`,`B::Array{Float64,2}`: first derivatives of [`Fsys()`](@ref) with respect to arguments `X` [`B`]
    and `XPrime` [`A`]
- `State2Control::Array{Float64,2}`: observation equation
- `LOMstate::Array{Float64,2}`: state transition equation
"""
function linearize_full_model(sr::SteadyResults, m_par::ModelParameters)
    A = zeros(sr.n_par.ntotal, sr.n_par.ntotal)
    B = zeros(sr.n_par.ntotal, sr.n_par.ntotal)
    if sr.n_par.verbose
      println("Initial linearization")
    end
    State2Control, LOMstate, SolutionError, nk, A, B = SGU(sr.XSS, copy(A), copy(B), m_par, sr.n_par, sr.indexes, sr.compressionIndexes, sr.distrSS; estim = false);  
    #   indexes               = produce_indexes(sr.n_par, sr.compressionIndexes[1], sr.compressionIndexes[2], sr.compressionIndexes[3])
    #   indexes_aggr          = produce_indexes_aggr(sr.n_par)
    # @timev State2Control, LOMstate, SolutionError, nk, A, B = SGU_estim(sr.XSSaggr, copy(A), copy(B), sr.m_par, sr.n_par, sr.indexes, sr.indexes_aggr, sr.distrSS; estim = true)
    #BLAS.set_num_threads(1)
    #BLAS.set_num_threads(Threads.nthreads())
    
    return LinearResults(State2Control, LOMstate, A, B, SolutionError, nk)
end

@doc raw"""
    find_mode(sr,mr)

Find parameter that maximizes likelihood of data given linearized model `mr`.

# Arguments
- `sr::SteadyResults`
- `lr::LinearResults`

# Returns
`struct` `EstimResults`, containing all returns of [`mode_finding()`](@ref)
"""
function find_mode(sr::SteadyResults, lr::LinearResults, m_par::ModelParameters)
    if sr.n_par.verbose
      println("Started mode finding. This might take a while...")
    end
    par_final, hessian_final, posterior_mode, meas_error, meas_error_std, parnames, Data, Data_missing, H_sel, priors =
        mode_finding(sr.XSSaggr, lr.A, lr.B, sr.indexes, sr.indexes_aggr,sr.distrSS, sr.compressionIndexes, m_par, sr.n_par, e_set)
    if sr.n_par.verbose
        println("Mode finding finished.")
    end
    return EstimResults(par_final, hessian_final, meas_error, meas_error_std, parnames, Data, Data_missing, H_sel, priors)
end

function load_mode(sr::SteadyResults;file::String = e_set.mode_start_file)
    @load file par_final hessian_final meas_error meas_error_std parnames Data Data_missing H_sel priors

    # Load data
    Data_temp = DataFrame(CSV.File(e_set.data_file; missingstring = "NaN"))
    data_names_temp = propertynames(Data_temp)
    for i in data_names_temp
      name_temp = get(e_set.data_rename, i, :none)
      if name_temp != :none
        rename!(Data_temp, Dict(i=>name_temp))
      end
    end

    observed_vars = e_set.observed_vars_input
    Data = Matrix(Data_temp[:, observed_vars])
    Data_missing = ismissing.(Data)
    nobservables = size(Data)[2]

    H_sel = zeros(nobservables, sr.n_par.nstates + sr.n_par.ncontrols)
    for i in eachindex(observed_vars)
      H_sel[i,   getfield(sr.indexes,(observed_vars[i]))]  = 1.0
    end

    return EstimResults(par_final, hessian_final, meas_error, meas_error_std, parnames, Data, Data_missing, H_sel, priors)
end

@doc raw"""
    montecarlo(mr,er;file=e_set.save_posterior_file)

Sample posterior of parameter vector with [`rwmh()`](@ref), take sample mean as
parameter estimate, and save all results in `file`.

# Arguments
- `sr::SteadyResults`
- `mr::LinearResults`
- `er::EstimResults`
"""
function montecarlo(sr::SteadyResults,lr::LinearResults, er::EstimResults, m_par::ModelParameters; file::String = e_set.save_posterior_file)
    hessian_sym = Symmetric(nearest_spd(inv(er.hessian_final)))
    if sr.n_par.verbose
        println("Started MCMC. This might take a while...")
    end
    if e_set.multi_chain_init == true
        init_draw, init_success, init_iter = multi_chain_init(er.par_final, hessian_sym, sr.n_par, er.Data, er.Data_missing,
        er.H_sel, sr.XSSaggr, lr.A, lr.B, sr.indexes,sr.indexes_aggr, m_par, e_set,
        sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std)

        par_final = init_draw
        if init_success == false
            error("Couldn't find initial value that produces posterior")
        end
    else 
        par_final = copy(er.par_final)
    end

    draws_raw, posterior, accept_rate = rwmh(par_final, hessian_sym, sr.n_par, er.Data, er.Data_missing,
    er.H_sel, sr.XSSaggr,lr.A, lr.B, sr.indexes,sr.indexes_aggr, m_par, e_set, sr.distrSS, sr.compressionIndexes,
    er.priors, er.meas_error, er.meas_error_std)

    draws = draws_raw[e_set.burnin+1:end, :]

    ##
    parnames_ascii = collect(metaflatten(m_par, label))
    if e_set.me_treatment != :fixed
    for i in eachindex(e_set.meas_error_input)
        push!(parnames_ascii, string("sigma_me_", e_set.meas_error_input[i]))
    end
    end

    chn = Chains(reshape(draws, (size(draws)...,1)), [string(parnames_ascii[i]) for i = 1:length(parnames_ascii)])
    chn_summary = summarize(chn)
    par_final = chn_summary[:,:mean];

    ##
    if e_set.me_treatment != :fixed
    m_par = Flatten.reconstruct(m_par, par_final[1:length(par_final) - length(er.meas_error)])
    else
    m_par = Flatten.reconstruct(m_par, par_final)
    end
    State2Control, LOMstate, alarm_sgu, nk = SGU_estim(sr.XSSaggr,  lr.A, lr.B, m_par, sr.n_par, sr.indexes, sr.indexes_aggr, sr.distrSS; estim = true)

    smoother_output = likeli(par_final, er.Data, er.Data_missing, er.H_sel, sr.XSSaggr, lr.A, lr.B, sr.indexes, sr.indexes_aggr, m_par, sr.n_par, e_set, 
                             sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std; smoother = true)
    
    if sr.n_par.verbose
        println("MCMC finished.")
    end
    ##
    # @save file draws_raw accept_rate aggr_names_ascii aggr_names state_names_ascii state_names control_names control_names_ascii e_set er.Data State2Control LOMstate er.parnames sr.indexes sr.indexes_aggr par_final m_par lr.A lr.B hessian_sym sr.n_par draws posterior er.hessian_final er.meas_error er.meas_error_std er.Data_missing er.H_sel er.priors smoother_output
    save(file,Dict(
                    "draws_raw" => draws_raw,
                    "accept_rate" => accept_rate,
                    "aggr_names_ascii" => aggr_names_ascii,
                    "aggr_names" => aggr_names,
                    "state_names_ascii" => state_names_ascii,
                    "state_names" => state_names,
                    "control_names" => control_names,
                    "control_names_ascii" => control_names_ascii,
                    "e_set" => e_set,
                    "Data" => er.Data,
                    "State2Control" => State2Control,
                    "LOMstate" => LOMstate,
                    "parnames" => er.parnames,
                    "indexes" => sr.indexes,
                    "indexes_aggr" => sr.indexes_aggr,
                    "par_final" => par_final,
                    "m_par" => m_par,
                    "A" => lr.A,
                    "B" => lr.B,
                    "hessian_sym" => hessian_sym,
                    "n_par" => sr.n_par,
                    "draws" => draws,
                    "posterior" => posterior,
                    "hessian_final" => er.hessian_final,
                    "meas_error" => er.meas_error,
                    "meas_error_std" => er.meas_error_std,
                    "Data_missing" => er.Data_missing,
                    "H_sel" => er.H_sel,
                    "priors" => er.priors,
                    "smoother_output" => smoother_output
    ))
end

end # module HANKEstim