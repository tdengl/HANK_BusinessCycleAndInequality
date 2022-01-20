# __precompile__(false)
# Code runs on Julia 1.6
module HANKrank
using Plots, Distributions, BenchmarkTools, JLD2, FileIO, DataFrames, ForwardDiff
using SparseArrays, LinearAlgebra, Random, LaTeXStrings, MatrixEquations
using KrylovKit; using SpecialFunctions: erf; using FFTW: dct
using Parameters, Setfield, MCMCChains, StatsPlots, Optim, CSV, OrderedCollections
using Flatten, FieldMetadata; import Flatten: flattenable
using JSON, Roots

# export LinearResults, linearize_full_model, EstimResults, find_mode, load_mode, montecarlo,
#         find_steadystate, prepare_linearization, @writeXSS, @make_fn, @make_fnaggr,
#         @make_struct, @make_struct_aggr, compute_steadystate, SteadyResults,
#         mylinearinterpolate3, Ksupply, Fsys, SGU, Fsys_agg, SGU_estim, SolveDiffEq,
#         mode_finding, likeli, kalman_filter, kalman_filter_smoother, measurement_error, rwmh,
#         ModelParameters, NumericalParameters, EstimationSettings,
#         Tauchen, EGM_policyupdate, Kdiff, distrSummaries, @generate_equations,
#         @make_deriv, @make_deriv_estim, prioreval, FSYS_MIT_ALL_SHOCKS
include("8_Rank/input_aggregate_names_rank.jl")
#include("8_rank/PreprocessInputs_rank.jl")

# ------------------------------------------------------------------------------
## Define Functions
# ------------------------------------------------------------------------------
#include("1_Model/Parameters.jl")
include("8_Rank/Parameters_rank.jl")
include("3_NumericalBasics/Structs.jl")
include("6_Estimation/prior.jl")

e_set = EstimationSettings()
@make_struct_aggr IndexStructAggr
@make_struct_aggr IndexStruct # will not be used, just defined here because FSYS_agg_generated needs to know it

include("3_NumericalBasics/BrentsMethod.jl")
include("3_NumericalBasics/LinInterpols.jl")
include("3_NumericalBasics/LocateFcn.jl")
include("3_NumericalBasics/GCintegration.jl")
include("3_NumericalBasics/DCT.jl")
include("3_NumericalBasics/FastRoot.jl")
include("3_NumericalBasics/MarkovChain.jl")
include("3_NumericalBasics/Schur_and_DualUtils.jl")

include("8_Rank/MacroUtils_rank.jl")
include("4_HetAgentsFcns/basic_fcns.jl")
include("8_Rank/prepare_linearization_generated_rank.jl")
include("5_LinearizationFunctions/SolveDiffEq.jl")
include("5_LinearizationFunctions/Shuffle.jl")
include("5_LinearizationFunctions/SGU_estim.jl")
include("8_Rank/fcn_incomes_rank.jl")
include("8_Rank/FSYS_agg_generated.jl")
include("8_Rank/SGU_rank.jl")
include("2_includeLists/include_Estimation.jl")

@make_fnaggr produce_indexes_aggr # macro makefn_aggr writes function produce_indexes_aggr which fills instance of IndexStructAggr with correct values of indexes
@make_fn produce_indexes # macro makefn writes function produce_indexes which fills instance of IndexStruct with correct values of indexes (no difference to indexes_aggr in RANK)

# # m_par = ModelParameters(δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, 
# # θ_Y = 0.176, γ_B = 0.4,  γ_π = -1.05, γ_Y = -0.852)
# # m_par = @set m_par.RB = 1.0./m_par.β # Rank interest rate on bonds is simply 1/beta
# # m_par = ModelParameters()
# # m_par = @set m_par.RB = 1.0./m_par.β
m_par = ModelParameters()
@load "7_Saves/parameter_example.jld2" par_final
m_par = Flatten.reconstruct(m_par, par_final[1:length(par_final)-length(e_set.meas_error_input)])
@set! m_par.RB = 1.0./m_par.β # Rank interest rate on bonds is simply 1/beta

@doc raw"""
    compute_steadystate()

Compute steady state.

# Returns
`struct` `SteadyResults`, containing returns of [`find_SS()`](@ref)
"""
function compute_steadystate(m_par)
  # Prepare steady state information for linearization
  XSS, XSSaggr, indexes, indexes_aggr, n_par, m_par  = prepare_linearization_rank(m_par)
  
  return SteadyResults(XSS, XSSaggr, indexes, indexes_aggr, zeros(1,1), n_par, m_par, zeros(1,1), zeros(1,1), zeros(1,1), zeros(1,1), zeros(1,1))
end

# sr = compute_steadystate(m_par)
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
  if sr.n_par.verbose
      println("Initial linearization")
  end
    State2Control, LOMstate, SolutionError, nk, A, B = SGU_rank(sr.XSS,m_par, sr.n_par, sr.indexes, sr.indexes_aggr; estim = false);  
    # @timev State2Control, LOMstate, SolutionError, nk, A, B = SGU_estim(sr.XSSaggr, copy(A), copy(B), sr.m_par, sr.n_par, sr.indexes, sr.indexes_aggr, sr.distrSS; estim = true)
    return LinearResults(State2Control, LOMstate, A, B, SolutionError, nk)
end
#lr = linearize_full_model(sr, m_par)
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
      if name_temp == :none

      else
        rename!(Data_temp, Dict(i=>name_temp))
      end
    end

    # reorder variables so that first-differences come first
    observed_vars = [e_set.observed_vars_input[e_set.growth_rate_select]; e_set.observed_vars_input[.!e_set.growth_rate_select]]
    Data = Matrix(Data_temp[:, observed_vars])
    Data_missing = ismissing.(Data)
    nobservables = size(Data)[2]

    n_growth_rates = sum(e_set.growth_rate_select) # number of growth rate variables
    if e_set.fd_flag == true
      H_sel = zeros(nobservables, 2 * (sr.n_par.nstates + sr.n_par.ncontrols))
      for i in eachindex(observed_vars)
        if i <= n_growth_rates
          H_sel[i,  getfield(indexes_,(observed_vars[i]))]  = 1.0
        else
          H_sel[i, sr.n_par.nstates + sr.n_par.ncontrols + getfield(sr.indexes,(observed_vars[i]))]  = 1.0
        end
      end
    else
      H_sel = zeros(nobservables, sr.n_par.nstates + sr.n_par.ncontrols)
      for i in eachindex(observed_vars)
        H_sel[i,   getfield(sr.indexes,(observed_vars[i]))]  = 1.0
      end
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

    if e_set.multi_chain_init == true
        init_draw, init_success, init_iter = multi_chain_init(er.par_final, hessian_sym, sr.n_par, er.Data, er.Data_missing,
        er.H_sel, sr.XSSaggr, lr.A, lr.B, sr.indexes,sr.indexes_aggr, m_par, e_set,
        sr.Copula, sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std)

        par_final = init_draw
        if init_success == false
            error("Couldn't find initial value that produces posterior")
        end
    else 
        par_final = copy(er.par_final)
    end

    draws_raw, posterior, accept_rate = rwmh(par_final, hessian_sym, sr.n_par, er.Data, er.Data_missing,
    er.H_sel, sr.XSSaggr,lr.A, lr.B, lr.LOMstate, lr.State2Control, sr.indexes,sr.indexes_aggr, m_par, e_set,
    sr.Copula, sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std)

    draws = draws_raw[e_set.burnin+1:end, :]

    ##
    parnames_ascii = collect(metaflatten(m_par, label))
    if e_set.me_treatment != :fixed
    for i in eachindex(e_set.meas_error_input)
        push!(parnames_ascii, string("sigma_me_", e_set.meas_error_input[i]))
    end
    end

    # need to ] add MCMCChains#master to get the sorted keyword
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
    sr.Copula, sr.distrSS, sr.compressionIndexes, er.priors, er.meas_error, er.meas_error_std; smoother = true)
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
end # module HANKrank

