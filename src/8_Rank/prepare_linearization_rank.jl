
@doc raw"""
    prepare_linearization(KSS, VmSS, VkSS, distrSS, n_par, m_par)

Compute a number of equilibrium objects needed for linearization.

# Arguments
- `KSS`: steady-state housing capital stock
- `VmSS`, `VkSS`: marginal value functions
- `distrSS::Array{Float64,3}`: steady-state distribution of idiosyncratic states, computed by [`Ksupply()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`

# Returns
- `XSS::Array{Float64,1}`, `XSSaggr::Array{Float64,1}`: steady state vectors produced by [`@writeXSS()`](@ref)
- `indexes`, `indexes_aggr`: `struct`s for accessing `XSS`,`XSSaggr` by variable names, produced by [`@make_fn()`](@ref),
        [`@make_fnaggr()`](@ref)
- `compressionIndexes::Array{Array{Int,1},1}`: indexes for compressed marginal value functions (``V_m`` and ``V_h``)
- `Copula(x,y,z)`: function that maps marginals `x`,`y`,`z` to approximated joint distribution, produced by
        [`mylinearinterpolate3()`](@ref)
- `n_par::NumericalParameters`,`m_par::ModelParameters`
- `CDF_SS`, `CDF_m`, `CDF_h`, `CDF_y`: cumulative distribution functions (joint and marginals)
- `distrSS::Array{Float64,3}`: steady state distribution of idiosyncratic states, computed by [`Hsupply()`](@ref)
"""
function prepare_linearization_rank(m_par)
    # Calculate other equilibrium quantities
    KSS, incgrossSS, incSS, NSS, rSS, wSS, YSS, ProfitsSS, ISS, RBSS, taxrevSS, av_tax_rateSS, eff_int = incomes_rank(m_par)

    naggrstates             = length(state_names)                   # number of aggregate states
    naggrcontrols           = length(control_names)                 # number of aggregate controls
    naggr                   = length(aggr_names)                    # number of all aggregates
    # Write changed parameter values to n_par
    n_par                   = NumericalParameters(nstates = naggrstates,ncontrols = naggrcontrols, naggrstates = naggrstates, naggrcontrols = naggrcontrols,
                            aggr_names  = aggr_names, naggr = naggr
                             )
    BSS = KSS./10.0 # What pins down steady state bond level in RANK?
    #TD: Taxes change here
    # Calculate taxes and government expenditures
    TSS                 = taxrevSS  + av_tax_rateSS * ((1.0 .- 1.0 ./ m_par.μw).*wSS.*NSS)
    GSS                 = TSS - (m_par.RB./m_par.π-1.0)*BSS


    # aggregate steady state marker
    # @include "../3_Model/input_aggregate_steady_state_rank.jl

    # write to XSS vector
    @writeXSS
    
    # produce indexes to access XSS etc.
    indexes_aggr          = produce_indexes_aggr(n_par)
    indexes               = produce_indexes(n_par)
    ntotal                = indexes.profits # Convention: profits is the last control in the list of control variables
    @set! n_par.ntotal    = ntotal   
    @set! n_par.ncontrols = n_par.naggrcontrols
    @set! n_par.LOMstate_save = zeros(n_par.nstates, n_par.nstates)
    @set! n_par.State2Control_save = zeros(n_par.ncontrols, n_par.nstates)
    
    return XSS, XSSaggr, indexes, indexes_aggr, n_par, m_par       
    end