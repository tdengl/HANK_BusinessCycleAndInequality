@doc raw"""
Rank version of SGU_estim()
    SGU_estim(XSS,A,B,m_par,n_par,indexes_aggr,distrSS;estim)

Calculate the linearized solution to the non-linear difference equations defined
by function [`Fsys`](@ref), while only differentiating with respect to the
aggregate part of the model, [`Fsys_agg()`](@ref).

The partials of the Jacobian belonging to the heterogeneous agent part of the model
are taken from the full-model derivatives provided as arguments, `A` and `B` (computed
by [`SGU()`](@ref)).

# Arguments
- `XSS`: steady state around which the system is linearized
- `A`,`B`: derivative of [`Fsys()`](@ref) with respect to arguments `X` [`B`] and
    `XPrime` [`A`]
- `m_par::ModelParameters`, `n_par::NumericalParameters`: `n_par.sol_algo` determines
    the solution algorithm
- `indexes::IndexStruct`,`indexes_aggr::IndexStructAggr`: access aggregate states and controls by name
- `distrSS::Array{Float64,3}`: steady state joint distribution

# Returns
as in [`SGU()`](@ref)
"""
# XSSaggr=sr.XSS
# n_par = sr.n_par 
# indexes = sr.indexes 
# indexes_aggr = sr.indexes_aggr 
# estim = false
function SGU_rank(XSSaggr::Array, m_par::ModelParameters, n_par::NumericalParameters, indexes::IndexStruct,
    indexes_aggr::IndexStructAggr; estim=false)

    ############################################################################
    # Calculate dericatives of non-lineear difference equation
    ############################################################################
    # X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,tuple(zeros(5)...))
    # X0 = zeros(length_X0)
    # F  = Fsys_agg(X0,X0,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC, IDC, DCD, IDCD)
    # if maximum(abs.(F))/10>n_par.ϵ
    #     @warn  "F=0 is not at required precision"
    # indicators = findall(abs.(F).>n_par.ϵ)
    # end

    length_X0   = length(XSSaggr) 
    # # different function here
    # #Fsys_rank(X::AbstractArray, XPrime::AbstractArray, Xss::Array{Float64,1}, m_par::ModelParameters,
    # #n_par::NumericalParameters, indexes::IndexStructAggr)
    # # n_par::NumericalParameters, indexes::Union{IndexStructAggr,IndexStruct})
    # # Fsys_agg(X::AbstractArray, XPrime::AbstractArray, Xss::Array{Float64,1},distrSS::AbstractArray, m_par::ModelParameters,
    # # ::NumericalParameters, indexes::Union{IndexStructAggr,IndexStruct})
    # # # X0 = zeros(length_X0) .+ ForwardDiff.Dual(0.0,tuple(zeros(5)...))
    # X0 = zeros(length_X0)
    # F  = Fsys_agg(X0,X0,XSSaggr, zeros(1,1),m_par,n_par,indexes)
    # if maximum(abs.(F))/10>n_par.ϵ
    #     @warn  "F=0 is not at required precision"
    #     indicators = findall(abs.(F).>n_par.ϵ)
    # end





    BA          = ForwardDiff.jacobian(x-> Fsys_agg(x[1:length_X0],x[length_X0+1:end],XSSaggr, zeros(1,1),m_par,n_par,indexes),zeros(2*length_X0))
    A          = BA[:,length_X0+1:end]
    B          = BA[:,1:length_X0]

    ############################################################################
    # Solve the linearized model: Policy Functions and LOMs
    ############################################################################
    gx, hx, alarm_sgu, nk = SolveDiffEq(A,B, n_par, estim)
    return gx, hx, alarm_sgu, nk, A, B
end