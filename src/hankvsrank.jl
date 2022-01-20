# Run main files here, make sure pwd is set to folder containing mainscript and HANKEstim
include("2_NumericalBasics/PreprocessInputs.jl")

push!(LOAD_PATH, pwd())

using HANKEstim, LinearAlgebra, Plots
# initialize model parameters
m_par = ModelParameters() # log utility for now
# Could load and add estimated parameters here
HANKEstim.@load "7_Saves/parameter_example.jld2" par_final
m_par = HANKEstim.Flatten.reconstruct(m_par, par_final[1:length(par_final)-length(HANKEstim.e_set.meas_error_input)])

sr = compute_steadystate(m_par)
lr = linearize_full_model(sr, m_par)

# Now Rank
include("8_Rank/PreprocessInputs_rank.jl")

push!(LOAD_PATH, pwd())

using HANKrank, LinearAlgebra, Plots

# m_par_rank =  HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92,
# ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, 
# θ_π = 2.1, θ_Y = 0.176, γ_B = 0, γ_π = -1.05, γ_Y = -0.852, ωF = 0.01, ωU = 0.01)
# m_par =  ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92,
# ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, 
# θ_π = 2.1, θ_Y = 0.176, γ_B = 0.0, γ_π = -1.05, γ_Y = -0.852, ωF = 0.01, ωU = 0.01)
# m_par = @set m_par.RB = 1.0./m_par.β # Rank interest rate on bonds is simply 1/beta
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, θ_Y = 0.176)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, 
# θ_Y = 0.176, γ_B = 0.4)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, 
# θ_Y = 0.176, γ_B = 0.4,  γ_π = -1.05)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, 
# θ_Y = 0.176, γ_B = 0.4,  γ_π = -1.05, γ_Y = -0.852)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, 
# θ_Y = 0.176, γ_B = 0.4,  γ_π = -1.05, γ_Y = -0.852, ωF = 0.1)
# m_par_rank = HANKrank.ModelParameters(ξ = 1.0, hs = 0.5, δ_s = 1.92, ϕ = 1.267, κ = 0.083, κw = 0.085, ρ_R = 0.818, σ_Rshock = 0.0025, θ_π = 2.1, 
# θ_Y = 0.176, γ_B = 0.4,  γ_π = -1.05, γ_Y = -0.852, ωF = 0.1, ωU = 0.1)

m_par_rank = HANKrank.ModelParameters()

HANKrank.@load "7_Saves/parameter_example.jld2" par_final
m_par = HANKrank.Flatten.reconstruct(m_par_rank, par_final[1:length(par_final)-length(HANKrank.e_set.meas_error_input)])



m_par_rank = HANKrank.@set m_par_rank.RB = 1.0./m_par_rank.β # Rank interest rate on bonds is simply 1/beta
sr_rank = HANKrank.compute_steadystate(m_par_rank)
lr_rank = HANKrank.linearize_full_model(sr_rank, m_par_rank)


# irfs, hank
# Monetary policy shock, original
irf_horizon         = 16
x0                  = zeros(size(lr.LOMstate,1), 1)
x0[sr.indexes.Rshock]= 100 * m_par.σ_Rshock
MX                  = [I; lr.State2Control]
x                   = x0 * ones(1, irf_horizon + 1)
IRF_state_sparse    = zeros(size(MX)[1], irf_horizon)

for t = 1:irf_horizon
        IRF_state_sparse[:, t] = (MX * x[:, t])'
        x[:, t+1]              = lr.LOMstate * x[:, t]
end

# irfs, hank
# Monetary policy shock
irf_horizon         = 16
x0                  = zeros(size(lr_rank.LOMstate,1), 1)
x0[sr_rank.indexes.Rshock]= 100 * 0.0025
MX                  = [I; lr_rank.State2Control]
x                   = x0 * ones(1, irf_horizon + 1)
IRF_state_sparse_rank    = zeros(size(MX)[1], irf_horizon)

for t = 1:irf_horizon
        IRF_state_sparse_rank[:, t] = (MX * x[:, t])'
        x[:, t+1]              = lr_rank.LOMstate * x[:, t]
end

labelsize= 5
ticksize = 5
legendsize = 5
titlesize = 10
subtitlesize = 5

plt1 = plot(IRF_state_sparse[sr.indexes.RB,2:end], reuse = false, label = "hank", title = "Nominal interest rate", xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize, titlefontsize = subtitlesize)
plt1 = plot!(IRF_state_sparse_rank[sr_rank.indexes.RB,2:end], label = "rank", xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize)

plt2 = plot(IRF_state_sparse[sr.indexes.Y,1:end-1], reuse = false, title = "Output", xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize, titlefontsize = subtitlesize)
plt2 = plot!(IRF_state_sparse_rank[sr_rank.indexes.Y,1:end-1], legend = false, xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize)

plt3 = plot(IRF_state_sparse[sr.indexes.C,1:end-1], reuse = false, title = "Consumption", xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize, titlefontsize = subtitlesize)
plt3 = plot!(IRF_state_sparse_rank[sr_rank.indexes.C,1:end-1], legend = false, xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize)

plt4 = plot(IRF_state_sparse[sr.indexes.I,1:end-1], reuse = false, title = "Investment", xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize, titlefontsize = subtitlesize)
plt4 = plot!(IRF_state_sparse_rank[sr_rank.indexes.I,1:end-1], legend = false, xlabel = "quarters", legendfontsize=legendsize,xguidefontsize = labelsize, xtickfontsize = ticksize)
irf = plot(plt1, plt2, plt3, plt4,layout = (2,2))
labelsize  = 5
ticksize   = 5
legendsize = 5
titlesize  = 10
subtitlesize = 5
ylinewidth = 1
ylinecolor = :black
ylinestyle = :dash

rows = 3
cols = 4
vars = ["Rshock","RB", "Y", "C", "B", "G", "I", "π", "mc", "w", "N", "K"]
irfplot = plot(layout=(rows,cols))
if rows*cols != length(vars)
    error("subplots not equal to number of vars")
end
for i = 1:length(vars)
    varnamesymb = Symbol(vars[i])
    index       = eval(:(sr.indexes.$varnamesymb))
    index_rank  = eval(:(sr_rank.indexes.$varnamesymb))
    irf         = IRF_state_sparse[index, :]
    irf_rank    = IRF_state_sparse_rank[index_rank, :]
    ylims       = (0.95*minimum(irf), 1.05*maximum(irf))
    irfplot     = plot!(subplot = i, irf, legend = false,label = "hank", titlefontsize = subtitlesize, 
                  title =vars[i], reuse = false, lw = 2, )
    irfplot     = plot!(subplot = i, irf_rank, legend = false,label ="rank", titlefontsize = subtitlesize, 
                  title =vars[i], reuse = false, lw = 2, )
    irfplot     = plot!(subplot = i, [1; irf_horizon], [0; 0], lw=ylinewidth, ls = ylinestyle, lc=ylinecolor, label = "")
end
irfplot