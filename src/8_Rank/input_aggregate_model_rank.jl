#------------------------------------------------------------------------------
# THIS FILE CONTAINS THE "AGGREGATE" MODEL EQUATIONS, I.E. EVERYTHING  BUT THE 
# HOUSEHOLD PLANNING PROBLEM. THE lATTER IS DESCRIBED BY ONE EGM BACKWARD STEP AND 
# ONE FORWARD ITERATION OF THE DISTRIBUTION.
#
# AGGREGATE EQUATIONS TAKE THE FORM 
# F[EQUATION NUMBER] = lhs - rhs
#
# EQUATION NUMBERS ARE GENEREATED AUTOMATICALLY AND STORED IN THE INDEX STRUCT
# FOR THIS THE "CORRESPONDING" VARIABLE NEEDS TO BE IN THE LIST OF STATES 
# OR CONTROLS.
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# AUXILIARY VARIABLES ARE DEFINED FIRST
#------------------------------------------------------------------------------


# Elasticities and steepness from target markups for Phillips Curves
η                       = μ / (μ - 1.0)                                 # demand elasticity
κ                       = η * (m_par.κ / m_par.μ) * (m_par.μ - 1.0)     # implied steepness of phillips curve
ηw                      = μw / (μw - 1.0)                               # demand elasticity wages
κw                      = ηw * (m_par.κw / m_par.μw) * (m_par.μw - 1.0) # implied steepness of wage phillips curve

# Capital Utilization
MPK_SS                  = exp(Xss[indexes.rSS]) - 1.0 + m_par.δ_0       # stationary equil. marginal productivity of capital
δ_1                     = MPK_SS                                        # normailzation of utilization to 1 in stationary equilibrium
δ_2                     = δ_1 .* m_par.δ_s                              # express second utilization coefficient in relative terms
# Auxiliary variables
Kserv                   = K * u                                         # Effective capital
MPKserv                 = mc .* Z .* m_par.α .* (Kserv ./ N) .^(m_par.α - 1.0)      # marginal product of Capital
depr                    = m_par.δ_0 + δ_1 * (u - 1.0) + δ_2 / 2.0 * (u - 1.0)^2.0   # depreciation

Wagesum                 = N * w                                         # Total wages in economy t
WagesumPrime            = NPrime * wPrime                               # Total wages in economy t+1

YREACTION = Y ./ exp(Xss[indexes.YSS])                                  # Policy reaction function to Y

# tax progressivity variabels used to calculate e.g. total taxes
#tax_prog_scale          = (m_par.γ + m_par.τ_prog) / ((m_par.γ + τprog))  # scaling of labor disutility including tax progressivity
inc                      = τlev .* (mcw .*w .* N).^(1.0 - m_par.τ_prog)           # net labor income
incp                     = τlev .* profits.^(1.0 - m_par.τ_prog)                  # profit income net of taxes

incgross                 = mcw .* w .* N                                   # gross labor income
incgrossp                = profits                                         # gross profit income

taxrev                  = incgross - inc + incgrossp - incp                           # tax revenues
incgrossaux             = incgross


############################################################################
#           Error term calculations (i.e. model starts here)          #
############################################################################

#-------- States -----------#
# Error Term on exogeneous States
# Shock processes
F[indexes.Gshock]       = log.(GshockPrime)         - m_par.ρ_Gshock * log.(Gshock)     # primary deficit shock

F[indexes.Rshock]       = log.(RshockPrime)         - m_par.ρ_Rshock * log.(Rshock)     # Taylor rule shock
F[indexes.Sshock]       = log.(SshockPrime)         - m_par.ρ_Sshock * log.(Sshock)     # uncertainty shock

# Stochastic states that can be directly moved (no feedback)
F[indexes.A]            = log.(APrime)              - m_par.ρ_A * log.(A)               # (unobserved) Private bond return fed-funds spread (produces goods out of nothing if negative)
F[indexes.Z]            = log.(ZPrime)              - m_par.ρ_Z * log.(Z)               # TFP
F[indexes.ZI]           = log.(ZIPrime)             - m_par.ρ_ZI * log.(ZI)             # Investment-good productivity

F[indexes.μ]            = log.(μPrime./m_par.μ)     - m_par.ρ_μ * log.(μ./m_par.μ)      # Process for markup target
F[indexes.μw]           = log.(μwPrime./m_par.μw)   - m_par.ρ_μw * log.(μw./m_par.μw)   # Process for w-markup target

F[indexes.Ylag]         = log(YlagPrime)    - log(Y)
F[indexes.Blag]         = log(BlagPrime)    - log(B)
F[indexes.Ilag]         = log(IlagPrime)    - log(I)
F[indexes.wlag]         = log(wlagPrime)    - log(w)
F[indexes.Tlag]         = log(TlagPrime)    - log(T)
F[indexes.qlag]         = log(qlagPrime)    - log(q)
F[indexes.Clag]         = log(ClagPrime)    - log(C)
F[indexes.av_tax_ratelag] = log(av_tax_ratelagPrime) - log(av_tax_rate)

# Growth rates
F[indexes.Ygrowth]      = log(Ygrowth)      - log(Y/Ylag)
F[indexes.Tgrowth]      = log(Tgrowth)      - log(T/Tlag)
F[indexes.Bgrowth]      = log(Bgrowth)      - log(B/Blag)
F[indexes.Igrowth]      = log(Igrowth)      - log(I/Ilag)
F[indexes.wgrowth]      = log(wgrowth)      - log(w/wlag)
F[indexes.Cgrowth]      = log(Cgrowth)      - log(C/Clag)

#  Taylor rule and interest rates
F[indexes.RB]           = log(RBPrime) - Xss[indexes.RBSS] -
                            ((1 - m_par.ρ_R) * m_par.θ_π) .* log(π) -
                            ((1 - m_par.ρ_R) * m_par.θ_Y) .* log(YREACTION) -
                            m_par.ρ_R * (log.(RB) - Xss[indexes.RBSS])  - log(Rshock)

F[indexes.τlev]         = av_tax_rate - (taxrev ./ incgrossaux) # Union profits are taxed at average tax rate, this equation pins down taxrev because av_tax rate is exogenous

F[indexes.T]            = log(T) - log(taxrev + av_tax_rate * unionprofits)

F[indexes.av_tax_rate]  = log(av_tax_rate) - m_par.ρ_τ * log(av_tax_ratelag)  - 
                            (1.0 - m_par.ρ_τ) * Xss[indexes.av_tax_rateSS] -
                            (1.0 - m_par.ρ_τ) * m_par.γ_Yτ * log(YREACTION) -
                            (1.0 - m_par.ρ_τ) * m_par.γ_Bτ * (log(B) - Xss[indexes.BSS])
# --------- Controls ------------
# Deficit rule
F[indexes.π]            = log(BgrowthPrime) + m_par.γ_B * (log(B)- Xss[indexes.BSS])  -
                          m_par.γ_Y * log(YREACTION)  - m_par.γ_π * log(π) - log(Gshock)

F[indexes.G]            = log(G) - log(BPrime + T - RB / π * B)             # Government Budget Constraint

# Phillips Curve to determine equilibrium markup, output, factor incomes 
F[indexes.mc]           = (log.(π)- Xss[indexes.πSS]) - κ *(mc - 1 ./ μ ) -
                            m_par.β * ((log.(πPrime) - Xss[indexes.πSS]) .* YPrime ./ Y) 
                        

# Wage Phillips Curve 
F[indexes.mcw]          = (log.(πw)- Xss[indexes.πwSS]) - (κw * (mcw - 1 ./ μw) +
                            m_par.β * ((log.(πwPrime) - Xss[indexes.πwSS]) .* WagesumPrime ./ Wagesum))
# worker's wage = mcw * firm's wage
# Wage Dynamics
F[indexes.πw]           = log.(w ./ wlag) - log.(πw ./ π)                   # Definition of real wage inflation

# Capital utilisation
F[indexes.u]            = MPKserv  -  q * (δ_1 + δ_2 * (u - 1.0))           # Optimality condition for utilization

# Prices
F[indexes.r]            = log.(r) - log.(1 + MPKserv * u - q * depr )       # rate of return on capital

F[indexes.mcww]         = log.(mcww) - log.(mcw * w)                        # wages that workers receive

F[indexes.w]            = log.(w) - log.(wage(Kserv, Z * mc, N, m_par))     # wages that firms pay

F[indexes.unionprofits] = log.(unionprofits)  - log.(w.*N .* (1.0 - mcw))  # profits of the monopolistic unions

F[indexes.profits]      = log.(profits)  - log.(Y .* (1.0 - mc) .+ q .* (KPrime .- (1.0 .- depr) .* K) .- I)  # profits of the monopolistic resellers

F[indexes.q]            = 1.0 - ZI * q * (1.0 - m_par.ϕ / 2.0 * (Igrowth - 1.0)^2.0 - # price of capital investment adjustment costs
                          m_par.ϕ * (Igrowth - 1.0) * Igrowth)  -
                          m_par.β * ZIPrime * qPrime * m_par.ϕ * (IgrowthPrime - 1.0) * (IgrowthPrime)^2.0
# Asset market premia
F[indexes.LP]           = log.(LP)                  - (log((q + r - 1.0)/qlag) - log(RB / π))                   # Ex-post liquidity premium           
F[indexes.LPXA]         = log.(LPXA)                - (log((qPrime + rPrime - 1.0)/q) - log(RBPrime / πPrime))  # ex-ante liquidity premium



# Aggregate Quantities
F[indexes.I]            = KPrime .-  K .* (1.0 .- depr)  .- ZI .* I .* (1.0 .- m_par.ϕ ./ 2.0 .* (Igrowth -1.0).^2.0)           # Capital accumulation equation
F[indexes.N]            = log.(N) - log.(((1.0 - m_par.τ_prog) * τlev * (mcw .* w).^(1.0 - m_par.τ_prog)).^(1.0 / (m_par.γ + m_par.τ_prog)))   # labor supply
F[indexes.Y]            = log.(Y) - log.(Z .* N .^(1.0 .- m_par.α) .* Kserv .^m_par.α)                                          # production function
F[indexes.C]            = log.(Y .- G .- I .+ (A .- 1.0) .* RB .* B ./ π) .- log(C)                            # Resource constraint
# TD: Does not contain  - (δ_1 * (u - 1.0) + δ_2 / 2.0 * (u - 1.0)^2.0).*K) anymore

# Error Term on prices/aggregate summary vars (logarithmic, controls), here difference to SS value averages
F[indexes.BY]           = log.(BY)    - log.(B/Y)                                                               # Bond to Output ratio
F[indexes.TY]           = log.(TY)    - log.(T/Y)                                                               # Tax to output ratio

# Now Euler equations missing

F[indexes.K]            = 1.0 ./C.^m_par.ξ - m_par.β/q* 1.0 ./CPrime.^m_par.ξ*(qPrime + rPrime - 1)
F[indexes.B]            = 1.0 ./C.^m_par.ξ - m_par.β*1.0 ./CPrime.^m_par.ξ*RBPrime./πPrime

