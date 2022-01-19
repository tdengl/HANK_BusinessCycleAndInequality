function incomes_rank(m_par)

    RBSS      = m_par.RB
    rSS       = RBSS
    kn        = ((rSS - 1 + m_par.δ_0)./(1.0 / m_par.μ*m_par.α)).^(1.0/(m_par.α-1)) # productivity is effectively 1.0/m_par.μ
    KSS       = (kn*(1.0 ./ (m_par.μ * m_par.μw) .* (1.0 - m_par.α) .* (m_par.τ_lev .* 
                (1.0 - m_par.τ_prog)).^(1.0 / (1.0 - m_par.τ_prog))).^((1.0 - m_par.τ_prog)./(m_par.γ + m_par.τ_prog + (m_par.α) 
                .* (1 - m_par.τ_prog)))).^(1.0/(1.0 - m_par.α.*((1.0 - m_par.τ_prog)./(m_par.γ + m_par.τ_prog + (m_par.α) 
                .* (1 - m_par.τ_prog))))) # productivity here is effectively 1.0/1.0 ./ (m_par.μ * m_par.μw)
    NSS       = employment(KSS, 1.0 ./ (m_par.μ * m_par.μw), m_par) # same as NSS = KSS/kn
    wSS       = wage(KSS, 1.0 ./m_par.μ, NSS, m_par)                    # wages
    YSS       = output(KSS, 1.0, NSS, m_par)

    ProfitsSS = (1.0 -1.0 / m_par.μ) .* YSS
    ISS       = m_par.δ_0 * KSS
    RBSS      = m_par.RB
    GHHFA     = ((m_par.γ + m_par.τ_prog) / (m_par.γ + 1.0))              # transformation (scaling) for composite good

    
    eff_int   = RBSS # efffective rate

    mcw = 1.0 ./ m_par.μw

    incgrossSS              = mcw .* wSS .* NSS                                   # gross labor income
    incgrosspSS             = ProfitsSS                                         # gross profit income
    incpSS                  = m_par.τ_lev .* ProfitsSS.^(1.0 - m_par.τ_prog)                  # profit income net of taxes
    incSS                   = m_par.τ_lev .* (mcw .*wSS .* NSS).^(1.0 - m_par.τ_prog)     # net labor income
    taxrevSS                = incgrossSS - incSS + incgrosspSS - incpSS                           # tax revenues

    incgrossaux             = incgrossSS
    av_tax_rateSS           = taxrevSS/incgrossaux


    return KSS, incgrossSS, incSS, NSS, rSS, wSS, YSS, ProfitsSS, ISS, RBSS, taxrevSS, av_tax_rateSS, eff_int
end