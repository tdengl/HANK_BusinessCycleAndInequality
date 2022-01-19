# Variables not included in rank model:  σ, Ht, BD, "GiniC", "GiniX", "sdlogy", "I90share","I90sharenet",  "w90share"
# Variables changed to state: K
const shock_names = [
    :A, :Z, :ZI, :μ, :μw, :Gshock, :Rshock, :Sshock, :Tlevshock, :Tprogshock
]

const state_names = [
    "A","K","B", "Z", "ZI", "RB" , "μ", "μw", "union_Retained", "Retained", "Ylag",
    "Blag", "Tlag", "Ilag", "wlag", "qlag", "Clag", "av_tax_ratelag", "τproglag",
    "Gshock", "Tprogshock", "Rshock", "Sshock",
]

const control_names = [
    "r", "w", "π" ,"πw" ,"Y","C", "q",  "N", "mc", "mcw", "u","av_tax_rate", "T", "I", "BY","TY", "mcww", "G", "τlev", "τprog", "Ygrowth", "Bgrowth",
    "Igrowth", "wgrowth", "Cgrowth", "Tgrowth", "LP", "LPXA", "totRetainedY",
    "union_firm_profits", "unionprofits", "firm_profits", "profits"
]

# ascii names used for cases where unicode doesn't work, e.g., file saves
const state_names_ascii = [
    "A","K","B", "Z", "ZI", "RB" , "mu","muw", "union_Retained", "Retained", "Ylag",
    "Blag", "Tlag", "Ilag", "wlag", "qlag", "Clag", "av_tax_ratelag", "tauproglag",
    "Gshock", "Tprogshock", "Rshock", "Sshock"
]

const control_names_ascii = [ 
    "r", "w", "pi" , "piw","Y", "C", "q",  "N", "mc", "mcw", "u",
    "av_tax_rate", "T", "I", "BY", "TY", "mcww", "G", "taulev", "tauprog", "Ygrowth",
    "Bgrowth", "Igrowth", "wgrowth", "Cgrowth", "Tgrowth", "LP", "LPXA",
    "totRetainedY", "union_firm_profits", "unionprofits", "firm_profits", "profits"
]
# All names in one array
const aggr_names         = [state_names; control_names]
const aggr_names_ascii   = [state_names_ascii; control_names_ascii]
# Identify distributional summary variables
const distr_names         = ["GiniC", "GiniX", "sdlogy", "I90share","I90sharenet",  "w90share"]
