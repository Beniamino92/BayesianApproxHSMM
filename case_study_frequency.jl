using PyPlot; pygui(true)
using CSV
using ProgressMeter
using Dates
using DataFrames
using KernelDensity
using RCall



# - set your current directory
# cd("/Users/beniamino/My_PostDoc/Extended_BayesApproxHSMM/temp/")
cd("/Users/beniamino/Desktop/reviewed_BayesianApproxHSMM/")
# cd("/Users/beniamino/Desktop/Extended_BayesApproxHSMM/temp/")
utils = pwd()*"/include/julia/".*["util.jl",
                                  "util_gibbs.jl",
                                  "util_data.jl"]
[include(path) for path in utils]

# temperature and physical activity
dataset = CSV.File(pwd()*"/data/S16.csv")
dformat = Dates.DateFormat("yyyy-mm-dd HH:MM:SS")
time = Dates.DateTime.(dataset.time, dformat)
activityRaw = ConvertArray(dataset.activity)
N = length(activityRaw)


# pre-processing data
activity = Standardize(NaNMath.sqrt.(activityRaw))

# hyperparms gibbs sampler
ω_mixing = 0.1 # ∈ [0, 1], mixing for draw periodogram or random walk
ω_upper = 0.3 # ∈ [0, 0.5] upper bound frequencies
σ_RW = 1/(25*N) # sd proposal random walk frequencies
σ_β = 5 # prior variance for β,  β ∼ Normal(0, σ_β * I)
ν0 = 4 # prior σ², InverseGamma(ν0/2, η0/2)
γ0 = 1 # prior σ², InverseGamma(ν0/2, η0/2)
hp = HyperPars(σ_RW, σ_β, ν0, γ0, ω_mixing, ω_upper)

N_MCMC = Int64(5e3) # n of gibbs iterations
N_draw = 5 # n draws from posterior predictive
M = 1

# - get spectrum ("lomb" for NAs, if needed)
R"""
library("lomb")
spectrum = lsp($activity, plot = FALSE)
"""
spectrum = rcopy(R"spectrum")
ω_start = spectrum[:peak_at][1]

# MCMC - gibbs sampler
MCMC_act = GibbsSamplerOscillatory(activity, M, N_MCMC, hp; plt = true,
                               ω_start = ω_start)


# acceptance rate
AcceptanceRate(vec(MCMC_act[:sample].ω))

# plot posterior predictive
PlotOscillatory(activity, time, MCMC_act, 10)

# posterior mean
MCMC_act[:ω]
axhline()
