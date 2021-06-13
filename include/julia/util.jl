using Random
using StatsBase
using Distributions
using DSP
using LinearAlgebra
using NaNMath
using RCall


# - * true generating model parameters
struct ModelPars
    ω # frequencies
    β # Fourier coefficients
    σ # innovations
end

# - * hyperprior and hyperparameters
mutable struct HyperPars
    σ_RW # sd proposal random walk frequencies
    σ_β # prior variance for β,  β ∼ Normal(0, σ_β * I)
    ν0 # prior σ², InverseGamma(ν0/2, η0/2)
    γ0 # prior σ², InverseGamma(ν0/2, η0/2)
    ω_mixing # ∈ [0, 1], mixing for draw periodogram or random walk
    ω_upper # ∈ [0, 0.5] upper bound frequencies
end


# - fourier design matrix (with intercept)
function FourierBasis(ω, N::Int64)

    M = length(ω)
    X = ones(Float64, N)
    for j in 1:M
      X = hcat(X, cos.(2π*(1:N)*ω[j]),
                  sin.(2π*(1:N)*ω[j]))
    end
    return X[:, 2:end]
end

# - simulate oscillatory time series
function GenerateData(pars::ModelPars, N::Int64)

    ω = pars.ω
    β = pars.β
    σ = pars.σ
    M = length(ω)

    X = FourierBasis(ω, N::Int64)
    signal = X * β
    noise = rand(Normal(0, σ), N)
    data = signal + noise

    return Dict(:data => signal + noise,
                :signal => signal)
end

# - center by subtracting linear trend from data
function RemoveLinearTrend(data)

    N = length(data)
    X = hcat(ones(N), 1:N)
    β = X \ data
    trend = X * β
    return data - trend
end

# - standardize time series (accept NaN's)
function Standardize(data)
    return (data .- NaNMath.mean(data))./
            NaNMath.std(data)
end
