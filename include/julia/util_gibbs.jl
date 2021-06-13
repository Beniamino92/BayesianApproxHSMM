# - MCMC sample of ω, β and σ
struct PosteriorSample
    ω::Matrix{Float64}
    β::Matrix{Float64}
    σ::Array{Float64}
end

# - acceptance rate MCMC sample
function AcceptanceRate(x)
    (sum(diff(x) .!= 0) + 1)/length(x)
end

# - time series draw from posterior predictive sample?
function DrawFromPredictive(MCMC::PosteriorSample, N::Int64)
    N_sample = size(MCMC.ω, 2)
    ii = sample(1:N_sample)
    signal = FourierBasis(MCMC.ω[:, ii], N) *
             MCMC.β[:, ii]
    f = [rand(Normal(signal[n], MCMC.σ[ii]), 1)[1] for n in 1:N]
    return f
end

# - log likelihood data
function LogLikelihood(data, ω, β, σ)
    N = length(data)
    out = -0.5*N*log(2π) - N*log(σ) -
          NaNMath.sum((data .- FourierBasis(ω, N) * β).^2)
    return out
end

# - log posterior distr of ω
function LogPosteriorFrequencies(data, ω, β, σ)
    N = length(data)
    X = FourierBasis(ω, N)
    out = -NaNMath.sum((data - X*β).^2)/(2*(σ^2))
    return out
end

# - draw ω from the periodogram
function SamplePeriodogram(data, spectrum, ω_current, β_current, σ_current)

    M = length(ω_current)
    freq = spectrum[:scanned][2:(end-1)]
    p = spectrum[:power][2:(end-1)]
    p_norm = p ./sum(p)
    ω_current_aux = copy(ω_current)

    for j in 1:M
      ω_curr = copy(ω_current_aux)
      aux_temp = false # (avoiding design matrix linear dependent)
      while (aux_temp == false)
        # proposing frequencies
        global ω_star = sample(freq, Weights(p_norm))
        global ω_prop = copy(ω_curr)
        ω_prop[j] = ω_star
        if (! (any(vcat(ω_prop[1:(j-1)], ω_prop[(j+1):end]) .== ω_star)))
          aux_temp = true
        end
      end

      # acceptance probabilities
      llk_ratio = LogPosteriorFrequencies(data, ω_prop, β_current, σ_current) -
                  LogPosteriorFrequencies(data, ω_curr, β_current, σ_current)
      lprop_ratio = log(p_norm[searchsortedlast(freq, ω_curr[j])]) -
                 log(p_norm[searchsortedlast(freq, ω_star)])
      MH_ratio = exp(llk_ratio + lprop_ratio)

      U = rand()
      if (U <= min(1, MH_ratio))
        ω_current_aux = ω_prop
      else
        ω_current_aux = ω_curr
      end
    end

    ω_out = sort(ω_current_aux)
    return ω_out
end

# - draw ω using random walk
function RandomWalkFrequencies(data, ω_current, β_current, σ_current,
                               hp::HyperPars)

    M = length(ω_current)
    ω_current_aux = copy(ω_current)

    for j in 1:M

        ω_curr = copy(ω_current_aux)
        aux_temp = false
        while (aux_temp == false)
            global ω_star = rand(Normal(ω_current[j], hp.σ_RW), 1)[1]
            if !(ω_star <= 0 || ω_star >= hp.ω_upper)
              aux_temp = true
            end
        end

        global ω_prop = copy(ω_curr)
        ω_prop[j] = ω_star

        llk_ratio = LogPosteriorFrequencies(data, ω_prop, β_current, σ_current) -
                    LogPosteriorFrequencies(data, ω_curr, β_current, σ_current)

        MH_ratio = exp.(llk_ratio)[1]

        U = rand()
        if (U <= min(1, MH_ratio))
          ω_current_aux = ω_prop
        else
          ω_current_aux = ω_curr
        end
    end

    ω_out = sort(ω_current_aux)
    return ω_out
end

# -  gibbs iteration full posterior
function GibbsStep(data, spectrum, ω_current, β_current, σ_current,
                   hp::HyperPars)

    N = length(data)
    M = length(ω_current)

    # -- * sampling frequencies * --
    U = rand()
    if (U <= hp.ω_mixing)
         # draw from periodogram
        ω_out = SamplePeriodogram(data, spectrum, ω_current,
                                  β_current, σ_current)
    else  # random walk metropolis
        ω_out = RandomWalkFrequencies(data, ω_current, β_current,
                                      σ_current, hp)
    end

    # -- * sampling linear coefficients β * --
    idx_sub = findall(!isnan, data)
    X_post = FourierBasis(ω_out, N)
    X_post = X_post[idx_sub, :]
    β_var_post = inv(I(2*M)/(hp.σ_β^2) + (X_post'*X_post)/(σ_current^2))
    β_var_post = (β_var_post' + β_var_post)/2
    β_mean_post = β_var_post*((X_post'*data[idx_sub])/(σ_current^2))
    β_out = rand(MultivariateNormal(β_mean_post, β_var_post), 1)

    # --  * sampling residual variance σ  * --
    res_var = sum((data[idx_sub] - X_post*β_out).^2)
    ν_post = (N + hp.ν0)/2
    γ_post = (hp.γ0 + res_var)/2
    σ_out  = sqrt.(rand(InverseGamma(ν_post, γ_post), 1))[1]


    output = Dict(:β => β_out, :ω => ω_out, :σ => σ_out)
end

# -  gibbs sampler + bayes estimates + draws posterior predictive
function GibbsSamplerOscillatory(data, M, N_MCMC, hp; plt = false, N_draw = 5,
                                 ω_start = nothing)

    R"""
    library("lomb")
    spectrum = lsp($data, plot = FALSE)
    """
    spectrum = rcopy(R"spectrum")

    N = length(data)

    # objects MCMC

    ω_sample = Matrix{AbstractFloat}(undef, M, N_MCMC+1)
    β_sample = Matrix{AbstractFloat}(undef, 2*M, N_MCMC+1)
    σ_sample = Array{AbstractFloat}(undef, N_MCMC+1)

    # starting values
    if ω_start == nothing
        ω_sample[1:M, 1] = rand(truncated(Uniform(), 1e-5, hp.ω_upper), M)
    else
        ω_sample[1:M, 1] .= ω_start
    end
    β_sample[1:(2*M), 1] = zeros((2*M))
    σ_sample[1] = 10

    @showprogress for tt in 2:(N_MCMC+1)

        ω = ω_sample[1:M, tt-1]
        β = β_sample[1:(2*M), tt-1]
        σ = σ_sample[tt-1]

        MCMC = GibbsStep(data, spectrum, ω, β, σ, hp)
        ω_sample[1:M, tt] = MCMC[:ω]
        β_sample[1:(2*M), tt] = MCMC[:β]
        σ_sample[tt] = MCMC[:σ]
    end

    # bayes estimates
    burn_MCMC = Int64(0.4*N_MCMC)
    posterior_sample = PosteriorSample(ω_sample[:, burn_MCMC:end],
                                       β_sample[:, burn_MCMC:end],
                                       σ_sample[burn_MCMC:end])
    ω̂ = [mean(ω_sample[j, burn_MCMC:end])  for j in 1:M]
    β̂ = [mean(β_sample[j, burn_MCMC:end]) for j in 1:(2*M)]
    σ̂ = mean(σ_sample[burn_MCMC:end])
    llk = LogLikelihood(data, ω̂, β̂, σ̂)
    power = [sum(β̂[(2*j-1):(2*j)].^2) for j = 1:M]
    rate_ω = [AcceptanceRate(ω_sample[j, burn_MCMC:end]) for j in 1:M]

    # plots
    if plt

        # trace plot frequencies
        close()
        for j in 1:M
            subplot(M, 1, j)
            plot(ω_sample[j, burn_MCMC:end])
            axhline(ω̂, color = "red")
        end

        # posterior predictive plot
        println("press enter : posterior predictive check")
        idx_draw = sample(burn_MCMC:(N_MCMC - burn_MCMC + 2),
                          N_draw, replace = true)
        readline()
        close()
        scatter(1:N, data, s = 8, alpha = 0.7);
        plot(1:N, data, alpha = 0.3, color = "black")
        for ii in idx_draw
            plot(1:N, DrawFromPredictive(posterior_sample, N),
                      color = "grey", alpha = 0.1)
        end
        plot(1:N, FourierBasis(ω̂, N) * β̂,
             alpha = 0.8, color = "red", )
    end

    return Dict(:ω => ω̂, :β => β̂, :σ => σ̂,
                :power => power, :rate_ω => rate_ω,
                :sample => posterior_sample,
                :signal => FourierBasis(ω̂, N) * β̂,
                :llk => llk)
end


function PlotOscillatory(activity, time, MCMC_act, N_draw = 5)

    plt = PyPlot.plt
    N = length(activity)
    start = hour(time[1])

    if (start>0) && (start<=12)
        tick_start = findall(hour.(time) .== 12)[1]
    else
        tick_start = findall(hour.(time) .== 0)[1]
    end
    ticks = collect(tick_start:144:N)
    labels = string.(hour.(time)[ticks])
    labels[labels .== "12"] .= "12:00"
    labels[labels .== "0"] .= "00:00"
    start_night = findall(hour.(time) .== 20)[1]

    plt.scatter(1:N, activity, alpha=0.6, s=5, color = "black")
    xlim(-5, N+5);
    ylim((NaNMath.minimum(activity) - 1.5,
          NaNMath.maximum(activity) + 1.5))
    xlabel(""); ylabel("Sqrt(Physical Activity)", fontsize =  13)
    # suptitle("Subject $sbj ", fontsize = 18, fontweight="bold")
    xticks(ticks, labels)
    cA = plt.gca()
    start_night = findall(hour.(time) .== 20)[1]
    for i in 1:8
        cA.add_patch(matplotlib.patches.Rectangle((start_night,
                     NaNMath.minimum(activity)-1.4),144, 0.1,
                     color = "black"))
        start_night += 288
    end
    plt.show()
    [plot(1:N, DrawFromPredictive(MCMC_act[:sample], N),
          color = "grey", alpha = 0.05) for i in 1:N_draw]
    plot(1:N, MCMC_act[:signal], alpha = 0.65, color = "red")
    pow = round(MCMC_act[:power][1], digits = 2)
    p = round(1/((MCMC_act[:ω][1]) *12), digits = 2)
end
