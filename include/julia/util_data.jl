# - returns indexes for  individuals with at
#   least n_days_min of recordings
function GetIndexSubjects(n_days_min::Int64)

    # 12 obs span 1h
    n_obs_min = 12 * 24 * n_days_min
    # - subjects excluded from this study
    sbj_exclude = vcat([2, 3, 5, 8, 12, 13, 14, 19,
                        37, 39, 42, 43], collect(45:55))
    sbj_temp = setdiff(collect(1:55), sbj_exclude)
    sbj_study = Int64[]

    for sbj in sbj_temp
        println(pwd())
        csv_reader = CSV.File(pwd()*"/data/S$sbj.csv")
        temperatureRaw = convert(Array{Float64}, csv_reader.Temp)
        N = length(temperatureRaw)
        if (N >= n_obs_min)
            append!(sbj_study, sbj)
        end
    end

    return sbj_study
end


# - ?
function ConvertArray(x)
    N = length(x)
    idx_obs = findall(x .!= "NA")
    idx_NA = setdiff(1:N, idx_obs)
    if length(idx_NA) == 0
        return convert(Array{Float64}, x)
    else
        out = zeros(Float64, N)
        out[idx_obs] .= parse.(Float64, x[findall(x .!= "NA")])
        out[idx_NA] .= NaN
        return out
    end
end
