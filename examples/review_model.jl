### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils
using Random
Random.seed!(1234)

# ╔═╡ cbe3de4f-b52f-496d-b44f-85df58e74961
begin
    import Pkg
    cd(joinpath(@__DIR__, ".."))
    Pkg.activate(".")

    using PartialWavesFromScratch
    using PartialWavesFromScratch.Masses
    import PartialWavesFromScratch.amplitudes_compass: get_wavelist, get_wavenames, get_wavebasis

    using Plots
    using DelimitedFiles
    using DataFrames
end

# ╔═╡ bb5edd40-565c-47a8-ac32-3403c2cb48d4
randrange(N, (a, b)) = rand(N) .* (b - a) .+ a

# ╔═╡ 04fcd55a-8708-4d0d-b264-7a8142461135
df = let N = 100
    s = randrange(N, 1.5 .+ (0, 1) .* 0.04) .^ 2
    σ1 = randrange(N, (2mπ, 1.5 - mπ)) .^ 2
    cosθ1 = randrange(N, (-1, 1))
    ϕ1 = randrange(N, (-π, π))
    cosθ23 = randrange(N, (-1, 1))
    ϕ23 = randrange(N, (-π, π))

    select(
        DataFrame(; s, σ1, cosθ1, ϕ1, cosθ23, ϕ23),
        :σ1, :cosθ1, :ϕ1, :cosθ23, :ϕ23, :s)
end

# ╔═╡ a5960074-8a9d-434f-be80-58a7cff177b4
begin
    # parameters
    mass_bin_name = "1540_1560" #ARGS[1]#
    tslice = "t1"
end

# ╔═╡ 164997ec-478f-41d7-bece-41d6cac651a8
begin
    wavelist = get_wavelist(joinpath("src", "wavelist_formated.txt");
        path_to_thresholds=joinpath("src", "thresholds_formated.txt"),
        M3pi=Meta.parse(mass_bin_name[1:4]) / 1000)
    # 
    wavenames = get_wavenames(wavelist)
    wavebasis = get_wavebasis(wavelist)
end

# ╔═╡ 2c885f1a-3ef6-4471-85dd-0869e8bdc78a
let
    τ1 = (; df[1, :]...)
    wi = wavebasis[2]
    value = wi(τ1...)
    (; τ1, value)
end

# ╔═╡ a9f5f31c-97b0-4f25-864c-f40c30e747fd
let
    wn, name, J, P, M, ϵ, S, L = wavelist[2, :]
    (; wn, name, J, P, M, ϵ, S, L)
end


# ╔═╡ 1e63c1d0-57fa-4d4f-914e-9ab51b180a0a
let
    wavelist_df = DataFrame(wavelist, [:wn, :name, :J, :P, :M, :ϵ, :S, :L])
    writedlm("waves.txt", eachrow(wavelist_df))

    amplitudes = hcat([
        map(eachrow(df)) do τ
            w(τ...)
        end for w in wavebasis
    ]...)
    writedlm("amplitudes_point1.txt", reim.(amplitudes[1, :]))
end


# ╔═╡ Cell order:
# ╠═cbe3de4f-b52f-496d-b44f-85df58e74961
# ╠═bb5edd40-565c-47a8-ac32-3403c2cb48d4
# ╠═04fcd55a-8708-4d0d-b264-7a8142461135
# ╠═a5960074-8a9d-434f-be80-58a7cff177b4
# ╠═164997ec-478f-41d7-bece-41d6cac651a8
# ╠═2c885f1a-3ef6-4471-85dd-0869e8bdc78a
# ╠═a9f5f31c-97b0-4f25-864c-f40c30e747fd
# ╠═1e63c1d0-57fa-4d4f-914e-9ab51b180a0a
