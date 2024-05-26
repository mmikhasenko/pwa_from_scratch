# parameters
mass_bin_name = "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
path_to_SDM = "data/SDMs/0.100000-0.112853/";


const M3pi = Meta.parse(mass_bin_name[1:4]) / 1000
#########################################################

using PartialWavesFromScratch.SDMHelper
using PartialWavesFromScratch.PWAHelper
using PartialWavesFromScratch.amplitudes_compass
using PartialWavesFromScratch.PlotHelper

using DelimitedFiles
using Plots
using LinearAlgebra
using JSON

# ## parameters
# To get the parameters from the matrix, one needs:
# - SDM
# - matrix of integrals
# - model description as a ModelBlocks

# Model description
wavelist = get_wavelist(joinpath(path_wavelist, "wavelist_formated.txt");
    path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt"),
    M3pi)
# 
const ModelBlocks = let
    noϵ = [1] # flat wave
    posϵ = [i for (i, ϵ) in enumerate(wavelist[:, 6]) if ϵ == "+"]
    negϵ = [i for (i, ϵ) in enumerate(wavelist[:, 6]) if ϵ == "-"]
    [noϵ, posϵ, negϵ, negϵ]
end

# Published data
SDM_RD = let
    index = 100 - div(2480 - Meta.parse(split(mass_bin_name, "_")[1]), 20)
    SDM_mass_bin_name = string(index)
    read_compass_SDM(joinpath(path_to_SDM, "sdm$(SDM_mass_bin_name)."),
        path_to_wavelist=joinpath(path_wavelist, "wavelist_formated.txt"),
        path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt");
        M3pi)
end

# read matrix of integrals
BmatFU = read_cmatrix("data/integrmat_$(mass_bin_name)_$(tslice)_fu.txt");
# set names
normfact = let
    nlines_rd, nlines_mc, nlines_fu =
        map(["rd", "mc", "fu"]) do ext
            filename =
                joinpath(path_to_working_folder,
                    "variables_$(mass_bin_name)_$(tslice)_") *
                ext * ".bin"
            # calculate normfact by getting number of lines
            nlines = read(open(filename, "r"), Int32)
        end
    nlines_fu / nlines_mc * nlines_rd
end
BmatFU_norm = BmatFU .* normfact
pars0 = SDM_to_pars(SDM_RD, BmatFU_norm, ModelBlocks)

weights = let
    pars = pars0
    # 
    Tmap = get_parameter_map(ModelBlocks, size(BmatFU, 1))
    pblocks = make_pblock_inds(ModelBlocks)
    pars_bl = fill(0.0, length(pars))
    bl = pblocks[2]
    pars_bl[bl] .= pars[bl]
    shrnk(pars_bl, Tmap)
end

basis = get_wavebasis(wavelist)

# reference point

τ1_test = (
    σ1=0.6311001857724697,
    cosθ1=-0.36619233111451877, ϕ1=0.09298675596700612,
    cosθ23=-0.611301179735489, ϕ23=0.6244178754076133, s=2.3253174651821458)
# 
references = map(b -> b(τ1_test...), basis)

# 
waves_summary = DataFrame(wavelist, [:wn, :name, :J, :P, :M, :ϵ, :S, :L])
waves_summary.weights = weights
waves_summary.references = references

let
    folder = joinpath(@__DIR__, "..", "tests", "references")
    filename = joinpath(folder, mass_bin_name * "_ref.json")
    _waves_summary = transform(waves_summary,
        :references => ByRow(x -> string(x)) => :references,
        :weights => ByRow(x -> string(x)) => :weights
    )
    open(filename, "w") do io
        JSON.print(io, _waves_summary, 4)
    end
end
