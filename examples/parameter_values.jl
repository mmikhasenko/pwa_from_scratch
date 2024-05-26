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

# ## parameters
# To get the parameters from the matrix, one needs:
# - SDM
# - matrix of integrals
# - model description as a ModelBlocks



# set names
normfact = let
    for app in ["rd", "mc", "fu"]
        pwf = path_to_working_folder
        @eval $(Symbol("basisfunc_" * app)) = joinpath($pwf, "variables_$(mass_bin_name)_$(tslice)_" * $app * ".bin")
    end
    # calculate normfact by getting number of lines
    nlines_rd = read(open(basisfunc_rd, "r"), Int32)
    nlines_mc = read(open(basisfunc_mc, "r"), Int32)
    nlines_fu = read(open(basisfunc_fu, "r"), Int32)
    nlines_fu / nlines_mc * nlines_rd
end


# Model description
wavelist = get_wavelist(joinpath(path_wavelist, "wavelist_formated.txt");
    path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt"),
    M3pi)
# 
const ModelBlocks = let
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
pars0 = SDM_to_pars(SDM_RD ./ normfact, BmatFU, ModelBlocks)
pars_to_SDM

weights = let
    pars = pars0
    # 
    Tmap = get_parameter_map(ModelBlocks, size(BmatFU, 1))
    pblocks = make_pblock_inds(ModelBlocks)
    pars_bl = fill(0.0, length(pars))
    bl = pblocks[2]
    pars_bl[bl] .= pars[bl]
    shrnk(pars_bl, Tmap)
    # 
end


basis = get_wavebasis(wavelist)
weights



τ1_test = (
    σ1=0.6311001857724697,
    cosθ1=-0.36619233111451877, ϕ1=0.09298675596700612,
    cosθ23=-0.611301179735489, ϕ23=0.6244178754076133, s=2.3253174651821458)
# 
sum(zip(weights, basis)) do (w, b)
    w * b(τ1...)
end