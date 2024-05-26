# parameters
mass_bin_name = "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"

const M3pi = Meta.parse(mass_bin_name[1:4]) / 1000
const Nwave_extended = 88
#########################################################

using PartialWavesFromScratch.SDMHelper
using PartialWavesFromScratch.PWAHelper
using PartialWavesFromScratch.amplitudes_compass
using PartialWavesFromScratch.PlotHelper

using DelimitedFiles
using Plots
using LinearAlgebra

# MAIN FIT
##### build the model:

# set names
for app in ["rd", "mc", "fu"]
    pwf = path_to_working_folder
    @eval $(Symbol("basisfunc_" * app)) = joinpath($pwf, "variables_$(mass_bin_name)_$(tslice)_" * $app * ".bin")
end

# calculate normfact by getting number of lines
nlines_rd = read(open(basisfunc_rd, "r"), Int32)
nlines_mc = read(open(basisfunc_mc, "r"), Int32)
nlines_fu = read(open(basisfunc_fu, "r"), Int32)

normfact = nlines_fu / nlines_mc * nlines_rd

# read matrix of integrals
BmatFU = read_cmatrix("data/integrmat_$(mass_bin_name)_$(tslice)_fu.txt");
size(readdlm(joinpath(path_wavelist, "wavelist_formated.txt")), 1)
# Model description
wavelist = get_wavelist(joinpath(path_wavelist, "wavelist_formated.txt");
    path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt"),
    M3pi)
# 
const posϵ = [i for (i, ϵ) in enumerate(wavelist[:, 6]) if ϵ == "+"]
const negϵ = [i for (i, ϵ) in enumerate(wavelist[:, 6]) if ϵ == "-"]
const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]
const Nwaves = size(wavelist, 1)


println("Start calculating SDMs")
minpars = readdlm("data/llhfit_1540_1560_t1_87374.txt")
SDM = normfact * pars_to_SDM(minpars, BmatFU, ModelBlocks)

threshold_mask = get_threshold_mask(joinpath(path_wavelist, "thresholds_formated.txt"), M3pi, 88)
SDM_enlarged = enlarge_with_zeros!(SDM, threshold_mask)


# PUBLISHED DATA
SDM_mass_bin_name = string(trunc(Int, 100 - ((2480 - Meta.parse(split(mass_bin_name, "_")[1])) ./ 20)))                #ARGS[2]

#Compass SDM
path_to_SDM = "data/SDMs/0.100000-0.112853/";
SDM_RD = read_compass_SDM(joinpath(path_to_SDM, "sdm$(SDM_mass_bin_name)."),
    path_to_wavelist=joinpath(path_wavelist, "wavelist_formated.txt"),
    path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt");
    M3pi)

SDM_RD_err = read_compass_SDM(joinpath(path_to_SDM, "sdm$(SDM_mass_bin_name)-err."),
    path_to_wavelist=joinpath("src", "wavelist_formated.txt"),
    path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt");
    M3pi)

pars0 = SDM_to_pars(SDM_RD ./ normfact, BmatFU, ModelBlocks)
abs.(diag(pars_to_SDM(pars0, BmatFU, ModelBlocks) .* normfact) .- diag(SDM_RD))

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

# SDM MAIN FIT PLOTTING

let
    bar(real.(diag(SDM_enlarged)), lab="S.S", color="red", size=(800, 500))
    scatter!(real.(diag(SDM_RD)), yerror=real.(diag(SDM_RD_err)), lab="F.H.-D.R.", color="green",
        xlab="# wave", ylab="Magnitude", title="Diagonal of the SDM_$(mass_bin_name)_$(tslice)", size=(800, 500))
end






# BOOTSTRAP PLOTTING
# does not work
bootstrap_file = "data/llhfit_bootstrap_$(mass_bin_name)_$(tslice).txt"
BootstrapResults = readdlm(bootstrap_file)
Nbstrap_attempts = vcat(size(BootstrapResults[1, :])...)[1]
SDMs = [Array{Complex{Float64}}(undef, Nwave_extended, Nwave_extended) for idx in 1:Nbstrap_attempts]
@time for b in 1:Nbstrap_attempts
    _br = BootstrapResults[2:end, b]
    SDMs[b] .= enlarge_with_zeros!(normfact * pars_to_SDM(_br, BmatFU, ModelBlocks), threshold_mask)
end

plotBSTsummary(SDMs, SDM_enlarged, SDM_RD, SDM_RD_err) #  tosort=true, toannotate=true
savefig(joinpath("plots/bootstrap_combined", "bootstrap_combined_$(mass_bin_name)_$(tslice)_test.pdf"))
#saveBSTSplot(joinpath("plots/bootstrap_summary","bootstrap_summary_$(mass_bin_name)_$(tslice)_test.pdf"),
# SDMs, SDM, SDM_RD, SDM_RD_err)
