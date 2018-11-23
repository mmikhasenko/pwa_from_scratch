# parameters
mass_bin_name = ARGS[1] # "1280_1300"# "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
list_of_files = ARGS[2:end]

println("Starting with ARGS = ", ARGS)
######################################################
using DelimitedFiles
push!(LOAD_PATH,"src")
using SDMHelper
using FittingPWALikelihood
using amplitudes_compass
using PWAHelper

BmatMC = read_cmatrix(
    joinpath(path_to_working_folder,"integrmat_$(mass_bin_name)_$(tslice)_mc.txt"));

# read wavelist throw waves below threshol
wavelist = get_wavelist(joinpath(path_wavelist,"wavelist_formated.txt");
     path_to_thresholds=joinpath(path_wavelist,"thresholds_formated.txt"),
     M3pi=Meta.parse(mass_bin_name[1:4])/1000)

const noϵ =  [1] # flat wave
const posϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="+"]
const negϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="-"]

# Model description
const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]

# load precalculated data array
const PsiRD = read_precalc_basis(
    joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_rd.txt"));

# get functions to calculate LLH, derivative and hessian
LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiRD, BmatMC, ModelBlocks);

println("Start calculating LLHs")
llhs = [LLH(vcat(readdlm(path_and_filename)...)) for path_and_filename in list_of_files]
sorted_list = sort(collect(zip(list_of_files, llhs)), by=x->x[2])
writedlm(joinpath(path_to_working_folder,"llh_attmpts_$(mass_bin_name)_$(tslice).txt"), sorted_list)
println("Done")
