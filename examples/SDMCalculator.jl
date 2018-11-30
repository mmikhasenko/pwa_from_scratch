# parameters
mass_bin_name = ARGS[1] # "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
const M3pi = Meta.parse(mass_bin_name[1:4])/1000
#########################################################


@show ARGS
push!(LOAD_PATH,"src")
using SDMHelper
using PWAHelper
using amplitudes_compass
using DelimitedFiles
using Plots
using LinearAlgebra
using PlotHelper

# MAIN FIT
##### build the model:

# set names
for app in ["rd", "mc", "fu"]
    pwf = path_to_working_folder
    @eval $(Symbol("basisfunc_"*app)) = joinpath($pwf,"variables_$(mass_bin_name)_$(tslice)_"*$app*".txt")
end
# calculate normfact by getting number of lines
nlines_rd = Meta.parse(split(read(`wc -l $(basisfunc_rd)`,String)," ")[1])
nlines_mc = Meta.parse(split(read(`wc -l $(basisfunc_mc)`,String)," ")[1])
nlines_fu = Meta.parse(split(read(`wc -l $(basisfunc_fu)`,String)," ")[1])
normfact = nlines_fu/nlines_mc*nlines_rd
# read matrix of integrals
BmatFU = read_cmatrix("data/integrmat_$(mass_bin_name)_$(tslice)_fu.txt");

# Model description
wavelist = get_wavelist(joinpath(path_wavelist,"wavelist_formated.txt");
         path_to_thresholds=joinpath(path_wavelist,"thresholds_formated.txt"),
         M3pi=M3pi)
const noϵ =  [1] # flat wave
const posϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="+"]
const negϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="-"]
const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]
const Nwaves = size(wavelist,1)


println("Start calculating SDMs")
SDM = []
list_of_files = split(readline("data/llhfit/$(mass_bin_name)/llh_attmpts_$(mass_bin_name)_$(tslice).txt"),"\t")[1]#ARGS[3:end]
path_and_filename = list_of_files
minpars = readdlm(path_and_filename);
SDM = normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks)

threshold_mask = get_threshold_mask(joinpath(path_wavelist,"thresholds_formated.txt"), M3pi, 88)
enlarge_with_zeros!(SDM,threshold_mask)[2,:]


# PUBLISHED DATA

SDM_mass_bin_name =string(trunc(Int,100- ((2480 - Meta.parse(split(mass_bin_name,"_")[1])) ./ 20)))                #ARGS[2]
#Compass SDM
path_to_SDM = "/mnt/data/compass/2008/pwa_results/SDMs/0.100000-0.112853/";
SDM_RD = read_compass_SDM(joinpath(path_to_SDM,"sdm$(SDM_mass_bin_name)."),
    path_to_wavelist   = joinpath("src", "wavelist_formated.txt");
    M3pi = Meta.parse(mass_bin_name[1:4])/1000)

SDM_RD_err = read_compass_SDM(joinpath(path_to_SDM,"sdm$(SDM_mass_bin_name)-err."),
    path_to_wavelist   = joinpath("src", "wavelist_formated.txt");
    path_to_thresholds = joinpath("src","thresholds_formated.txt"),
    M3pi = Meta.parse(mass_bin_name[1:4])/1000)


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

# SDM MAIN FIT PLOTTING

bar(real.(diag(SDM)), lab="S.S",color = "red" ,size=(800,500))
scatter!(real.(diag(SDM_RD)),yerror=real.(diag(SDM_RD_err)), lab="F.H.-D.R.",color = "green",
                        xlab="# wave",ylab ="Magnitude", title = "Diagonal of the SDM_$(mass_bin_name)_$(tslice)", size=(800,500))
savefig(joinpath("plots","sdm_results_$(mass_bin_name)_$(tslice)_test.png"))

# BOOTSTRAP PLOTTING

bootstrap_file = "data/llhfit_bootstrap_$(mass_bin_name)_$(tslice).txt"
BootstrapResults = readdlm(bootstrap_file)
Nbstrap_attempts =vcat(size(BootstrapResults[1,:])...)[1]
SDMs = [Array{Complex{Float64}}(undef, Nwaves,Nwaves) for idx in 1:Nbstrap_attempts]
@time for b in 1:Nbstrap_attempts
    _br = BootstrapResults[2:end,b]
    SDMs[b] .= normfact*pars_to_SDM(_br, BmatFU, ModelBlocks)
end

plotBSTsummary(SDMs, SDM, SDM_RD, SDM_RD_err) #  tosort=true, toannotate=true
savefig(joinpath("plots","bootstrap_combined_$(mass_bin_name)_$(tslice)_test.pdf"))
# saveBSTSplot(joinpath("plots","bootstrap_summary_$(mass_bin_name)_$(tslice)_test.pdf"),
#             SDMs, SDM, SDM_RD, SDM_RD_err)
