# parameters
mass_bin_name =  ARGS[1] #"1540_1560"
SDM_mass_bin_name =ARGS[2]
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
list_of_files = ARGS[3:end] #Can i do this ?

#########################################################
@show ARGS
push!(LOAD_PATH,"src")

using SDMHelper
using PWAHelper
using amplitudes_compass
using DelimitedFiles
using Plots
using LinearAlgebra

#pyplot()
# set names
for app in ["rd", "mc", "fu"]
    pwf = path_to_working_folder
    @eval $(Symbol("basisfunc_"*app)) = joinpath($pwf,"functions_$(mass_bin_name)_$(tslice)_"*$app*".txt")
end
# get number of lines
nlines_rd = Meta.parse(split(read(`wc -l $(basisfunc_rd)`,String)," ")[1])
nlines_mc = Meta.parse(split(read(`wc -l $(basisfunc_mc)`,String)," ")[1])
nlines_fu = Meta.parse(split(read(`wc -l $(basisfunc_fu)`,String)," ")[1])
normfact = nlines_fu/nlines_mc*nlines_rd

# read matrix of integrals
BmatFU = read_cmatrix("data/integrmat_$(mass_bin_name)_$(tslice)_fu.txt");

# read wavelist throw waves below threshol
wavelist = get_wavelist(joinpath(path_wavelist,"wavelist_formated.txt");
         path_to_thresholds=joinpath(path_wavelist,"thresholds_formated.txt"),
         M3pi=Meta.parse(mass_bin_name[1:4])/1000)

const Nwaves = size(wavelist,1)

# Blocks in wavelist
const noϵ =  [1] # flat wave
const posϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="+"]
const negϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="-"]

# Model description
const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]

println("Start calculating SDMs")
#Loop to calculate SDMs
SDM = []
for path_and_filename in list_of_files
    minpars = readdlm(path_and_filename);
    # calculate
    SDM = normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks)
    #for i in 1:size(SDM,1)
    	#println(SDM[i,i])
    #end
    # push!(SDMs,normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks))
    # as you see you also need to load/calculate size(PsiRD) and minpars.
end

#Compass SDM
path_to_SDM = "/mnt/data/compass/2008/pwa_results/SDMs/0.100000-0.112853/";
SDM_RD = read_compass_SDM(joinpath(path_to_SDM,"sdm$(SDM_mass_bin_name)."),
    path_to_wavelist   = joinpath("src", "wavelist_formated.txt");
    path_to_thresholds = joinpath("src","thresholds_formated.txt"),
    M3pi = Meta.parse(mass_bin_name[1:4])/1000)


SDM_RD_err = read_compass_SDM(joinpath(path_to_SDM,"sdm$(SDM_mass_bin_name)-err."),
    path_to_wavelist   = joinpath("src", "wavelist_formated.txt");
    path_to_thresholds = joinpath("src","thresholds_formated.txt"),
    M3pi = Meta.parse(mass_bin_name[1:4])/1000)

#Plot (Assuming that there is only one file, can be put in the loop for SDM)
#bar(hcat(real.(diag(SDM_RD)), real.(diag(SDM))),yerror=hcat(real.(diag(SDM_RD_err)),real.(diag(zeros(size(SDM_RD_err))))), lab=["F.H.-D.R." "S.S."],
        #xlab="# wave", title = "Diagonal of the SDM_$(mass_bin_name)_$(tslice)", size=(800,500))
bar(real.(diag(SDM)), lab="S.S", size=(800,500))
scatter!(real.(diag(SDM_RD)),yerror=real.(diag(SDM_RD_err)), lab="F.H.-D.R.",
                        xlab="# wave",ylab ="Magnitude", title = "Diagonal of the SDM_$(mass_bin_name)_$(tslice)", size=(800,500))

savefig(joinpath("plots","sdm_results_$(mass_bin_name)_$(tslice)_test.png"))
#savefig(joinpath("plots","sdm_results$(mass_bin_name)_$(tslice)_test.pdf"))
