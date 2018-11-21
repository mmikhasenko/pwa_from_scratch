# parameters
mass_bin_name = ARGS[1]# "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
list_of_files = ARGS[2:end] #Can i do this ?

#########################################################
@show ARGS
push!(LOAD_PATH,"src")
print(typeof(mass_bin_name))
using SDMHelper
using PWAHelper
using amplitudes_compass
using DelimitedFiles

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
SDMs = []
for path_and_filename in list_of_files
    minpars = readdlm(path_and_filename);
    # calculate
    SDM = normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks)
    for i in 1:size(SDM,1)
    	println(SDM[i,i])
    end
    # push!(SDMs,normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks))
    # as you see you also need to load/calculate size(PsiRD) and minpars.
end
