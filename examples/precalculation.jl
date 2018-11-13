# parameters
mass_bin_name = ARGS[1]# "1540_1560"
tslice = "t1"
path_to_working_folder = "data"
path_wavelist = "src"

######################################################
@show ARGS
push!(LOAD_PATH,"src")
using amplitudes_compass
using PWAHelper

# set names
for app in ["rd", "mc", "fu"]
    pwf = path_to_working_folder
    @eval $(Symbol("kinvar_"*app))    = joinpath($pwf,"variables_$(mass_bin_name)_$(tslice)_"*$app*".txt")
    @eval $(Symbol("basisfunc_"*app)) = joinpath($pwf,"functions_$(mass_bin_name)_$(tslice)_"*$app*".txt")
end

# read wavelist throw waves below threshol
wavelist = get_wavelist(joinpath(path_wavelist,"wavelist_formated.txt");
     path_to_thresholds=joinpath(path_wavelist,"thresholds_formated.txt"),
     M3pi=Meta.parse(mass_bin_name[1:4])/1000)

# extract names and wavefunctions from the wavelist
wavenames = get_wavenames(wavelist)
wavebasis = get_wavebasis(wavelist)

# do precalculations
@time precalculate_compass_basis(wavebasis, kinvar_rd, basisfunc_rd)
 @time precalculate_compass_basis(wavebasis, kinvar_mc, basisfunc_mc)
 @time precalculate_compass_basis(wavebasis, kinvar_fu, basisfunc_fu)
