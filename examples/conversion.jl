# parameters
path_to_working_folder = "data"
mass_bin_name = ARGS[1]
tslice = "t1"
tmin = "0.1"
tmax = "0.112853"

######################################################
@show ARGS
# set names
for app in ["rd", "mc", "fu"]
    pwf = path_to_working_folder

    @eval $(Symbol("kinvar_"*app))    = joinpath($pwf,"variables_$(mass_bin_name)_$(tslice)_"*$app*".txt")
    @eval $(Symbol("basisfunc_"*app)) = joinpath($pwf,"functions_$(mass_bin_name)_$(tslice)_"*$app*".txt")
end
# Check that the files exist
#
name_of_the_file_rd = joinpath(path_to_working_folder,"$(mass_bin_name)_rd.root")
isfile(name_of_the_file_rd) || error("Cannot find the file $name_of_the_file_rd")
#
name_of_the_file_mc = joinpath(path_to_working_folder,"$(mass_bin_name)_mc.root")
isfile(name_of_the_file_mc) || error("Cannot find the file $name_of_the_file_mc")

# Converting data.root to data.txt
let
    @time run(pipeline(`src/Extract data/$(mass_bin_name)_rd.root D $(tmin) $(tmax)`, stdout=kinvar_rd))
   @time run(pipeline(`src/Extract data/$(mass_bin_name)_mc.root M $(tmin) $(tmax)`, stdout=kinvar_mc))
  @time run(pipeline(`src/Extract data/$(mass_bin_name)_mc.root F $(tmin) $(tmax)`, stdout=kinvar_fu))
end
