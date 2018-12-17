# parameters
path_to_working_folder = "data"
mass_bin_name = "0520_0540"#ARGS[1]
tslice = "t1"
tmin = "0.1"
tmax = "0.112853"
t_list = ["0.100000-0.112853","0.112853-0.127471","0.127471-0.144385","0.144385-0.164401","0.164401-0.188816"
                ,"0.188816-0.219907","0.219907-0.262177","0.262177-0.326380","0.326380-0.448588","0.448588-0.724294","0.724294-1.000000"]
######################################################
@show ARGS
# set names
for app in ["rd", "mc", "fu"]
    pwf = path_to_working_folder

    @eval $(Symbol("kinvar_"*app))    = joinpath($pwf,"variables_$(mass_bin_name)_$(tslice)_"*$app*".bin")
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
    @time run(`src/Extract data/$(mass_bin_name)_rd.root D $(tmin) $(tmax) $(kinvar_rd)`)
   @time run(`src/Extract data/$(mass_bin_name)_mc.root M $(tmin) $(tmax) $(kinvar_mc)`)
  @time run(`src/Extract data/$(mass_bin_name)_mc.root F $(tmin) $(tmax) $(kinvar_fu)`)
end
