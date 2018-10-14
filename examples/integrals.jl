# parameters
mass_bin_name = "2320_2340"# "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"


#########################################################
push!(LOAD_PATH,"src")
using SDMHelper
using PWAHelper

const BmatMC = let

    println("1. load data")
    @time PsiMC = read_precalc_basis(
        joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_mc.txt"));

    println("2. calculate sums")
    Nd, Nwaves = size(PsiMC)
    @time v = [sum(PsiMC[e,i]'*PsiMC[e,j] for e in 1:Nd)
        for i=1:Nwaves, j=1:Nwaves] ./ Nd;
    v
end

println("3. save the result")
write_cmatrix(BmatMC,
    joinpath(path_to_working_folder,"integrmat_$(mass_bin_name)_$(tslice)_mc.txt"))
