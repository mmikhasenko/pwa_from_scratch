# parameters
mass_bin_name = ARGS[1]# "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"


#########################################################
@show ARGS
push!(LOAD_PATH, "src")
using SDMHelper
using PWAHelper

const BmatFU = let

    println("1. load data")
    @time PsiFU = read_precalc_basis(
        joinpath(path_to_working_folder, "functions_$(mass_bin_name)_$(tslice)_fu.txt"))

    println("2. calculate sums")
    PsiTFU = transpose(PsiFU)
    @views v = sum(PsiTFU[:, e] * PsiTFU[:, e]' for e in 1:size(PsiTFU, 2)) ./ size(PsiFU, 1)
    conj(v)
    v
end

println("3. save the result")
write_cmatrix(BmatFU,
    joinpath(path_to_working_folder, "integrmat_$(mass_bin_name)_$(tslice)_fu.txt"))
