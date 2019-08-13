# parameters
mass_bin_name = "1540_1560" #ARGS[1]#
tslice = "t1"
path_to_working_folder = "/afs/cern.ch/user/m/mimikhas/work/private/pwa_data"#"data"
path_wavelist = "src"

######################################################
# @show ARGS
push!(LOAD_PATH,"src")
using amplitudes_compass
using PWAHelper
using SDMHelper

basisfunc_mc = joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_mc.bin")
const PsiMC = read_precalc_basis(basisfunc_mc);
Nd, Nwaves = size(PsiMC)
#
const Nbstrap_attempts = 100
const allBmatMC = [fill(0.0im,Nwaves,Nwaves) for b in 1:Nbstrap_attempts]
#
for b in 1:Nbstrap_attempts
    @show "Calculating integrals"
    randMCentries = rand(1:Nd,Nd) # here is the bootstrap
    @time allBmatMC[b] .= [sum(PsiMC[e,i]'*PsiMC[e,j] for e in randMCentries)
        for i=1:Nwaves, j=1:Nwaves] ./ Nd;
end
#

for (b,_BmatMC) in enumerate(allBmatMC)
    write_cmatrix(_BmatMC,
        joinpath(path_to_working_folder,"integrmat_$(mass_bin_name)_$(tslice)_mc_b$(b).txt"))
end

# using Plots
# using LinearAlgebra
# plot(real.([diag(BmatMC) diag(allBmatMC[1]) diag(allBmatMC[2])]))
# plot(real.([BmatMC[50,:] allBmatMC[1][50,:] allBmatMC[2][50,:]]))
