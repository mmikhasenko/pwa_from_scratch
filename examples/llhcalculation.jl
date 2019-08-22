# parameters
mass_bin_name = ARGS[1]
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
list_of_files = ARGS[2:end]

println("Starting with ARGS = ", ARGS)
######################################################
using DelimitedFiles
using LinearAlgebra

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
    joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_rd.bin"));

# get functions to calculate LLH, derivative and hessian
LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiRD, BmatMC, ModelBlocks);

println("---> Start calculating LLHs")
llhs = [LLH(vcat(readdlm(path_and_filename)...)) for path_and_filename in list_of_files]
sorted_list = sort(collect(zip(list_of_files, llhs)), by=x->x[2])
writedlm(joinpath(path_to_working_folder,"llh_attmpts_$(mass_bin_name)_$(tslice).txt"), sorted_list)


###########################################################################
println("---> Trying to improve minimum")
# fit the best minimum again
## normalize B
Bscale = [BmatMC[i,i] for i in 1:size(BmatMC,1)];
for i=1:size(BmatMC,1), j=1:size(BmatMC,2)
    BmatMC[i,j] *= 1/sqrt(Bscale[i]*Bscale[j])
end
## normalize Psi
for i in 1:size(PsiRD,2)
    PsiRD[:,i] .*= 1.0/sqrt(Bscale[i])
end
## get functions to calculate LLH, derivative and hessian
LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiRD, BmatMC, ModelBlocks);
## load the best parameters
best_parameters = vcat(readdlm(sorted_list[1][1])...)
parscale = abs.(extnd(sqrt.(Bscale), get_parameter_map(ModelBlocks, size(BmatMC,1))))
best_parameters_scalled = best_parameters .* parscale
## fit
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_LBFGS, verbose=1, starting_pars=best_parameters_scalled,
    llhtolerance = 1e-6)

max_reldiff_of_parameters = max(abs.(1 .- minpars ./ best_parameters_scalled)...)
@show max_reldiff_of_parameters
(max_reldiff_of_parameters > 0.1) && warn("Too big improvement in the fit while fine tuning!")

@time minpars2 = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_SLSQP, verbose=1, starting_pars=minpars,
    llhtolerance = 1e-6)

max_reldiff_of_parameters2 = max(abs.(1 .- minpars2 ./ minpars)...)
@show max_reldiff_of_parameters2
(max_reldiff_of_parameters2 > 0.1) && warn("Too big improvement in the fit while fine tuning!")

###########################################################################
println("---> Calculating Hessian")
# get derivative
hes = HES(minpars2)
singular = (det(hes) ≈ 0.0); singular && warn("Hessian matrix is singular!")
invhes = !(singular) ? inv(hes) : fill(0.0,size(hes)...)
invhes_scaled_back = invhes ./ transpose(parscale) ./ parscale
writedlm(joinpath(path_to_working_folder,"invhes_$(splitdir(sorted_list[1][1])[end])"), invhes_scaled_back)
println("Done")
