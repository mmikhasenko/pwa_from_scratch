# parameters
mass_bin_name = ARGS[1]
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
Nbstrap_attempts = 5 # No. of bootstrap attempts

######################################################
using DelimitedFiles
push!(LOAD_PATH,"src")
using SDMHelper
using FittingPWALikelihood
using amplitudes_compass
using PWAHelper
using Random

###############################################################################
# Model setting

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
const Npar = get_npars(ModelBlocks)

# integrals from MC data
BmatMC = read_cmatrix(
    joinpath(path_to_working_folder,"integrmat_$(mass_bin_name)_$(tslice)_mc.txt"));

# normalize B
Bscale = [BmatMC[i,i] for i in 1:size(BmatMC,1)];
for i=1:Nwaves, j=1:Nwaves
    BmatMC[i,j] *= 1/sqrt(Bscale[i]*Bscale[j])
end
# load precalculated data array
const PsiRD = read_precalc_basis(
    joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_rd.txt"));
# normalize Psi
for i in 1:size(PsiRD,2)
    PsiRD[:,i] .*= 1.0/sqrt(Bscale[i])
end
const Nd = size(PsiRD,1)

###############################################################################
# load the best minimum
best_minimum = rand(Npar) # this line to be changed!

###############################################################################
final_pars_all_att = Array{Float64}(Npar+1, Nbstrap_attempts)
const pseudoPsiRD = copy(PsiRD)
#attempt at bootstrap
for b in 1:Nbstrap_attempts
    @show "Doing something"
    for e in 1:Nd
        pseudoPsiRD[e,:] .= PsiRD[rand(1:Nd),:]
    end
    # get functions to calculate LLH, derivative and hessian
    LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(pseudoPsiRD, BmatMC, ModelBlocks);

    # fit
    @time minpars = minimize(LLH, LLH_and_GRAD!;
        algorithm = :LD_LBFGS, verbose=1, starting_pars=best_minimum,
        llhtolerance = 1e-4)

    # scale parameters back
    parscale = abs.(extnd(sqrt.(Bscale), get_parameter_map(ModelBlocks, Nwaves)))
    final_pars = minpars ./ parscale
    # set up
    final_pars_all_att[1,b] = LLH(final_pars)
    final_pars_all_att[2:end,b] .= final_pars
end

# save results
output_name = joinpath(path_to_working_folder,"llhfit_bootstrap_$(mass_bin_name)_$(tslice).txt")
println("Done! Saving to $(output_name) ...")
writedlm(output_name, final_pars_all_att)
