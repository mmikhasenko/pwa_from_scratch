# parameters
mass_bin_name = ARGS[1]
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
Nbstrap_attempts = 100 # No. of bootstrap attempts

######################################################
using DelimitedFiles
push!(LOAD_PATH, "src")
using SDMHelper
using FittingPWALikelihood
using amplitudes_compass
using PWAHelper
using Random

###############################################################################
# Model setting

# read wavelist throw waves below threshol
wavelist = get_wavelist(joinpath(path_wavelist, "wavelist_formated.txt");
    path_to_thresholds=joinpath(path_wavelist, "thresholds_formated.txt"),
    M3pi=Meta.parse(mass_bin_name[1:4]) / 1000)

const Nwaves = size(wavelist, 1)

# Blocks in wavelist
const noϵ = [1] # flat wave
const posϵ = [i for (i, ϵ) in enumerate(wavelist[:, 6]) if ϵ == "+"]
const negϵ = [i for (i, ϵ) in enumerate(wavelist[:, 6]) if ϵ == "-"]

# Model description
const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]
const Npar = get_npars(ModelBlocks)

# integrals from MC data
BmatMC = read_cmatrix(
    joinpath(path_to_working_folder, "integrmat_$(mass_bin_name)_$(tslice)_mc.txt"));

# get scale
const sqrtBscale = [sqrt(real(BmatMC[i, i])) for i in 1:size(BmatMC, 1)];
const parscale = abs.(extnd(sqrtBscale, get_parameter_map(ModelBlocks, Nwaves)))
# normalize B
scaled_BmatMC = BmatMC ./ sqrtBscale ./ transpose(sqrtBscale)

# load precalculated data array
const PsiRD = read_precalc_basis(
    joinpath(path_to_working_folder, "functions_$(mass_bin_name)_$(tslice)_rd.bin"));

const Nd = size(PsiRD, 1)
# normalize Psi
const scaled_PsiRD = PsiRD
scaled_PsiRD ./= transpose(sqrtBscale)
# for i in 1:size(PsiRD,2)
#     scaled_PsiRD[:,i] .*= 1.0/sqrtBscale[i]
# end

###############################################################################
# load the best minimum
path_and_filename = split(readline("data/llh_attmpts_$(mass_bin_name)_$(tslice).txt"), "\t")[1]
bestminpars = readdlm(path_and_filename)
best_minimum = vcat(bestminpars...)
best_minimum .*= parscale

###############################################################################
final_pars_all_att = Matrix{Float64}(undef, Npar + 1, Nbstrap_attempts)
const scaled_pseudoPsiRD = copy(scaled_PsiRD)
#attempt at bootstrap
for b in 1:Nbstrap_attempts
    @show "Doing something"
    for e in 1:Nd
        scaled_pseudoPsiRD[e, :] .= scaled_PsiRD[rand(1:Nd), :]
    end
    # get functions to calculate LLH, derivative and hessian
    scaled_LLH, scaled_GRAD, scaled_LLH_and_GRAD!, scaled_HES =
        createLLHandGRAD(scaled_pseudoPsiRD, scaled_BmatMC, ModelBlocks)

    # fit
    @time minpars = minimize(scaled_LLH, scaled_LLH_and_GRAD!;
        algorithm=:LD_LBFGS, verbose=1, starting_pars=best_minimum,
        llhtolerance=1e-4)
    # calculate LLH
    final_pars_all_att[1, b] = scaled_LLH(minpars)
    # scale parameters back
    final_pars = minpars ./ parscale
    final_pars_all_att[2:end, b] .= final_pars
end

# save results
output_name = joinpath(path_to_working_folder, "llhfit_bootstrap_$(mass_bin_name)_$(tslice).txt")
println("Done! Saving to $(output_name) ...")
writedlm(output_name, final_pars_all_att)
