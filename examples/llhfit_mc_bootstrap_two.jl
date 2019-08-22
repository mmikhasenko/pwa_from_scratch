# parameters
mass_bin_name = ARGS[1]
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
Nbstrap_attempts = 100 #ARGS[2]

######################################################
using DelimitedFiles
push!(LOAD_PATH,"src")
using SDMHelper
using FittingPWALikelihood
using amplitudes_compass
using PWAHelper
using Random
rng = MersenneTwister()
serial= round(Int,rand(Random.seed!(rng))*100000)
print(serial)

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

# load precalculated data array
const PsiRD = read_precalc_basis(
    joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_rd.bin"));
const _PsiRD = copy(PsiRD)

#  _|          _|_|      _|_|    _|_|_|
#  _|        _|    _|  _|    _|  _|    _|
#  _|        _|    _|  _|    _|  _|_|_|
#  _|        _|    _|  _|    _|  _|
#  _|_|_|_|    _|_|      _|_|    _|

###############################################################################
# load the best minimum
path_and_filename = split(readline("data/llh_attmpts_$(mass_bin_name)_$(tslice).txt"),"\t")[1]
best_minimum = vcat(readdlm(path_and_filename)...)

###############################################################################
final_pars_all_att = Matrix{Float64}(undef, Npar+1, Nbstrap_attempts)
#attempt at bootstrap
for b in 1:Nbstrap_attempts
    @show "Doing something", b
    # integrals from MC data
    _BmatMC = read_cmatrix(
        joinpath(path_to_working_folder,"integrmat_$(mass_bin_name)_$(tslice)_mc_b$(b).txt"));

    # normalize B
    _Bscale = [_BmatMC[i,i] for i in 1:size(_BmatMC,1)];
    for i=1:Nwaves, j=1:Nwaves
        _BmatMC[i,j] *= 1/sqrt(_Bscale[i]*_Bscale[j])
    end

    # normalize Psi
    _PsiRD .= PsiRD
    for i in 1:size(_PsiRD,2)
        _PsiRD[:,i] .*= 1.0/sqrt(_Bscale[i])
    end

    # get functions to calculate LLH, derivative and hessian
    LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(_PsiRD, _BmatMC, ModelBlocks);
    # parameters
    _parscale = abs.(extnd(sqrt.(_Bscale), get_parameter_map(ModelBlocks, Nwaves)))
    _best_minimum = best_minimum .* _parscale
    #
    @show LLH(_best_minimum)
    # fit
    # start vast algorithm which goes directly to the minimum
    @time minpars = minimize(LLH, LLH_and_GRAD!;
        algorithm = :LD_LBFGS, verbose=1, starting_pars=_best_minimum,
        llhtolerance = 1e-4)
    # scale parameters back
    final_pars = minpars ./ _parscale
    #
    final_pars_all_att[1,b] = LLH(minpars)
    final_pars_all_att[2:end,b] .= final_pars
    #
end


# save results
output_name = joinpath(path_to_working_folder,"llhfit_mc_bootstrap_$(mass_bin_name)_$(tslice).txt")
println("Done! Saving to $(output_name) ...")
writedlm(output_name, final_pars_all_att)
