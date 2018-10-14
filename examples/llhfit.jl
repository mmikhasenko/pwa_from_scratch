# parameters
mass_bin_name = "2320_2340"# "1540_1560"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"


######################################################
using DelimitedFiles
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
    joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_rd.txt"));

# get functions to calculate LLH, derivative and hessian
LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiRD, BmatMC, ModelBlocks);

# generate random parameters on a surface of intensity constraint
const Npar = get_npars(ModelBlocks)
pars0 = rand(Npar);
pars0 .*= get_parameter_ranges(BmatMC, ModelBlocks)
normalize_pars!(pars0, BmatMC, ModelBlocks)

# fit
# start vast algorithm which goes directly to the minimum
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_LBFGS, verbose=1, starting_pars=pars0)
# start more precise algorithm
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_SLSQP, verbose=1, starting_pars=minpars)

writedlm("llhfit_$(mass_bin_name)_$(tslice).txt", minpars)
