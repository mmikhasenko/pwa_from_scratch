# parameters
mass_bin_name = ARGS[1]
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
    joinpath(path_to_working_folder,"functions_$(mass_bin_name)_$(tslice)_rd.bin"));
# normalize Psi
for i in 1:size(PsiRD,2)
    PsiRD[:,i] .*= 1.0/sqrt(Bscale[i])
end

# get functions to calculate LLH, derivative and hessian
LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiRD, BmatMC, ModelBlocks);

# generate random parameters on a surface of intensity constraint
const Npar = get_npars(ModelBlocks)
const pars0 = rand(Npar);
pars0 .*= get_parameter_ranges(BmatMC, ModelBlocks)
normalize_pars!(pars0, BmatMC, ModelBlocks)

# fit
# start vast algorithm which goes directly to the minimum
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_LBFGS, verbose=1, starting_pars=pars0,
    llhtolerance = 1e-4)

# scale parameters back
parscale = abs.(extnd(sqrt.(Bscale), get_parameter_map(ModelBlocks, Nwaves)))
final_pars = minpars ./ parscale

println("Done! Saving...")
writedlm(joinpath(path_to_working_folder,"llhfit_$(mass_bin_name)_$(tslice)_$(serial).txt"), final_pars)
