using Plots
theme(:default)

using DelimitedFiles

push!(LOAD_PATH,"src")
using DalitzPlotAnalysis
using amplitudes_compass
using PWAHelper
using SDMHelper
using FittingPWALikelihood

mass_bin_name = "2320_2340"# "1540_1560"
# set names
for app in ["rd", "mc", "fu"]
    @eval $(Symbol("kinvar_"*app)) = joinpath("data","variables_$(mass_bin_name)_t1_"*$app*".txt")
    @eval $(Symbol("basisfunc_"*app)) = joinpath("data","functions_$(mass_bin_name)_t1_"*$app*".txt")
end

##### Converting data data #####

isfile("data/$(mass_bin_name)_rd.root") || error("Cannot find the file")
isfile("data/$(mass_bin_name)_mc.root") || error("Cannot find the file")
@time let tmin = "0.1", tmax = "0.112853"
    run(pipeline(`src/Extract data/$(mass_bin_name)_rd.root D $(tmin) $(tmax)`, stdout=kinvar_rd))
    run(pipeline(`src/Extract data/$(mass_bin_name)_mc.root M $(tmin) $(tmax)`, stdout=kinvar_mc))
    run(pipeline(`src/Extract data/$(mass_bin_name)_mc.root F $(tmin) $(tmax)`, stdout=kinvar_fu))
end

################################

# buildbasis(joinpath(pwd(),"src/wavelist_formated.txt"); M3pi=parse(mass_bin_name[1:4])/1000,
#     path_to_thresholds=joinpath(pwd(),"src/thresholds_formated.txt"))
wavelist = get_wavelist(joinpath(pwd(),"src/wavelist_formated.txt");
    path_to_thresholds=joinpath(pwd(),"src/thresholds_formated.txt"),
    M3pi=parse(mass_bin_name[1:4])/1000)

Nwaves = size(wavelist,1)

wavenames = get_wavenames(wavelist)
wavebasis = get_wavebasis(wavelist)

@time for i=1:10000
    compass_jmels_basis_psi(QNs=(1,0,1*true,0,1),τ1=(1.1,rand(),rand(),rand(),rand()),s=1.5)
end

@profiler for i=1:10000
    compass_jmels_basis_psi(QNs=(1,0,1*true,0,1),τ1=(1.1,rand(),rand(),rand(),rand()),s=1.5)
end

@time precalculate_compass_basis(wavebasis, kinvar_rd, basisfunc_rd)
@time precalculate_compass_basis(wavebasis, kinvar_mc, basisfunc_mc)
@time precalculate_compass_basis(wavebasis, kinvar_fu, basisfunc_fu)

isfile("src/wavelist_formated.txt")

@time const PsiMC = read_precalc_basis(basisfunc_mc);

@time const BmatMC = let
    v = [sum(PsiMC[e,i]'*PsiMC[e,j] for e in 1:size(PsiMC,1))
        for i=1:Nwaves, j=1:Nwaves] /size(PsiMC,1);
    v
end
write_SDM(BmatMC, "BmatMC_$(mass_bin_name).txt")
BmatMC = read_SDM("BmatMC_$(mass_bin_name).txt")

# 10.593515 seconds (101.00 k allocations: 4.752 MiB)

let BmatMC_n = [BmatMC[i,j]/sqrt(BmatMC[i,i]*BmatMC[j,j]) for i=1:Nwaves, j=1:Nwaves];
    heatmap(real(BmatMC_n))
end

const noϵ = [i==1 for i=1:size(wavesfile[:,6],1)]
const posϵ = [ϵ=="+" for ϵ in wavesfile[:,6]]
const negϵ = [ϵ=="-" for ϵ in wavesfile[:,6]]

const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]
const Npar = size(get_parameter_map(ModelBlocks),2)

const PsiRD = read_precalc_basis(basisfunc_rd);

LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiRD, BmatMC, ModelBlocks);

test_t = rand(Npar)
@time @show LLH(test_t)
@time @show LLH(test_t)

let scale = get_intesity(test_t, BmatMC, ModelBlocks)
    get_intesity(test_t ./ sqrt(scale), BmatMC, ModelBlocks)
end

minpars0 = vcat(readdlm("minpars_compass_$(mass_bin_name).txt")...)[1:Npar];
# start nice algorithm which goes directly to the minimum
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_LBFGS, verbose=1, starting_pars=minpars0)
# start more precise algorithm
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_SLSQP, verbose=1, starting_pars=minpars)
writedlm("minpars_compass_$(mass_bin_name).txt", minpars)

#####################################################################################
#####################################################################################

@time hes_anal = HES(minpars)
inv_hes = inv(hes_anal)
diag(inv_hes)

writedlm("invhes_compass_$(mass_bin_name).txt", inv_hes)

diag_error = sqrt.(diag(inv_hes))

@show wavenames[82]
plot(abs.(minpars[1:end-26]), yerr = diag_error[1:end-13], ylim=(0,.1))

let tog = [[minpars[i],diag_error[i]] for i in 1:(length(minpars)-26)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,.1))
end
savefig(joinpath("plots","official_comass_fit_parameters_nonegref.pdf"))

#####################################################################################
#####################################################################################

intensity_ranges_per_parameter = [begin
        pars = fill(0.0,Npar)
        pars[i] = 1.0
        get_intesity(pars, BmatMC, ModelBlocks)
    end for i=1:Npar]
limits = 1./sqrt.(intensity_ranges_per_parameter)

function get_rand_starting_values(limits)
    return [(2*rand()-1)*l for l in limits]
end

# start from random points on the contrained surface
for e in 13:100
    rand_pars =  get_rand_starting_values(1.5.*limits)
    normalize_pars!(rand_pars, BmatMC, ModelBlocks)
    @time pars_in_min = minimize(LLH, LLH_and_GRAD!;
        algorithm = :LD_LBFGS, verbose=1, starting_pars=rand_pars)
    writedlm("data/studies_of_minimas/m$(e).txt", pars_in_min)
end

minima = hcat([readdlm("data/studies_of_minimas/m$(e).txt")[:,1] for e in 1:95]...)
minimaLLH = [LLH(minima[:,i]) for i=1:size(minima,2)]

collect(1:95)[minimaLLH.-LLH(minpars0) .> 5]
histogram(minimaLLH.-LLH(minpars0), bins=linspace(-2,5,10))

#####################################################################################
#####################################################################################


@time weights = let mpars = minpars,
    Tmap = get_parameter_map(ModelBlocks),
    pblocks = make_pblock_masks(ModelBlocks)
        [cohsq(extnd(PsiMC[i,:],Tmap).*mpars,pblocks) for i in 1:size(PsiMC,1)];
end;

const MC1 = readdlm(kinvar_mc);
@time const MC3 = hcat([swap_kin_parameters(MC1[e,:]...) for e in 1:size(MC1,1)]...)';

histogram(sqrt.(vcat(MC1[:,2],MC3[:,2])), weights=vcat(weights,weights), bins=(linspace(0.3,2.2,100)))

const DT1 = readdlm(kinvar_rd);
@time const DT3 = hcat([swap_kin_parameters(DT1[e,:]...) for e in 1:size(DT1,1)]...)';

histogram2d(DT1[:,2], DT3[:,2], bins=linspace(0,6,100))

stephist(sqrt.(vcat(MC1[:,2],MC3[:,2])), weights=vcat(weights,weights)/size(MC1,1)*size(DT1,1),
    bins=linspace(0.3,2.2,100),lab="weighted MC", size=(800,500), xlab="M3pi (GeV)")
let v = stephist!(sqrt.(vcat(DT1[:,2],DT3[:,2])), bins=linspace(0.3,2.2,100),lab="data")
    yarr = v.subplots[1][2][:y]
    xarr = v.subplots[1][2][:x]
    @show length(xarr)
    x2 = [(xarr[2i]+xarr[2(i+1)])/2 for i in 1:(div(length(xarr),2)-1)];
    y2 = [yarr[2i] for i in 1:(div(length(yarr),2)-1)];
    plot!(x2, y2, err=sqrt.(y2), lab="", marker = (0, stroke(0.5, :gray)), l=nothing)
end
savefig(joinpath("plots","data_with_errors.png"))
savefig(joinpath("plots","data_with_errors.pdf"))

#####################################################################################
#####################################################################################



@time const BmatFU = let
    @time const PsiFU = read_precalc_basis(basisfunc_fu);
    [sum(PsiFU[e,i]'*PsiFU[e,j] for e in 1:size(PsiFU,1))
        for i=1:Nwaves, j=1:Nwaves] /size(PsiMC,1)
end
write_SDM(BmatFU, "BmatFU_$(mass_bin_name).txt")
BmatFU = read_SDM("BmatFU_$(mass_bin_name).txt")

SDM = size(PsiRD,1)*pars_to_SDM(minpars, BmatFU, ModelBlocks)
SDM1 = size(PsiRD,1)*pars_to_SDM(minima[:,1], BmatFU, ModelBlocks)
SDM9 = size(PsiRD,1)*pars_to_SDM(minima[:,9], BmatFU, ModelBlocks)
SDM26 = size(PsiRD,1)*pars_to_SDM(minima[:,26], BmatFU, ModelBlocks)

minpars./SDM_to_pars(SDM/size(PsiRD,1), BmatFU, ModelBlocks)
SDM_RD = let path = "/localhome/mikhasenko/cernbox/tmp/pwa_from_scratch_data"
    v = readdlm(joinpath(path,"0.100000-0.112853","sdm53.re")) +
        1im*readdlm(joinpath(path,"0.100000-0.112853","sdm53.im"))
    v[thresholds_filter,thresholds_filter]
end

SDM_RD_err = let path = "/localhome/mikhasenko/cernbox/tmp/pwa_from_scratch_data"
    v = readdlm(joinpath(path,"0.100000-0.112853","sdm53-err.re")) +
        1im*readdlm(joinpath(path,"0.100000-0.112853","sdm53-err.im"))
    v[thresholds_filter,thresholds_filter]
end

minpars_rd = SDM_to_pars(SDM_RD/size(PsiRD,1), BmatFU, ModelBlocks)
@time minpars_from_rd = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_LBFGS, verbose=1, starting_pars=minpars_rd)
minpars = minpars_from_rd

plot(hcat(real.(diag(SDM_RD)), real.(diag(SDM))), lab=["F.H.-D.R." "M.M."],
    xlab="# wave", title = "Diagonal of the SDM", size=(800,500))
savefig(joinpath("plots","sdm_results.png"))
savefig(joinpath("plots","sdm_results.pdf"))

let rd=real.(diag(SDM_RD)), wn = wavenames
    sort([(wn[i],rd[i]) for i in 1:length(rd)], by=x->x[2],rev=true)
end
let
    plot(real.(SDM_RD[2,:]), frame=:origin)
    plot!(real.(SDM[2,:]))
end
plot(hcat(real.(SDM_RD[2,:]), real.(SDM[2,:])), lab=["F.H.-D.R." "M.M."],
    xlab="# wave", title = "Real part of the second row of SDM", size=(800,500))
savefig(joinpath("plots","sdm2_results.png"))
savefig(joinpath("plots","sdm2_results.pdf"))

sum(diag(SDM))
sum(diag(SDM_RD))

plotPolarSDM(SDM, cutoff_scale=0.55)
savefig(joinpath("data","arrow_intensities.pdf"))

#####################################################################################
#####################################################################################
v = let Nb = 500
    # path to save
    path_to_tmp_res = joinpath("data","bootstrap_tmp")
    minpars0 = vcat(readdlm("minpars_compass_$(mass_bin_name).txt")...);
    # res = Matrix{Float64}(Nb,length(minpars0))
    # llh = Vector{Float64}(Nb)
    # b = 16
    @parallel for b in 1:Nb
        @show "Doing something"
        Nd = size(PsiRD,1)
        const _PsiRD = hcat([PsiRD[rand(1:Nd),:] for i in 1:Nd]...).'
        # @show size(_PsiRD)
        _LLH, _GRAD, _LLH_and_GRAD!, _HES = FittingPWALikelihood.createLLHandGRAD(_PsiRD, BmatMC, ModelBlocks);
        try
            @time _minpars = FittingPWALikelihood.minimize(_LLH, _LLH_and_GRAD!;
                verbose=0, starting_pars=minpars0) #,algorithm=:LD_MMA
            writedlm(joinpath(path_to_tmp_res,"BootstrapResults-$(b).txt"), _minpars)

            writedlm(joinpath(path_to_tmp_res,"llh-$(b).txt"), _LLH(_minpars))
        catch TheException
            @show TheException
        end
    end
end
# fetch(v[2])
# BootstrapResults = let Nb = 500, path_to_tmp_res = joinpath("data","bootstrap_tmp")
#     res = []
#     for b in 1:Nb
#         file = joinpath(path_to_tmp_res,"BootstrapResults-$(b).txt")
#         isfile(file) && push!(res, readdlm(file))
#     end
#     hcat(res...)'
# end
# writedlm(joinpath("data","bootstrap_tmp","bootstrap_main_$(mass_bin_name).txt"), BootstrapResults)
# llh = let Nb = 578, path_to_tmp_res = joinpath("data","bootstrap_tmp")
#     res = []
#     for b in 1:Nb
#         file = joinpath(path_to_tmp_res,"llh-$(b).txt")
#         isfile(file) && push!(res, readdlm(file))
#     end
#     hcat(res...)'[:,1]
# end
# writedlm(joinpath("data","bootstrap_tmp","bootstrap_llh_$(mass_bin_name).txt"), llh)
@time for b in 1:size(BootstrapResults,1)
    _br = BootstrapResults[b,:]
    _sdm = size(PsiRD,1)*pars_to_SDM(_br, BmatFU, ModelBlocks)
    path_to_tmp_res = joinpath("data","bootstrap_tmp","SDMs-1540_1560")
    writedlm(joinpath(path_to_tmp_res,"rb-$(b)-sdm.txt"), [real(_sdm) imag(_sdm)])
end
SDMs = [read_SDM(joinpath("data","bootstrap_tmp","SDMs_$(mass_bin_name)","rb-$(i)-sdm.txt"))
    for i in 1:497]

#####################################################################################
#####################################################################################

histogram(BootstrapResults[:,2], bins=100)
histogram(BootstrapResults[:,1], bins=100)

using PlotHelper
plotBSTsample(2, SDMs, SDM, SDM_RD, SDM_RD_err)
arg(z) = atan2(imag(z),real(z))
plotBSTsample(x->180/π*arg(x[2,15]*cis(π/2))-90, SDMs, SDM, SDM_RD, SDM_RD_err)
plotBSTsample(x->180/π*arg(x[2,15]), SDMs, SDM, SDM_RD, SDM_RD_err)

writedlm("/tmp/phase215.txt",broadcast(x->180/π*arg(x[2,15]*cis(π/2))-90, SDMs))

combres = vcat([constract_values(i, SDMs, SDM, SDM_RD, SDM_RD_err) for i in 1:Nwaves]...);
plotBSTsummary(combres) #  tosort=true, toannotate=true
savefig(joinpath("plots","bootstrap_combined_$(mass_bin_name)_t1.pdf"))

saveBSTSplot(joinpath("plots","bootstrap_summary_2300_t1.pdf"),
    SDMs, SDM, SDM_RD, SDM_RD_err)

##############################################################################
##############################################################################

plot(hcat(real.(diag(SDM_RD)), real.(diag(SDM1)), real.(diag(SDM9))),
    lab=["official" "min1" "min9"], lc=[:red :orange :green],
    xlab="# wave", ylab = "Diagonal of the SDM")
savefig("plots/two_minima_diagonal.pdf")

histogram([real(s[13,13]) for s in SDMs], bins=100, xlab=wavenames[13][3:end],
    c=:grey, lab="bootstrap")
vline!(broadcast(x->real(x[13,13]),[SDM_RD, SDM1, SDM9])',
    lab=["official" "min1" "min9"], lw = 2, lc=[:red :orange :green])
savefig("plots/two_minimas_sdm13.pdf")

histogram(broadcast(x->180/π*arg(x[2,15]*cis(π/2))-90+360,SDMs), bins=100,
    c=:grey, xlab="1++0+: f0(980)piP - rhopiS (deg.)", ylab = "Bootstrap Entries",
    lab="bootstrap")
vline!(broadcast(x->180/π*arg(x[2,15]*cis(π/2))-90+360,[SDM_RD, SDM1, SDM9])',
    lab=["official" "min1" "min9"], lw = 2, lc=[:red :orange :green])
savefig("plots/two_minimas_sdm215.pdf")
