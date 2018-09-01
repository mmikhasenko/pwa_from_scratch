using Plots
theme(:default)

push!(LOAD_PATH,"src")

using DalitzPlotAnalysis
using amplitudes_compass
using PWAHelper
using SDMHelper
using FittingPWALikelihood

@time precalculate_compass_basis(joinpath("data","2300_2320_t1_rd.txt"), "rd.jld")
@time precalculate_compass_basis(joinpath("data","2300_2320_t1_mc.txt"), "mc.jld")
@time precalculate_compass_basis(joinpath("data","2300_2320_t1_fu.txt"), "fu.jld")

@time const PsiMC = read_precalc_basis("mc.jld");

@time const BmatMC = [sum(PsiMC[e,i]'*PsiMC[e,j] for e in 1:size(PsiMC,1))
    for i=1:Nwaves, j=1:Nwaves] /size(PsiMC,1);
# 10.593515 seconds (101.00 k allocations: 4.752 MiB)

let BmatMC_n = [BmatMC[i,j]/sqrt(BmatMC[i,i]*BmatMC[j,j]) for i=1:Nwaves, j=1:Nwaves];
    heatmap(real(BmatMC_n))
end

const noϵ = [i==1 for i=1:size(wavesfile[:,6],1)]
const posϵ = [ϵ=="+" for ϵ in wavesfile[:,6]]
const negϵ = [ϵ=="-" for ϵ in wavesfile[:,6]]

const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]
const Npar = size(get_parameter_map(ModelBlocks),2)

const PsiDT = read_precalc_basis("rd.jld");

LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiDT, BmatMC, ModelBlocks);

test_t = rand(Npar)
@time @show LLH(test_t)

minpars0 = vcat(readdlm("minpars_compass.txt")...);
# start nice algorithm which goes directly to the minimum
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_LBFGS, verbose=1, starting_pars=minpars0)
# start more precise algorithm
@time minpars = minimize(LLH, LLH_and_GRAD!;
    algorithm = :LD_SLSQP, verbose=1, starting_pars=minpars)
writedlm("minpars_compass.txt", minpars)

#####################################################################################
#####################################################################################

@time hes_anal = HES(minpars)
inv_hes = inv(hes_anal)
diag(inv_hes)

writedlm("invhes_compass.txt", inv_hes)

diag_error = sqrt.(diag(inv_hes))

@show wavenames[82]
plot(abs.(minpars[1:end-26]), yerr = diag_error[1:end-13], ylim=(0,.1))

let tog = [[minpars[i],diag_error[i]] for i in 1:(length(minpars)-26)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,.1))
end
savefig("plots/official_comass_fit_parameters_nonegref.pdf")

#####################################################################################
#####################################################################################

@time weights = let mpars = minpars,
    Tmap = get_parameter_map(ModelBlocks),
    pblocks = make_pblock_masks(ModelBlocks)
        [cohsq(extnd(PsiMC[i,:],Tmap).*mpars,pblocks) for i in 1:size(PsiMC,1)];
end;

const MC1 = readdlm(joinpath("data","2300_2320_t1_mc.txt"));
@time const MC3 = hcat([swap_kin_parameters(MC1[e,:]...) for e in 1:size(MC1,1)]...)';

histogram(sqrt.(vcat(MC1[:,2],MC3[:,2])), weights=vcat(weights,weights), bins=(linspace(0.3,2.2,100)))

const DT1 = readdlm(joinpath("data","2300_2320_t1_rd.txt"));
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

@time const PsiFU = read_precalc_basis("fu.jld");

@time const BmatFU = [sum(PsiFU[e,i]'*PsiFU[e,j] for e in 1:size(PsiFU,1))
    for i=1:Nwaves, j=1:Nwaves] /size(PsiMC,1);

SDM = size(PsiDT,1)*pars_to_SDM(minpars, BmatFU, ModelBlocks)

minpars./SDM_to_pars(SDM/size(PsiDT,1), BmatFU, ModelBlocks)

SDM_RD = let path = "/localhome/mikhasenko/cernbox/tmp/pwa_from_scratch_data"
    readdlm(joinpath(path,"0.100000-0.112853","sdm91.re")) +
    1im*readdlm(joinpath(path,"0.100000-0.112853","sdm91.im"))
end

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

diag(SDM)

minpars_rd = SDM_to_pars(SDM_RD/size(PsiDT,1), BmatFU, ModelBlocks)
minpars_cv = SDM_to_pars(SDM   /size(PsiDT,1), BmatFU, ModelBlocks)


arg(z) = atan2(imag(z), real(z))
SDM_with_sign = diag(SDM).*cis.([arg(SDM[2,i]) for i in 1:size(SDM,1)]);


gr()
let v = plot(size=(500,400))
    cutoff = 0.2*max(abs.(SDM_with_sign)...)
    for (i,v) in enumerate(SDM_with_sign)
        plot!([0, real(v)], [0, imag(v)], lab="", l=(1))
        (abs(v) > cutoff) && annotate!(
            [(0.75real(v), 0.75imag(v), text(wavenames[i][3:end],10))])
    end
    plot!()
end
savefig("/tmp/arrow_intensities.pdf")
writedlm("/tmp/SDM.2300.re",real(SDM))
writedlm("/tmp/SDM.2300.im",imag(SDM))

#####################################################################################
#####################################################################################

const testPsiDT = hcat([PsiDT[rand(1:size(PsiDT,1)),:] for i in 1:size(PsiDT,1)]...).'

BootstrapResults = let Nb = 200
    minpars0 = vcat(readdlm("minpars.txt")...);
    res = Matrix{Float64}(Nb,length(minpars0))
    llh = Vector{Float64}(Nb)
    @progress for b in 1:Nb
        Nd = size(PsiDT,1)
        const _PsiDT = hcat([PsiDT[rand(1:Nd),:] for i in 1:Nd]...).'
        @show size(_PsiDT)
        _LLH, _GRAD, _LLH_and_GRAD!, _HES = createLLHandGRAD(_PsiDT, BmatMC, ModelBlocks);
        @time minpars = minimize(_LLH, _LLH_and_GRAD!; verbose=1, starting_pars=minpars0)
        # @show minpars
        res[b,:] = minpars
        llh[b] = _LLH(minpars)
        writedlm("BootstrapResults-interm.txt", res)
        writedlm("llh.txt", llh)
    end
    res
end
BootstrapResults = readdlm("BootstrapResults.txt")

minpars
histogram(BootstrapResults[:,1])
scatter(BootstrapResults[:,10], BootstrapResults[:,40])

btp_error = sqrt.([cov(BootstrapResults[:,i]) for i in 1:Npar])
btp_mean = [mean(BootstrapResults[:,i]) for i in 1:Npar]

plot(hcat(minpars,btp_mean), lab=["main fit" "bootstrap mean"], xlab = "# parameter")

let tog = [[minpars[i],diag_error[i]] for i in 1:length(minpars)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,1),
    lab="", title="PWA results and the errors Hessian", xlab="# parameter", size=(800,500))
end
savefig(joinpath("plots","parameters_with_hessian_errors.png"))
savefig(joinpath("plots","parameters_with_hessian_errors.pdf"))

let tog = [[btp_mean[i],btp_error[i]] for i in 1:length(minpars)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,1),
    lab="", title="PWA results and the errors Bootstrap", xlab="# parameter", size=(800,500))
end
savefig(joinpath("plots","parameters_with_bootstrap_errors.png"))
savefig(joinpath("plots","parameters_with_bootstrap_errors.pdf"))
