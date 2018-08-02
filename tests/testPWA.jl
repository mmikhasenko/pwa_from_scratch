using Plots
theme(:juno)

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

@time const sum_mat = [sum(PsiMC[e,i]'*PsiMC[e,j] for e in 1:size(PsiMC,1))
    for i=1:88, j=1:88] /size(PsiMC,1);
# 10.593515 seconds (101.00 k allocations: 4.752 MiB)

let sum_mat_n = [sum_mat[i,j]/sqrt(sum_mat[i,i]*sum_mat[j,j]) for i=1:88, j=1:88];
    heatmap(real(sum_mat_n))
end
# waves = readdlm(pwd()*"/src/wavelist_formated.txt");
waves = readdlm(joinpath("src","wavelist_formated.txt"));

const noϵ = [i==1 for i=1:size(waves[:,6],1)]
const posϵ = [ϵ=="+" for ϵ in waves[:,6]]
const negϵ = [ϵ=="-" for ϵ in waves[:,6]]

const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]

const PsiDT = read_precalc_basis("rd.jld");

LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiDT, sum_mat, ModelBlocks);

test_t = rand(186)
@time @show LLH(test_t)

minpars0 = vcat(readdlm("minpars.txt")...);
@time minpars = minimize(LLH, LLH_and_GRAD!; verbose=1, starting_pars=minpars0)
writedlm("minpars.txt", minpars)

#####################################################################################
#####################################################################################

@time hes_anal = HES(minpars)
inv_hes = inv(hes_anal)

diag_error = sqrt.(diag(inv_hes))

let tog = [[minpars[i],diag_error[i]] for i in 1:length(minpars)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,1))
end

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
    for i=1:88, j=1:88] /size(PsiMC,1);

SDM = size(PsiDT,1)*pars_to_SDM(minpars, BmatFU, ModelBlocks)

minpars./SDM_to_pars(SDM/size(PsiDT,1), BmatFU, ModelBlocks)

SDM_RD = readdlm(joinpath("SDMs","0.100000-0.112853","sdm91.re")) +
    1im*readdlm(joinpath("SDMs","0.100000-0.112853","sdm91.im"))

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
        _LLH, _GRAD, _LLH_and_GRAD!, _HES = createLLHandGRAD(_PsiDT, sum_mat, ModelBlocks);
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

btp_error = sqrt.([cov(BootstrapResults[:,i]) for i in 1:186])
btp_mean = [mean(BootstrapResults[:,i]) for i in 1:186]

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
