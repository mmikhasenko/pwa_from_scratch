using Plots
theme(:juno)
theme(:default)

push!(LOAD_PATH,"src")

using DalitzPlotAnalysis
using amplitudes_compass
using PWAHelper
using SDMHelper
using FittingPWALikelihood


mass_bin_name = "2320_2340"

for app in ["rd", "mc", "fu"]
    @eval $(Symbol("kinvar_"*app)) = joinpath("data","$(mass_bin_name)_t1_"*$app*".txt")
    @eval $(Symbol("basisfunc_"*app)) = joinpath("data","$(mass_bin_name)_t1_"*$app*".jld")
end

##### Converting data data #####

iffile("data/$(mass_bin_name)_rd.root") || error("Cannot find the file")
iffile("data/$(mass_bin_name)_mc.root") || error("Cannot find the file")
@time let tmin = "0.1", tmax = "0.112853"
    run(pipeline(`src/Extract data/$(mass_bin_name)_rd.root D $(tmin) $(tmax)`, stdout=kinvar_rd))
    run(pipeline(`src/Extract data/$(mass_bin_name)_mc.root M $(tmin) $(tmax)`, stdout=kinvar_mc))
    run(pipeline(`src/Extract data/$(mass_bin_name)_mc.root F $(tmin) $(tmax)`, stdout=kinvar_fu))
end

################################


@time precalculate_compass_basis(kinvar_rd, basisfunc_rd)
@time precalculate_compass_basis(kinvar_mc, basisfunc_mc)
@time precalculate_compass_basis(kinvar_fu, basisfunc_fu)

@time const PsiMC = read_precalc_basis(basisfunc_mc);

@time const sum_mat = [sum(PsiMC[e,i]'*PsiMC[e,j] for e in 1:size(PsiMC,1))
    for i=1:Nwaves, j=1:Nwaves] /size(PsiMC,1);
# 10.593515 seconds (101.00 k allocations: 4.752 MiB)

let sum_mat_n = [sum_mat[i,j]/sqrt(sum_mat[i,i]*sum_mat[j,j]) for i=1:Nwaves, j=1:Nwaves];
    heatmap(real(sum_mat_n))
end

const noϵ = [i==1 for i=1:size(wavesfile[:,6],1)]
const posϵ = [ϵ=="+" for ϵ in wavesfile[:,6]]
const negϵ = [ϵ=="-" for ϵ in wavesfile[:,6]]

const ModelBlocks = [noϵ, posϵ]
const Npar = size(get_parameter_map(ModelBlocks),2)

const PsiDT = read_precalc_basis(basisfunc_rd);

LLH, GRAD, LLH_and_GRAD!, HES = createLLHandGRAD(PsiDT, sum_mat, ModelBlocks);

test_t = rand(Npar)-0.5
@time @show LLH(test_t)
@time minpars = minimize(LLH, LLH_and_GRAD!; verbose=1, starting_pars=test_t)

minpars0 = vcat(readdlm("minpars_r1.txt")...);
@time minpars = minimize(LLH, LLH_and_GRAD!; verbose=1, starting_pars=minpars0)
writedlm("minpars_r1.txt", minpars)

@time minpars_tries = let Natt = 5
    tries = [];
    for i in 1:Natt
        _minpars = minimize(LLH, LLH_and_GRAD!; verbose=1, starting_pars=rand(Npar))
        _llh = LLH(_minpars)
        push!(tries,(_llh,_minpars))
    end
    tries
end

# three well pronounced minimas
@show [v[1] for v in minpars_tries]

minpars_first = minpars_tries[3][2]
writedlm("minpars_r2.txt", minpars_first)
minpars_second = minpars_tries[1][2]
writedlm("minpars_r2.txt", minpars_second)
minpars_third = minpars_tries[5][2]
writedlm("minpars_r3.txt", minpars_third)

#####################################################################################
#####################################################################################
@time hes_anal = HES(minpars)
inv_hes = inv(hes_anal)

diag_error = sqrt.(diag(inv_hes))
# diag(inv_hes);
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
# savefig(joinpath("plots","data_with_errors.png"))
# savefig(joinpath("plots","data_with_errors.pdf"))

#####################################################################################
#####################################################################################

@time const PsiFU = read_precalc_basis(basisfunc_fu);

@time const BmatFU = [sum(PsiFU[e,i]'*PsiFU[e,j] for e in 1:size(PsiFU,1))
    for i=1:Nwaves, j=1:Nwaves] /size(PsiMC,1);
# BmatFU
SDM = size(PsiDT,1)*pars_to_SDM(minpars, BmatFU, ModelBlocks)

minpars./SDM_to_pars(SDM/size(PsiDT,1), BmatFU, ModelBlocks)

SDM_RD = let path_to_SDM = "/localhome/mikhasenko/cernbox/tmp/pwa_from_scratch_data"
    readdlm(joinpath(path_to_SDM,"0.100000-0.112853","sdm92.re")) +
    1im*readdlm(joinpath(path_to_SDM,"0.100000-0.112853","sdm92.im"))
end

SDM_RD_err = let path_to_SDM = "/localhome/mikhasenko/cernbox/tmp/pwa_from_scratch_data"
    readdlm(joinpath(path_to_SDM,"0.100000-0.112853","sdm92-err.re")) +
    1im*readdlm(joinpath(path_to_SDM,"0.100000-0.112853","sdm92-err.im"))
end

minpars_FHDR = SDM_to_pars(SDM_RD/size(PsiDT,1), BmatFU, ModelBlocks)

plot([real.(diag(SDM_RD)) real.(diag(SDM))], lab=["F.H.-D.R." "M.M."],
    xlab="# wave", title = "Diagonal of the SDM", size=(800,500))
savefig(joinpath("plots","sdm_results_$(mass_bin_name)_t1.png"))
savefig(joinpath("plots","sdm_results_$(mass_bin_name)_t1.pdf"))
plot(-real.(diag(SDM))./real.(diag(SDM_RD)) .+ 1, lab=["1 - M.M. / F.H.-D.R."],
    xlab="# wave", title = "Diagonal of the SDM", size=(1800,500),
    xticks=[i for i in 1:88])
savefig(joinpath("plots","sdm_diff_$(mass_bin_name)_t1.pdf"))

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
@show [minpars_cv minpars_rd]
minpars_rd2 = [copysign(v,minpars[i]) for (i,v) in enumerate(minpars_rd)]
minpars

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

writedlm(joinpath("data","bootstrap_tmp","minpars_r1.txt"), minpars)
writedlm(joinpath("data","bootstrap_tmp","minpars_r2.txt"), minpars_second)
writedlm(joinpath("data","bootstrap_tmp","minpars_r3.txt"), minpars_third)
BootstrapResults = let Nb = 100
    # path to save
    path_to_tmp_res = joinpath("data","bootstrap_tmp")
    # path to save
    minpars0 = vcat(readdlm(joinpath(path_to_tmp_res,"minpars_r3.txt"))...);
    res = Matrix{Float64}(Nb,length(minpars0))
    llh = Vector{Float64}(Nb)
    @progress for b in 15:Nb
        Nd = size(PsiDT,1)
        const _PsiDT = hcat([PsiDT[rand(1:Nd),:] for i in 1:Nd]...).'
        # @show size(_PsiDT)
        _LLH, _GRAD, _LLH_and_GRAD!, _HES = createLLHandGRAD(_PsiDT, sum_mat, ModelBlocks);
        @time _minpars = minimize(_LLH, _LLH_and_GRAD!; verbose=1, starting_pars=minpars0)
        # @show minpars
        res[b,:] = _minpars
        llh[b] = _LLH(_minpars)
        writedlm(joinpath(path_to_tmp_res,"BootstrapResults-$(b).txt"), res[b,:])
        writedlm(joinpath(path_to_tmp_res,"llh-$(b).txt"), llh[b])
    end
    res
end
# BootstrapResults = let Nb = 100, path_to_tmp_res = joinpath("data","bootstrap_tmp")
#     res = Matrix{Float64}(Nb,Npar)
#     for b in 1:Nb
#         res[b,:] = readdlm(joinpath(path_to_tmp_res,"BootstrapResults-$(b).txt"))
#     end
#     res
# end
# writedlm(joinpath("data","bootstrap_tmp","bootstrap_r3.txt"), BootstrapResults)
# llh = let Nb = 100, path_to_tmp_res = joinpath("data","bootstrap_tmp")
#     res = Vector{Float64}(Nb)
#     for b in 1:Nb
#         res[b] = readdlm(joinpath(path_to_tmp_res,"llh-$(b).txt"))[1]
#     end
#     res
# end
# writedlm(joinpath("data","bootstrap_tmp","llh_r3.txt"), llh)
# histogram!(llh, bins=30, lab = "-LogLH")
# vline!([LLH(minpars_third)], lab = "main fit")
# savefig("plots/llh_three_minima.pdf")

vline!([minpars[1]])
histogram(BootstrapResults[:,1])
scatter(BootstrapResults[:,103], BootstrapResults[:,102])

btp_error = sqrt.([cov(BootstrapResults[:,i]) for i in 1:Npar])
btp_mean = [mean(BootstrapResults[:,i]) for i in 1:Npar]

plot([minpars btp_mean], lab=["main fit" "bootstrap mean"], xlab = "# parameter")

let m = max(diag_error...)
	for (i,v) in enumerate(diag_error)
		(v == m) && @show i
	end
end

let tog = [[minpars[i],diag_error[i]] for i in 1:length(minpars)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,1),
    lab="", title="PWA results and the errors Hessian", xlab="# parameter", size=(800,500))
end
savefig(joinpath("plots","parameters_with_hessian_errors.png"))
savefig(joinpath("plots","parameters_with_hessian_errors.pdf"))

let tog = [[btp_mean[i], btp_error[i]] for i in 1:length(minpars)]
    stog = sort(tog, by=x->abs(x[1]), rev=true)
    hcat_stog = hcat(stog...)
    plot!(abs.(hcat_stog[1,:]), yerr = hcat_stog[2,:], ylim=(0,1),
    lab="", title="PWA results and the errors Bootstrap", xlab="# parameter", size=(800,500))
end
savefig(joinpath("plots","parameters_with_bootstrap_errors.png"))
savefig(joinpath("plots","parameters_with_bootstrap_errors.pdf"))

###
SDM_BootstrapResults = hcat(
	[real.(diag(size(PsiDT,1)*pars_to_SDM(BootstrapResults[i,:], BmatFU, ModelBlocks)))
		for i in 1:size(BootstrapResults,1)]...);

plot(SDM_BootstrapResults, lab="")

let v = SDM_BootstrapResults
	p = sort([v[i,:] for i in 1:size(v,1)], by=x->mean(x), rev=true)
	quantiles = [quantile(p[i],[0.16,0.84]) for i in 1:length(p)]
	plot([quantiles[i][1] for i in 1:length(p)], fillrange=[quantiles[i][2] for i in 1:length(p)], lab="16%-84% quantiles", c=:lightgrey)
	plot!([mean(p[i]) for i in 1:length(p)], lab="Bootstrap mean")
	pm = sort(real.(diag(SDM)), rev=true)
	plot!(pm, lab="main fit", size=(1000,400))
end
let v = real.(diag(SDM_RD)), ve = real.(diag(SDM_RD_err))
	vs = sort([[v[i],ve[i]] for i in 1:length(v)], by=x->x[1], rev=true);
	plot!([i[1] for i in vs], yerr = [i[2] for i in vs], lab = "Florain-Dima")
end
savefig(joinpath("plots","compare_everything.pdf"))


let i = 8
	histogram(SDM_BootstrapResults[i,:], bins=20)
	vline!([(real(SDM[i,i]) + sqrt(cov(SDM_BootstrapResults[i,:]))*[-1,0,1])...])
end


max(abs.(minpars_second./minpars_first)...)

vcat([1,3,4]...,[4])

let v = BootstrapResults.'
	p0 = sort([vcat(v[i,:]..., minpars_first[i], minpars_second[i]) for i in 1:size(v,1)], by=x->mean(x), rev=true)
	# p = [v[i,:] for i in 1:size(v,1)]
	p = [p0[i][1:end-2] for i in 1:length(p0)]
	quantiles = [quantile(p[i][1:end-2],[0.16,0.84]) for i in 1:length(p)]
	plot([quantiles[i][1] for i in 1:length(p)], fillrange=[quantiles[i][2] for i in 1:length(p)], lab="16%-84% quantiles", c=:lightgrey)
	plot!([mean(p[i]) for i in 1:length(p)], lab="Bootstrap mean", size=(2000,800))
	plot!([p0[i][end-1] for i in 1:length(p0)], lab="deepest min")
	plot!([p0[i][end] for i in 1:length(p0)], lab="second min")
end

savefig(joinpath("plots","minimias.pdf"))
