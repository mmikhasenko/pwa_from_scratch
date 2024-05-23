# In this file I compare three values for the model integrals

push!(LOAD_PATH, "src")
using DalitzPlotAnalysis
using amplitudes_compass
using PWAHelper

# calculate integrals with the acceptance
@time precalculate_compass_basis(joinpath("data", "2300_2320_t1_mc.txt"), "mc.jld")
@time const PsiMC = read_precalc_basis("mc.jld");
@time const BmatMC = [sum(PsiMC[e, i]' * PsiMC[e, j] for e in 1:size(PsiMC, 1))
                      for i = 1:Nwaves, j = 1:Nwaves] / size(PsiMC, 1);

# calculate integrals on phase space
@time precalculate_compass_basis(joinpath("data", "2300_2320_t1_fu.txt"), "fu.jld")
@time const PsiFU = read_precalc_basis("fu.jld");
@time const BmatFU = [sum(PsiFU[e, i]' * PsiFU[e, j] for e in 1:size(PsiFU, 1))
                      for i = 1:Nwaves, j = 1:Nwaves] / size(PsiFU, 1);

# same integral analytically
using Cuba
d2_analy = let s = 2.31^2
    function diagint(ind)
        function integrand(x, f)
            σ1 = 4mπ2 + x[1] * ((√s - mπ)^2 - 4mπ2)
            cosθ1, ϕ1 = 2 * x[2] - 1, π * (2 * x[3] - 1)
            cosθ23, ϕ23 = 2 * x[4] - 1, π * (2 * x[5] - 1)
            val = COMPASS_wave(ind, s, σ1, cosθ1, ϕ1, cosθ23, ϕ23)
            f[1] = abs2(val) * sqrt(λ(σ1, mπ2, mπ2) * λ(s, σ1, mπ2)) / (2π * s * σ1) * ((√s - mπ)^2 - 4mπ2)
        end
        cuhre(integrand, 5, 1)[1][1]
    end
    diagint(2) / diagint(1)
end
d2_mc = diag(BmatMC)[2]
d2_fu = diag(BmatFU)[2]

@show d2_analy, d2_mc, d2_fu
