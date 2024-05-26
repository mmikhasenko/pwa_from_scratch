# using QuadGK
# using Plots
# using GSL
using Cuba

push!(LOAD_PATH, ENV["HOME"] * "/Documents/amplitude_analysis/modules")

using isobars
using Masses
using DalitzPlotAnalysis

function phsp(s)
    function integrand(x, f)
        # complex path
        #  end points:
        s1_i = 4mπ2
        s1_f = (sqrt(s) - mπ)^2 - 1e-8 # otherwise, rounding error
        s1_m = real(s1_f) * one(s)
        #  variables:
        s1 = (x[1] < 0.5) ? s1_i + 2x[1] * (s1_m - s1_i) : s1_m + 2(x[1] - 0.5) * (s1_f - s1_m)
        #  jacobian
        res = one(1.0im)
        res *= (x[1] < 0.5) ? 2(s1_m - s1_i) : 2(s1_f - s1_m)
        # angular variables
        τ1 = (2 * x[2] - 1, 0.0, 2 * x[3] - 1, π * (2 * x[4] - 1))
        τ3 = fill(zero(s), 4)
        s3, τ3[1], τ3[2], τ3[3], τ3[4] = change_basis(s1, τ1[1], τ1[2], τ1[3], τ1[4], mπ2, mπ2, mπ2, s)
        [
            if !isfinite(τi)
                print("s = ", s, ":\n", s1, ", ", τ1, "\n ", s3, ", ", τ3, "\n\n")
            end for τi in τ3
        ]
        # the angular part of the rho-pi S-wave
        # For the rho pi S-wave, it reads:
        #    cosθ1*cosθ23-sinθ1*sinθ23*cosϕ23
        ang = [
            begin
                cI = τ[1]
                sI = sqrt(1.0 - cI^2)
                cV = τ[3]
                sV = sqrt(1.0 - cV^2)
                sqrt(3.0) * (cI * cV - sI * sV * cos(τ[4]))
            end for τ in [τ1, τ3]]
        # value
        amp = f1_I(s1) * ang[1] - f1_I(s3) * ang[2]
        ampc = f1_II(s1) * ang[1] - f1_II(s3) * ang[2]
        res *= amp * ampc / 2
        # intagra jacobian factors
        res *= sqrt(λ(s, s1, mπ2) * λ(s1, mπ2, mπ2)) / (s * s1)
        f[1], f[2] = reim(res)
    end
    result = cuhre(integrand, 4, 2)
    complex(result[1]...) / (2 * π) / (8 * π)^2
end

function phsp_sim(s)
    function integrand(x, f)
        # complex path
        #  end points:
        s1_i = 4mπ2
        s1_f = (sqrt(s) - mπ)^2 - 1e-8 # otherwise, rounding error
        s1_m = real(s1_f) * one(s)
        #  variables:
        s1 = (x[1] < 0.5) ? s1_i + 2x[1] * (s1_m - s1_i) : s1_m + 2(x[1] - 0.5) * (s1_f - s1_m)
        #  jacobian
        res = one(1.0im)
        res *= (x[1] < 0.5) ? 2(s1_m - s1_i) : 2(s1_f - s1_m)
        # value
        amp = f1_I(s1)
        ampc = f1_II(s1)
        res *= amp * ampc
        # intagra jacobian factors
        res *= sqrt(λ(s, s1, mπ2) * λ(s1, mπ2, mπ2)) / (s * s1)
        f[1], f[2] = reim(res)
    end
    result = cuhre(integrand, 1, 2)
    complex(result[1]...) / (2 * π) / (8 * π)^2
end

# phsp_sim(1.2+0.1im)

let fout_name = "/tmp/phsp_sym_50.txt", func = phsp, nVx = 50, nVy = 50
    sxv = linspace(0.8, 2.0, nVx)
    syv = linspace(-1.0, 0.1, nVy)[end:-1:1]
    close(open(fout_name, "w")) # to clean the file
    for (j, sy) in enumerate(syv), (i, sx) in enumerate(sxv)
        v = func(sx + 1im * sy)
        file = open(fout_name, "a")
        writedlm(file, [sx sy real(v) imag(v)])
        close(file)
        perc = 100(nVx * (j - 1) + i) / (nVx * nVy)
        println("\%$perc done: ($i/$nVx, $j/$nVy), y = $sy")
    end
end

# data = readdlm("/tmp/phsp.txt")
#
# matr = hcat([data[(50(i-1)+1):(50i),3]+1im*data[(50(i-1)+1):(50i),4] for i in 1:50]...)
#
# contour(real(matr))
# contour!(imag(matr))
#
# heatmap(real(matr))
# heatmap(imag(matr))
#
