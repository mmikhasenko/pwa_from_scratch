## definitions
push!(LOAD_PATH,"src")
using DalitzPlotAnalysis
using amplitudes_compass

using Cuba
using Plots

m3π = collect(0.5:0.5:2.5)

Φ(s, m1sq, m2sq) = (s > (sqrt(m1sq)+sqrt(m2sq))^2) ? sqrt(λ(s, m1sq, m2sq))/(8 * π * s) : 0
function everything(i::Int64, folder = "output")
  function total(s, i::Int64)
    function integrand(x,f)
      s3min = 4.0 * mπ2
      s3max = (sqrt(s) - mπ)^2

      s3 =  s3min + x[1] * (s3max - s3min)
      cosθ1 = -1.0 + 2.0*x[2]
      ϕ1 = -π + 2.0*π*x[3]
      cosθ23 = -1.0 + 2.0*x[4]
      ϕ23 = -π + 2.0*π*x[5]
      a = COMPASS_wave(i, s, s3, cosθ1, ϕ1, cosθ23, ϕ23)
      f[1] = abs2(a) * Φ(s, s3, mπ2) * Φ(s3, mπ2, mπ2)*(s3max-s3min)
    end
    return cuhre(integrand, 5, 1)[1][1]
  end

  result = Array{Float64}(length(m3π));
  @progress for (p,msq) in enumerate(m3π.^2)
    result[p] = total(msq, i)
  end

  writedlm("$(folder)/h$(i).txt",[m3π result])
  return result
end
## calculating integrals


for index in [2,11]
  everything(index, "output_new")
end

# let s = 0.257672, τ1 = [0.100016 -0.997405 -0.694468 0.0442564 -1.24848],
#   τ3 = [0.109153 0.529428 -2.63108 0.14344 -π+0.0805736]
#   τ3p = change_basis(τ1...,mπ2,mπ2,mπ2,s)
#   # τ2 = change_basis(τ3...,mπ2,mπ2,mπ2,s)
#   # τ1p = change_basis(τ2...,mπ2,mπ2,mπ2,s)
#   [τ3[i]-τ3p[i] for i in 1:5]
# end

## interference
function interference(i::Int64, j::Int64, folder = "output")
  function total2(s, i::Int64, j::Int64)
    function integrand2(x,f)
      s3min = 4.0 * mπ2
      s3max = (sqrt(s) - mπ)^2

      s3 =  s3min + x[1] * (s3max - s3min)
      cosθ1 = -1.0 + 2.0*x[2]
      ϕ1 = -π + 2.0*π*x[3]
      cosθ23 = -1.0 + 2.0*x[4]
      ϕ23 = -π + 2.0*π*x[5]
      a = COMPASS_wave(i, s, s3, cosθ1, ϕ1, cosθ23, ϕ23)
      b = COMPASS_wave(j, s, s3, cosθ1, ϕ1, cosθ23, ϕ23)
      f[1], f[2] = reim(a * conj(b) * Φ(s, s3, mπ2) * Φ(s3, mπ2, mπ2)*(s3max-s3min))
    end
    return complex(cuhre(integrand2, 5, 2)[1]...)
  end

  preresult = Array{Complex{Float64}}(length(m3π));
  @progress for (p,msq) in enumerate(m3π.^2)
    preresult[p] = total2(msq, i, j)
  end
  resultre = real.(preresult)
  resultim = imag.(preresult)

  writedlm("$(folder)/h$(i)_$(j).txt",[m3π resultre resultim])
  return resultre, resultim
end

interference(2,11, "output_new")
