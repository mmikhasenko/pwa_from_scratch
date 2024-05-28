push!(LOAD_PATH, "src")
using DalitzPlotAnalysis
using Cuba
using Plots


const mπ = 0.14;
const mπ2 = mπ^2;

for v in [("I", "-"), ("II", "+")]
  @eval function $(Symbol("fρ_" * v[1]))(σ)
    mρ = 0.7755
    Γρ = 0.156
    return 1.0 / (mρ^2 - σ + $(Symbol(v[2]))(1.0im * mρ * Γρ))
  end
end
fρ_I(1.1)

function qtb(s)
  function integrand(x, f)
    s1min = 4.0 * mπ2
    s1max = (sqrt(s) - mπ)^2

    s1 = s1min + x[1] * (s1max - s1min)
    f[1] = fρ_I(s1) * fρ_II(s1) * sqrt(λ(s, s1, mπ2) * λ(s1, mπ2, mπ2)) / (s * s1) * (s1max - s1min)
  end
  return cuhre(integrand, 1, 1)[1][1]
end

let ev = linspace(0.5, 2.5, 100)
  cal = [qtb(e^2) for e in ev]
  plot(ev, cal)
end


function angular_pow2(s, s1, z, m1sq, m2sq, m3sq)
  θ23 = acos(z)
  θ3 = acos(cross_basis_cosθ3(s1, z, m1sq, m2sq, m3sq, s))
  θ12 = acos(cross_basis_cosθ12(s1, z, m1sq, m2sq, m3sq, s))
  res = Wignerd(1, 0, 0, θ3 + θ12 - θ23)
end


function interference(s)
  function integrand(x, f)
    s1min = 4.0 * mπ2
    s1max = (sqrt(s) - mπ)^2

    s1 = s1min + x[1] * (s1max - s1min)
    z = -1.0 + 2.0 * x[2]
    s3 = cross_basis_s3(s1, z, mπ2, mπ2, mπ2, s)
    ampl2 = fρ_I(s1) * fρ_II(s3) * angular_pow2(s, s1, z, mπ2, mπ2, mπ2)
    f[1], f[2] = reim(ampl2 * sqrt(λ(s, s1, mπ2) * λ(s1, mπ2, mπ2)) / (s * s1) * (s1max - s1min))
  end
  return complex(cuhre(integrand, 2, 2)[1]...)
end

@time interference(1.1)

@profile interference(1.1)

let ev = linspace(0.5, 2.5, 100)
  cal = [interference(e^2) for e in ev]
  plot!(ev, real.(cal))
end

Juno.Profile.clear()
Juno.profiler()
