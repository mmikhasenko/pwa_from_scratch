using Cuba

push!(LOAD_PATH, "src")
using DalitzPlotAnalysis

using amplitudes_compass

function the_big_integral(s, i::Int64)
  function integrand(x, f)
    s1min = 4.0 * mπ2
    s1max = (sqrt(s) - mπ)^2

    s1 = s1min + x[1] * (s1max - s1min)
    cosθ1 = -1.0 + 2.0 * x[2]
    ϕ1 = -π + 2.0 * π * x[3]
    cosθ23 = -1.0 + 2.0 * x[4]
    ϕ23 = -π + 2.0 * π * x[5]
    a = COMPASS_wave(i, s, s1, cosθ1, ϕ1, cosθ23, ϕ23)
    λλ = λ(s, s1, mπ2) * λ(s1, mπ2, mπ2)
    (λλ < 0) && error("λλ < 0: check masses")
    f[1] = abs2(a) * sqrt(λ(s, s1, mπ2)) * sqrt(λ(s1, mπ2, mπ2)) / (s * s1) * (s1max - s1min)
  end
  return cuhre(integrand, 5, 1)
end

# index 2 for ρπ S-wave
the_big_integral(1.1, 2)


function the_sml_integral(s, i::Int64)
  function integrand(x, f)
    s1min = 4.0 * mπ2
    s1max = (sqrt(s) - mπ)^2

    s1 = s1min + x[1] * (s1max - s1min)
    z = -1.0 + 2.0 * x[2]
    ampl2 = 0.0
    for lm = -1:1
      ampl2 += abs2(COMPASS_wave_short(i, lm, s, s1, z))
    end
    f[1] = ampl2 * sqrt(λ(s, s1, mπ2)) * sqrt(λ(s1, mπ2, mπ2)) / (s * s1) * (s1max - s1min)
  end
  return cuhre(integrand, 2, 1)
end

the_small_integral(1.1, 2)[1][1]

##########################
@time [the_big_integral(1.1, 2)[1][1] for i in 1:3]
@time [the_sml_integral(1.1, 2)[1][1] for i in 1:3]

@profile [the_big_integral(1.1, 2)[1][1] for i = 1:2]
Juno.profiler()

#########################
using Plots
theme(:juno)
let ev = linspace(0.5, 2.5, 10)
  big_cal = [the_big_integral(e^2, 2)[1][1] for e in ev]
  sht_cal = [the_sml_integral(e^2, 2)[1][1] for e in ev]
  plot(ev, big_cal, lab="the 5d-integral")
  plot!(ev, sht_cal, lab="the 2d-integral")
end

#########################
let ev = linspace(0.5, 2.5, 50)
  sht_cal = [the_sml_integral(e^2, 2)[1][1] for e in ev]
  plot(ev, sht_cal, lab="the 2d-integral")
end
