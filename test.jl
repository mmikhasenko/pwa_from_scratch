## definitions
using DalitzPlotAnalysis
using amplitudes_compass

using Cuba
using Plots

m3π = collect(0.5:0.01:2.5)

Φ(s, m1sq, m2sq) = (s > (sqrt(m1sq)+sqrt(m2sq))^2) ? sqrt(λ(s, m1sq, m2sq))/(8 * π * s) : 0
function everything(i::Int64)
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
    return cuhre(integrand, 5, 1)
  end

  preresult = total.(m3π.^2, i)
  result = []
  for index in 1:length(preresult)
        push!(result,preresult[index][1][1])
      end


  writedlm("output_new/h$(i).txt",[m3π result])
  return result
end
## calculating integrals
misha = []
for index in [2,3,4,5,6,7,9,11,16,18]
  push!(misha, everything(index))
end
## plotting
function CreatePlot(i::Int64)
  stephan = readdlm("/mnt/data/compass/2008/integrals_stefan/h$(i).tc.out")
  misha = readdlm("output_new/h$(i).txt")
  stephan[1:201,2] ./= (stephan[1:201,1])
  let p = sum(stephan[1:201,2])
          plot(stephan[1:201,1],stephan[1:201,2]./p,label="stephan")
  end
  let p = sum(misha[:,2])
          plot!(m3π, misha[:,2]./p,label="misha")
  end
end

## interference 
function interference(i::Int64, j::Int64)
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
    return cuhre(integrand2, 5, 2)
  end

  preresult = total2.(m3π.^2, i, j)
  resultre = []
  resultim = []
  for index in 1:length(preresult)
    push!(resultre,preresult[index][1][1])
    push!(resultim,preresult[index][1][2])
  end


  writedlm("output_new/h$(i)_$(j).txt",[m3π resultre resultim])
  return resultre, resultim
end

inter211 = interference(2,11)
## test stuff
CreatePlot(1)
result
