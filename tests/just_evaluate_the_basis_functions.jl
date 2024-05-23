push!(LOAD_PATH, "src")
using DalitzPlotAnalysis
using amplitudes_compass

event = let s = 1.5, s1 = 0.8
    [s, s1, 0.3, 0.1π, 0.3, 0.1π]
end
# = [COMPASS_wave(i,event...) for e=1:size(mm,1), i=1:88];

cal = [COMPASS_wave(i, event...) for i in 1:88]
for i in 1:88
    println(cal[i])
end

@show Z(1, 0, true, 0, 1, event[3:end]...)

@show ClebschGordon((2 * [1, 1, 1, 0, 1, 1])...)
@show WignerDϵ(true, 1, 0, 0, event[3:4]..., 0.0)
@show Wignerd(1, 1, 0, event[5])
