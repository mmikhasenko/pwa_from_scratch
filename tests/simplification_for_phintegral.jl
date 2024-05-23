push!(LOAD_PATH, "src")
using DalitzPlotAnalysis

let pars = change_basis(0.3^2, 1.0, 0.0, 0.5, 0.0, 0.14^2, 0.14^2, 0.14^2, 1.5)
    names = ["s3", "cosθ3", "ϕ3", "cosθ12", "ϕ12"]
    println("----------------")
    for (i, v) in enumerate(names)
        println("$v = ", pars[i])
    end
end

using amplitudes_compass
#
COMPASS_wave(26, 1.5, 0.3^2, 1.0, 0.0, 0.3, 0.0)
COMPASS_wave_short(26, 0, 1.5, 0.3^2, 0.3)

@time for i = 1:100000
    COMPASS_wave(2, 1.5, 0.3^2, 1.0, 0.0, 0.3, 0.0)
end

@time for i = 1:100000
    COMPASS_wave_short(2, 0, 1.5, 0.3^2, 0.3)
end

@profile for i = 1:100000
    COMPASS_wave_short(2, 0, 1.5, 0.3^2, 2rand() - 1)
    # COMPASS_wave(2,1.5,0.3^2,1.0,0.0,0.3,0.0)
end

@profile for i = 1:100000
    COMPASS_wave(2, 1.5, 0.3^2, 2rand() - 1, π * (2rand() - 1), 2rand() - 1, π * (2rand() - 1))
    # COMPASS_wave(2,1.5,0.3^2,1.0,0.0,0.3,0.0)
end


Juno.profiler()
