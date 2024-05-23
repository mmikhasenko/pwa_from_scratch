push!(LOAD_PATH, "src")

using amplitudes_compass
using DelimitedFiles
function getReIm(func)
    ev = LinRange(2mπ, 2.5, 500)
    cal = func.(ev .^ 2)
    return [ev real.(cal) imag.(cal)]
end

function writeEverything2txt(path_to_write="/tmp")
    for (func, foutname) in zip([fρ, ff2, fρ3, fσ, ff0_980, ff0_1500],
        ["rho", "f2", "rho3", "pipiS", "f0980", "f01500"])
        writedlm(joinpath(path_to_write, "isobar_" * foutname * ".txt"), getReIm(func))
    end
end
writeEverything2txt()

using Plots
pyplot()
function plotEverything()
    plot(layout=grid(1, 2), size=(500, 350))
    for f in [fρ, ff2, fρ3]
        mtrx = getReIm(f)
        plot!(mtrx[:, 1], sqrt.(mtrx[:, 2] .^ 2 + mtrx[:, 3] .^ 2), subplot=1, lab="")
    end
    for f in [fσ, ff0_980, ff0_1500]
        mtrx = getReIm(f)
        plot!(mtrx[:, 1], sqrt.(mtrx[:, 2] .^ 2 + mtrx[:, 3] .^ 2), subplot=2, lab="")
    end
    plot!()
end

plotEverything()
