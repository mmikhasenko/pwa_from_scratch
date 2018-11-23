using Plots
theme(:default)
pyplot()

using DelimitedFiles

push!(LOAD_PATH,"src")
using DalitzPlotAnalysis
using amplitudes_compass
using PWAHelper
using SDMHelper
using FittingPWALikelihood

################################
mass_bin_name = "1540_1560" # "2320_2340"
wavelist = get_wavelist(joinpath(pwd(),"src/wavelist_formated.txt");
    path_to_thresholds=joinpath(pwd(),"src/thresholds_formated.txt"),
    M3pi=Meta.parse(mass_bin_name[1:4])/1000)
Nwaves = size(wavelist,1)
wavenames = get_wavenames(wavelist)
wavebasis = get_wavebasis(wavelist)
################################

function readnlines(f,n)
    lines = ""
    for i=1:n
        lines=lines*"\n"*readline(f)
    end
    return lines
end

readmatrix(f, nlines) = readdlm(IOBuffer(readnlines(f,nlines)))

my_old = open("data/functions_1540_1560_t1_rd.txt") do io
    readmatrix(io, 10)
end

srijan = open("data/srijan_calcs/functions_1540_1560_t1_rd.txt") do io
    readmatrix(io, 10)
end

sort(collect(zip(1:170,(my_old .- srijan)[1,:])), by=x->x[2])

my_dta = open("data/variables_1540_1560_t1_rd.txt") do io
    readmatrix(io, 10)
end

current = [wf(my_dta[1,2:end]..., my_dta[1,1]) for wf in wavebasis]

max(abs.((real.(current) - my_old[1,1:85]))...)

old_comt = [COMPASS_wave(i,my_dta[1,1:end]...) for i in 1:88]

max(abs.((real.(old_comt)[1:4] - my_old[1,1:4]))...)
max(abs.((real.(old_comt)[1:4] - srijan[1,1:4]))...)
max(abs.((real.(  my_old)[1,1:4] - srijan[1,1:4]))...)

my_dta[1:4]




srijan[1,1:4]




my_old[1,1:4]




real.(old_comt)[1:4]





COMPASS_wave(4,my_dta[1,1:end]...)#, my_dta[1,1]

wavenames[1:4]



my_dta[1,1]
my_dta[1,2]
λ(my_dta[1,1], mπ2, my_dta[1,2])
