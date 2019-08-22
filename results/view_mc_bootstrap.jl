path_to_working_folder = "data"
mass_bin_name = "1540_1560"
tslice = "t1"

using DelimitedFiles
using Statistics

using Plots
pyplot()

bmc = readdlm(joinpath(path_to_working_folder,"llhfit_mc_bootstrap_1540_1560_t1.txt"))
hess = readdlm(joinpath(path_to_working_folder,"invhes_llhfit_1540_1560_t1_95122.txt"))
bpars = let
    path_to_attemts = joinpath(path_to_working_folder,"llh_attmpts_1540_1560_t1.txt")
    path_to_best = split(readline(path_to_attemts),"\t")[1] # first in sorted
    vcat(readdlm(path_to_best)...)
end
brd = readdlm(joinpath(path_to_working_folder,"llhfit_bootstrap_1540_1560_t1.txt"))

function plot_bootstrap_results(i)
    σ = sqrt(abs(hess[i,i]))
    m = bpars[i]
    range = extrema(vcat(bmc[i+1,:],brd[i+1,:]))
    plot(layout=grid(2,1,heights=(0.85,0.15)), link=:x, size=(500,500))
    histogram!(sp=1, [bmc[i+1,:] brd[i+1,:]], bins=LinRange(range...,30), lab=["mc bootstrap" "rd bootstrap"], xaxis=nothing, α=0.4, title="parameter #$(i)")
    #
    plot!(sp=2, quantile(bmc[i+1,:],[0.16,0.84]), [0.25, 0.25], l=(5,0.2), lab="mc quantile")
    plot!(sp=2, quantile(brd[i+1,:],[0.16,0.84]), [0.75, 0.75], l=(5,0.2), lab="rd quantile")
    #
    plot!(sp=2, [m], [0.5], xerr=[σ], ylim=(0,1), m=(5,:red), yaxis=false, lab="hessian")
end

plot_bootstrap_results(3)

for i=1:180
    plot_bootstrap_results(i)
    savefig(joinpath("plots","bootstrap","par$(i).pdf"))
end
