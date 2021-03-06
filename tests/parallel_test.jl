# addprocs(3)

push!(LOAD_PATH,"src")
using DalitzPlotAnalysis

v = @parallel for i in 1:50
    sleep(1)
    val = DalitzPlotAnalysis.λ(1,3,4.0)
    writedlm("/tmp/test-$(i).txt", fill(val, 20))
end

fetch(v[1])

BootstrapResults = let Nb = 500
    # path to save
    path_to_tmp_res = "/tmp"
    # path to save
    minpars0 = vcat(readdlm("minpars_compass_$(mass_bin_name).txt")...);
    res = Matrix{Float64}(Nb,length(minpars0))
    llh = Vector{Float64}(Nb)
    # b = 16
    @parallel for b in 1:Nb
        @show b, "progress is ", b/Nb
        writedlm(joinpath(path_to_tmp_res,"BootstrapResults-$(b).txt"), res[b,:])
        writedlm(joinpath(path_to_tmp_res,"llh-$(b).txt"), llh[b])
    end
    res
end
