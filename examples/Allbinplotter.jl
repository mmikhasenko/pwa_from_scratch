tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
const Nwave_extended = 88

################################################################################
using DelimitedFiles
push!(LOAD_PATH,"src")
using SDMHelper
using FittingPWALikelihood
using amplitudes_compass
using PWAHelper
using Random
using Plots
using LinearAlgebra
using DelimitedFiles
using Statistics

#global Diagonal_element = Array{Float64}(100)
global Diagonal_element = Array{Float64}(100, Nwave_extended)
global m_3pi = Array{Float64}(100)
global v_qnt = [Array{Complex{Float64}}(undef,88,500) for idx in 1:100]
#global bts = Array{Float64}(undef,88,500)
###############################################################################
global z = 1
for i in 0.5:0.02:2.48
    k= Int(i *1000)
    #Defining the mass bin
    global mass_bin_name = (string(string(k,pad=4),"_",string(k+20 , pad=4)))
    M3pi = Meta.parse(mass_bin_name[1:4])/1000

    for app in ["rd", "mc", "fu"]
        pwf = path_to_working_folder
        @eval $(Symbol("basisfunc_"*app)) = joinpath($pwf,"functions_$(mass_bin_name)_$(tslice)_"*$app*".txt")
    end
    # get number of lines
    nlines_rd = Meta.parse(split(read(`wc -l $(basisfunc_rd)`,String)," ")[1])
    nlines_mc = Meta.parse(split(read(`wc -l $(basisfunc_mc)`,String)," ")[1])
    nlines_fu = Meta.parse(split(read(`wc -l $(basisfunc_fu)`,String)," ")[1])
    normfact = nlines_fu/nlines_mc*nlines_rd

    # read matrix of integrals
    BmatFU = read_cmatrix("data/integrmat_$(mass_bin_name)_$(tslice)_fu.txt");

    # read wavelist throw waves below threshol
    wavelist = get_wavelist(joinpath(path_wavelist,"wavelist_formated.txt");
             path_to_thresholds=joinpath(path_wavelist,"thresholds_formated.txt"),
             M3pi=M3pi)

    const Nwaves = size(wavelist,1)
    # Blocks in wavelist
    const noϵ =  [1] # flat wave
    const posϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="+"]
    const negϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="-"]
    # Model description
    const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]

    #SDM Calculation
    SDM = []
    list_of_files = split(readline("data/llhfit/$(mass_bin_name)/llh_attmpts_$(mass_bin_name)_$(tslice).txt"),"\t")[1]#ARGS[3:end]
    path_and_filename = list_of_files
    minpars = readdlm(path_and_filename);
    SDM = normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks)

    threshold_mask = get_threshold_mask(joinpath(path_wavelist,"thresholds_formated.txt"), M3pi, Nwave_extended)
    SDM_enlarged=enlarge_with_zeros!(SDM,threshold_mask)
    for g in 1:88
        Diagonal_element[z,g]=SDM_enlarged[g,g]
    end

    bootstrap_file = "data/llhfit_bootstrap_$(mass_bin_name)_$(tslice).txt"
    BootstrapResults = readdlm(bootstrap_file)
    Nbstrap_attempts =vcat(size(BootstrapResults[1,:])...)[1]
    SDMs = [Array{Complex{Float64}}(undef, Nwave_extended,Nwave_extended) for idx in 1:500]
    @time for b in 1:500
        _br = BootstrapResults[2:end,b]
        SDMs[b] .= enlarge_with_zeros!(normfact*pars_to_SDM(_br, BmatFU, ModelBlocks),threshold_mask)
    end
    bts = Array{Float64}(undef,88,500)
    for g in 1:88
        bts[g,:] .= [real(s[g,g]) for s in SDMs]
    end
    v_qnt[z] .= bts
    m_3pi[z]= i .+ 0.015
    z = z+1
    @show z

end

wavelist = readdlm(joinpath("src","wavelist_formated.txt"))
Diagonal_element

#global err = Array{Float64}(100, 500)

for j in 2:88
        global err = Array{Float64}(100, 500)
        for i in 1:100
            err[i,:]=real.(v_qnt[i][j,:])
        end

        qtls = Array{Float64}(100, 2)
        for i in 1:100
            qtls[i,:] .= quantile(err[i,:],[0.16,0.84])
        end

        # Plotting the results
        x = m_3pi - 0.015
        dxh = (x[2]-x[1])/2  # 10 MeV
        new_x = [x[1]-dxh, (x+dxh)...]
        #for i in 1:100
            #plot(xlab="# wave", ylab="Wave Intensity, SDM[#,#]",size(800,500))
            plot(new_x, Diagonal_element[:,j], seriestype=:stepbins,fill_between=fill(0,size(Diagonal_element,1)),
                                    lab="$(wavelist[j,2])", lw=0.25, lc=:black,xticks = 0:0.2:10,size(800,500)
                                    ,legend = :best,legendfontsize = 5,m=(1.5,:black, stroke(0.0)))
            # bar!(new_x,err, lab="",fillrange=err[:,2],fillalpha=0.02, c=:orange ,l=nothing)
            bar!(x,qtls[:,1], lab="Bootstrap Quantiles",fillrange=qtls[:,2],fillalpha=0.8, c=:orange
                    ,l=nothing,xlab="M3pi (GeV)", ylab="Events/20MeV",size(800,500),title = "$(wavelist[j,2])")
            # bar!(new_x,err[:,1], lab="bstrap quantiles",fillrange=err[:,2],fillalpha=0.2, c=:orange ,l=nothing,size(800,500))
            #plot(new_x, Diagonal_element[:,2], seriestype=:stepbins,
                                    #            lab="Bootstrap errors", lw=0.5, lc=:black,xticks = 0:0.2:10,size(800,500)
                                    #            ,legend = :best,legendfontsize = 5,m=(1.5,:orange, :d, stroke(0.0)) )
        #end
        #bar(transpose(massbin),transpose(second_element),xticks = 0:0.2:10,size=(800,500),
        #            xlab="M_(3pi)",ylab ="Magnitude", title = "SDM[2,2] for all bins")
        savefig(joinpath("plots/sdm_results/final/allwaves/with_quantiles","sdm_[$(j),$(j)]_$(tslice)_with_quantile.pdf"))
end
