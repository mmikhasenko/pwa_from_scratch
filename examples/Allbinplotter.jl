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

global Diagonal_element = Array{Float64}(100)
global m_3pi = Array{Float64}(100)
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

    Diagonal_element[z]=SDM_enlarged[2,2]
    m_3pi[z]= i + 0.015
    z = z+1
    @show z

end

# Plotting the results

x = m_3pi
dxh = (x[2]-x[1])/2  # 10 MeV
new_x = [x[1]-dxh, (x+dxh)...] #
plot(new_x, Diagonal_element, seriestype=:stepbins, fill_between=fill(0,size(transpose(second_element),1)),
            lab="1-(1++)0+rhopiS ", lw=0.5, lc=:black,xticks = 0:0.2:10)

#bar(transpose(massbin),transpose(second_element),xticks = 0:0.2:10,size=(800,500),
#            xlab="M_(3pi)",ylab ="Magnitude", title = "SDM[2,2] for all bins")
savefig(joinpath("plots/sdm_results/final","sdm_[2,2]_$(tslice)_test.png"))
