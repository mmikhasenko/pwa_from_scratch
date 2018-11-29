#mass_bin_name_list =[,"1220_1240",] #ARGS[1] #"2320_2340"
tslice = "t1"
path_wavelist = "src"
path_to_working_folder = "data"
Nbstrap_attempts = 50 # No. of bootstrap attempts
SDM_mass_bin_name =string(trunc(Int,100- ((2480 - Meta.parse(split(mass_bin_name,"_")[1])) ./ 20)))                #ARGS[2]

using DelimitedFiles
push!(LOAD_PATH,"src")
using SDMHelper
using FittingPWALikelihood
using amplitudes_compass
using PWAHelper
using Random
using Plots
using LinearAlgebra

global second_element = Array{Float64}(1,100)
global starting = 0.515
global massbin = Array{Float64}(1,100)
###############################################################################
global z = 1
for i in 5:24

    for j in ["00_20","20_40","40_60","60_80","80_"]

            k=split(j,"_")[1]
            l=split(j,"_")[2]
            if i <9
                if j == "80_"
                    m=i+1
                    global mass_bin_name = string("0$(i)$(k)_0$(m)00")

                else
                    global mass_bin_name = string("0$(i)$(k)_0$(i)$(l)")
                end
            elseif i == 9
                if j == "80_"
                    m=i+1
                    global mass_bin_name = string("0$(i)$(k)_$(m)00")

                else
                    global mass_bin_name = string("0$(i)$(k)_0$(i)$(l)")
                end
            else
                if j == "80_"
                    m=i+1
                    global mass_bin_name = string("$(i)$(k)_$(m)00")
                else
                    global mass_bin_name = string("$(i)$(k)_$(i)$(l)")
                end
            end



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
                     M3pi=Meta.parse(mass_bin_name[1:4])/1000)

            const Nwaves = size(wavelist,1)

            # Blocks in wavelist
            const noϵ =  [1] # flat wave
            const posϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="+"]
            const negϵ = [i for (i,ϵ) in enumerate(wavelist[:,6]) if ϵ=="-"]

            # Model description

            const ModelBlocks = [noϵ, posϵ, negϵ, negϵ]
            list_of_files = split(readline("data/llhfit/$(mass_bin_name)/llh_attmpts_$(mass_bin_name)_$(tslice).txt"),"\t")[1]
            println("Start calculating SDMs")
            #Loop to calculate SDMs
            SDM = []
            path_and_filename = list_of_files
            #for path_and_filename in list_of_files
                minpars = readdlm(path_and_filename);
                # calculate
                SDM = normfact*pars_to_SDM(minpars, BmatFU, ModelBlocks)
            second_element[1,z]=SDM[17,17]
            massbin[1,z]= starting
            starting += 0.02
            z = z+1
            @show z
    end
end

second_element
massbin
x = transpose(massbin)
dxh = (x[2]-x[1])/2  # 10 MeV
new_x = [x[1]-dxh, (x+dxh)...] #
plot(new_x, transpose(second_element), seriestype=:stepbins, fill_between=fill(0,size(transpose(second_element),1)),
            lab="1++0+ RhoPi S", lw=0.5, lc=:black)

#bar(transpose(massbin),transpose(second_element),xticks = 0:0.2:10,size=(800,500),
#            xlab="M_(3pi)",ylab ="Magnitude", title = "SDM[2,2] for all bins")
savefig(joinpath("plots/sdm_results/final","sdm_[17,17]_$(tslice)_test.png"))
