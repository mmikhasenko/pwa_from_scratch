module PartialWavesFromScratch

include("isobars.jl")

include("DalitzPlotAnalysis.jl") # need LinearAlgebra, GSL

include("amplitudes_compass.jl") # depends on DalitzPlotAnalysis.jl

include("PWAHelper.jl") # needs amplitudes_compass, DelimitedFiles, NLopt, ForwardDiff

include("FittingPWALikelihood.jl") # needs PWAHelper, LinearAlgebra, NLopt, ForwardDiff

include("SDMHelper.jl") # needs PWAHelper, DelimitedFiles, LinearAlgebra

include("PlotHelper.jl") # needs Statistics, Plots

end # module PartialWavesFromScratch
