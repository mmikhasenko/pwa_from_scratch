using Plots

function CreatePlot(i::Int64, fileextra = "", filefolder = "output", output::Bool = true)
  """
  i: index of the wave,
  fileextra: additional text in the output-filename (default "")
  folder: name of the folder where the data files are located (default "output")
  output: boolean to decide if the plot should be saved in the outputfile (default true)
  """
  stephan = readdlm("/mnt/data/compass/2008/integrals_stefan/h$(i).n.out")
  misha = readdlm("$(filefolder)/h$(i).txt")
  stephan[2:202,2] ./= (stephan[2:202,1])
  m3π = collect(0.5:0.01:2.5)
  let p = sum(stephan[2:202,2])
    plot(stephan[2:202,1],stephan[2:202,2]./p,label="stephan", title=stephan[1,1]*" "*stephan[1,2]*" "*stephan[1,3]*" "*stephan[1,4])
  end
  let p = sum(misha[:,2])
    plot!(m3π, misha[:,2]./p,label="misha")
  end
  if output 
    savefig("misha_vs_stephan$(fileextra)$(i).pdf")
  end
  plot!()
end

for i in [2]
  CreatePlot(i, "", "output_new", false)
end

function CreatePlotInterference(i_misha::Int64, j_misha::Int64, i_florian::Int64, j_florian::Int64, fileextra = "", filefolder = "output", output::Bool = true)
  mishai = readdlm("$(filefolder)/h$(i_misha).txt")
  mishaj = readdlm("$(filefolder)/h$(j_misha).txt")
  mishaij = readdlm("$(filefolder)/h$(i_misha)_$(j_misha).txt")
  
  mishaij[:,2:3] ./= sqrt.(mishai[:,2].*mishaj[:,2])
  plot(mishaij[1:200,1], -mishaij[1:200,2], label="-re(misha)")
  plot!(mishaij[1:200,1], mishaij[1:200,3], label="im(misha)")
  
  florian = readdlm("diag.txt")
  florianij = readdlm("offdiag_$(i_florian)_$(j_florian).txt")
  florianij[:,2:3] = florianij[:,2:3]./sqrt.(florian[i_florian+1,2:201].*florian[j_florian+1,2:201])
  plot!(florianij[:,1], florianij[:,2], label="re(florian)", title="Interference:  $(florian[i_florian+1,1])...$(florian[i_florian+1,1])")
  plot!(florianij[:,1], florianij[:,3], label="im(florian)")
  if output
    savefig("interference_$(fileextra)$(i_misha)_$(j_misha).pdf")
  end
  plot!()
end

CreatePlotInterference(2,11,13,17,"test", "output_new")
