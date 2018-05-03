using Plots

function CreatePlot(i::Int64)
  stephan = readdlm("/mnt/data/compass/2008/integrals_stefan/h$(i).n.out")
  misha = readdlm("output_new/h$(i).txt")
  stephan[2:202,2] ./= (stephan[2:202,1])
  m3π = collect(0.5:0.01:2.5)
  let p = sum(stephan[2:202,2])
          plot(stephan[2:202,1],stephan[2:202,2]./p,label="stephan", title=stephan[1,1]*" "*stephan[1,2]*" "*stephan[1,3]*" "*stephan[1,4])
  end
  let p = sum(misha[:,2])
          plot!(m3π, misha[:,2]./p,label="misha")
  end
end

CreatePlot(5)
for i in [2,3,4,5,6,7,9,11,16,18]
  CreatePlot(i)
  savefig("misha_vs_stephan_new$(i).pdf")
end

function CreatePlotInterference()
  misha2 = readdlm("output_new/h2.txt")
  misha11 = readdlm("output_new/h11.txt")
  misha211 = readdlm("output_new/h2_11.txt")
  
  misha211[:,2:3] ./= sqrt.(misha2[:,2].*misha11[:,2])
  plot(misha211[1:200,1], -misha211[1:200,2], label="-re(misha)")
  plot!(misha211[1:200,1], misha211[1:200,3], label="im(misha)")
  
  florian = readdlm("diag.txt")
  florian1317 = readdlm("offdiag_13_17.txt")
  florian1317[:,2:3] = florian1317[:,2:3]./sqrt.(florian[14,2:201].*florian[18,2:201])
  plot!(florian1317[:,1], florian1317[:,2], label="re(florian)", title="Interference:  $(florian[14,1])...$(florian[18,1])")
  plot!(florian1317[:,1], florian1317[:,3], label="im(florian)")
  savefig("interference_new.pdf")
  plot!()
end

CreatePlotInterference()
