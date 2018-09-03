module PlotHelper
using Plots

export plotBSTsample, saveBSTSplot, plotBSTsummary

function plotBSTsample(ind, SDMs, SDM, SDM_RD, SDM_RD_err)
    bts = [real(s[ind,ind]) for s in SDMs]
    stephist(bts, bins=20, lab="bootstrap", title="SDM[$(ind),$(ind)]")
    vline!([real(SDM[ind,ind])], l=(:red, 2), lab="main fit")
    vline!(real.([SDM_RD[ind,ind]] .+ SDM_RD_err[ind,ind] .* [-1,1]), l=(:orange, 1.5), lab="official fit")
    vline!(quantile(bts,[0.16,0.84]), lab="quantile")
    plot!()
end

function saveBSTSplot(fout, SDMs, SDM, SDM_RD, SDM_RD_err)
    Nw = size(SDM,2)
    onames = Array{String}(Nw)
    for i in 1:Nw
        plotBSTsample(i, SDMs, SDM, SDM_RD, SDM_RD_err)
        onames[i] = fout * "-tmp-$(i).pdf"
        savefig(onames[i])
    end
    run(`pdfunite $onames $fout`)
    run(`rm $onames`)
end

function plotBSTsummary(combres; tosort=false, toannotate=false)
    Nw = size(combres,1)
    mat = tosort ? sortrows(combres, by=x->x[2], rev=true) : combres
    bar(mat[:,5], fillrange=mat[:,6], c=:green, l=nothing, lab="bstrap quantiles",
        xlab="# wave", ylab="Wave Intensity, SDM[#,#]")
    tosort && plot!(xaxis=nothing)
    toannotate && annotate!([(i+0.4,mat[i,6]+50,text(Int64(mat[i,1]),4)) for i in 1:Nw])
    scatter!(mat[:,3], yerr=mat[:,4], xerr=[0.5 for i in 1:size(mat,1)],
        m=(1, stroke(0.0)), lab="official fit")
    scatter!(mat[:,2], m=(1.5, :red, :d, stroke(0.0)), lab="main fit")
end


end
