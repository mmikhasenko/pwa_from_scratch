module PlotHelper
using Plots
using Statistics

export plotBSTsample, saveBSTSplot, plotBSTsummary, plotPolarSDM
export format_wavename

function plotBSTsample(ind::Int, SDMs, SDM, SDM_RD, SDM_RD_err)
    bts = [real(s[ind,ind]) for s in SDMs]
    stephist(bts, bins=20, lab="bootstrap", title="SDM[$(ind),$(ind)]")
    vline!([real(SDM[ind,ind])], l=(:red, 2), lab="main fit")
    vline!(real.([SDM_RD[ind,ind]] .+ SDM_RD_err[ind,ind] .* [-1,1]), l=(:orange, 1.5), lab="official fit")
    vline!(quantile(bts,[0.16,0.84]), lab="quantile")
    plot!()
end

function plotBSTsample(func, SDMs, SDM, SDM_RD, SDM_RD_err)
    bts = [func(s) for s in SDMs]
    stephist(bts, bins=20, lab="bootstrap", title="")
    vline!([func(SDM)], l=(:red, 2), lab="main fit")
    vline!([func(SDM_RD)], l=(:orange, 1.5), lab="official fit")
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


function construct_values(ind, SDMs, SDM, SDM_RD, SDM_RD_err)
    bts = [real(s[ind,ind]) for s in SDMs]
    v_main = real(SDM[ind,ind])
    v_off = real(SDM_RD[ind,ind])
    v_off_err = real(SDM_RD_err[ind,ind])
    v_qnt = quantile(bts,[0.16,0.84]) # sigma to every direction
    [ind v_main v_off v_off_err v_qnt...]
end

function plotBSTsummary(SDMs, SDM, SDM_RD, SDM_RD_err; tosort=false, toannotate=false)
    combres = vcat([construct_values(i, SDMs, SDM, SDM_RD, SDM_RD_err) for i in 1:size(SDM,1)]...);
    plotBSTsummary(combres; tosort=tosort, toannotate=toannotate)
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


function format_wavename(name)
    name == "flat" && return "\$\\mathrm{FLAT}\$"
    tmps = "\$$(name[4])^{$(name[5:6])}$(name[8])^{$(name[9])}\\,$(name[10:end-1])\\,$(name[end])\$"
    for (k,v) in Dict("rho" => "\\rho",
                       "f2" => "f_2",
                       "pi" => "\\pi",
                       "(\\pi\\pi)_S" => "\\sigma")
        tmps = replace(tmps,k,v)
    end
    tmps
end

function plotPolarSDM(SDM; cutoff_scale=0.4)
    plot(proj=:polar)
    arg(z) = atan2(imag(z), real(z))
    SDM_with_sign = sqrt.(diag(SDM)).*cis.([arg(SDM[2,i]) for i in 1:size(SDM,1)]);
    cutoff = cutoff_scale*max(abs.(SDM_with_sign)...)
    for (i,v) in enumerate(SDM_with_sign)
        plot!([0, arg(v)], [0, abs(v)], lab="", l=(1))
        (abs(v) > cutoff) && annotate!(
            [(arg(v), 1.1*abs(v), text(format_wavename(wavenames[i]),10))])
    end
    plot!()
end


end
