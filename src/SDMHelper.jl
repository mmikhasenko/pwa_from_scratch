module SDMHelper

using DelimitedFiles
using PWAHelper

export SDM_to_pars, pars_to_SDM
export write_cmatrix, read_cmatrix
export constract_values

function pars_to_SDM(pars, Bmat, block_masks)
    Tmap = get_parameter_map(block_masks, size(Bmat,1))
    pblocks = make_pblock_inds(block_masks)
    Nw = size(Bmat,1)
    nonormSDM = fill(0.0im, Nw, Nw)
    for bl in pblocks
        pars_bl = fill(0.0, length(pars))
        pars_bl[bl] .= pars[bl]
        sp = shrnk(pars_bl,Tmap)
        nonormSDM .+= [conj(sp[i])*sp[j] for i in 1:Nw, j in 1:Nw]
    end
    diagSHBM = real.(diag(Bmat))
    lSDM = [nonormSDM[i,j]*sqrt(diagSHBM[i]*diagSHBM[j]) for i in 1:Nw, j in 1:Nw]
    return lSDM
end

function SDM_to_pars(SDM, Bmat, block_masks)
    pars = sqrt.(diag(SDM))
    for i in 2:size(pars,1)
        el = SDM[i-1,i]
        (el != 0.0) && (pars[i] *= el/abs(el) * pars[i-1]/abs(pars[i-1]))
    end
    pars = pars ./ sqrt.(diag(Bmat))
    # convert to re im
    Tmap = get_parameter_map(block_masks, size(Bmat,1))
    fpars = fill(0.0, size(Tmap,2))
    for i in 1:size(Tmap,2)
        if Tmap[1,i] != 0
            (fpars[i] == 0.0) && (fpars[i] = real(pars[Tmap[1,i]]))
        else
            (fpars[i] == 0.0) && (fpars[i] = imag(pars[Tmap[2,i]]))
        end
    end
    return fpars
end

function write_cmatrix(sdm, fout)
    writedlm(fout, [real(sdm) imag(sdm)])
end

function read_cmatrix(fin)
    ld = readdlm(fin)
    Nh = div(size(ld,2),2)
    ld[:,1:Nh] + 1im .* ld[:,(Nh+1):end]
end


function constract_values(ind, SDMs, SDM, SDM_RD, SDM_RD_err)
    bts = [real(s[ind,ind]) for s in SDMs]
    v_main = real(SDM[ind,ind])
    v_off = real(SDM_RD[ind,ind])
    v_off_err = real(SDM_RD_err[ind,ind])
    v_qnt = quantile(bts,[0.16,0.84]) # sigma to every direction
    [ind v_main v_off v_off_err v_qnt...]
end


end
