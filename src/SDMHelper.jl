module SDMHelper

using DelimitedFiles
using PWAHelper

using LinearAlgebra

export SDM_to_pars, pars_to_SDM
export write_cmatrix, read_cmatrix
export constract_values
export enlarge_with_zeros!

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

function enlarge_with_zeros!(SDM_new,SDM_old,thresholds_mask)
    size(SDM_new,1) != length(thresholds_mask) && error("Check size of SDM_new")
    size(SDM_old,1) != sum(thresholds_mask) && error("Check size of SDM_old")
    SDM_new[thresholds_mask,thresholds_mask] .= SDM_old
    return
end

function enlarge_with_zeros!(SDM_old,thresholds_mask)
    size(SDM_old,1) != sum(thresholds_mask) && error("Check size of SDM_old")
    SDM_new = fill(0.0im, length(thresholds_mask), length(thresholds_mask))
    SDM_new[thresholds_mask,thresholds_mask] .= SDM_old
    return SDM_new
end

end
