module SDMHelper

using PWAHelper

export SDM_to_pars, pars_to_SDM

function pars_to_SDM(pars, Bmat, block_masks)
    Tmap = get_parameter_map(block_masks)
    pblocks = make_pblock_masks(block_masks)
    Nw = size(Bmat,1)
    nonormSDM = fill(0.0im, Nw, Nw)
    for bl in pblocks
        sp = shrnk(pars.*bl,Tmap)
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
    Tmap = get_parameter_map(block_masks)
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

end
