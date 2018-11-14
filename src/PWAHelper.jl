module PWAHelper

using amplitudes_compass

using DelimitedFiles
using Juno

export precalculate_compass_basis, read_precalc_basis
# export precalculate_compass_basis_txt, read_precalc_basis_txt
export get_npars
export get_parameter_map, make_pblock_inds
export extnd, extnd!, shrnk, cohsq, cohts, cohts!
export contract_to_intensity, get_intesity
export get_parameter_ranges, normalize_pars!

function precalculate_compass_basis(basis,fin,fout)
    mm = readdlm(fin); Nd = size(mm,1)
    m2 = fill(0.0im, Nd, length(basis))
    @progress for e in 1:Nd
        for (i,b) in enumerate(basis)
            m2[e,i] = b(mm[e,2:end]..., mm[e,1]);
        end
    end
    writedlm(fout,[real(m2) imag(m2)])
end
function read_precalc_basis(fname)
    ld = readdlm(fname)
    Nh = div(size(ld,2),2)
    return ld[:,1:Nh]+1im*ld[:,(Nh+1):end]
end

function get_parameter_map(block_inds, Nw)
    temp = []; numb = []
    for (i,bl) in enumerate(block_inds)
        # push false for the first parameter, true for others
        push!(temp,false,[true for i=1:(length(bl)-1)]...);
        push!(numb,bl...)
    end
    Tmap = fill(0,2,sum(temp .+ 1))
    count=1
    for (i,b) in enumerate(temp)
        Tmap[1,count] = numb[i]
        count+=1
        if b
            Tmap[2,count] = numb[i];
            count+=1
        end
    end
    return Tmap
end

function get_npars(block_inds)
    sum(2*length(bl)-1 for bl in block_inds)
end

function make_pblock_inds(block_inds)
    Nps = [2*length(bl)-1 for bl in block_inds]
    counter = 1
    pbls = [let v = collect(counter:(counter+Np-1))
                counter += Np
                v
            end for Np in Nps]
    return pbls
end

extnd(Ψ, tmap)  = [((tmap[2,i]==0) ? Ψ[tmap[1,i]] : 1im*Ψ[tmap[2,i]]) for i in 1:size(tmap,2)]
function extnd!(X, Ψ, tmap)
    for i in 1:size(tmap,2)
        @inbounds X[i] = (tmap[2,i]==0) ? Ψ[tmap[1,i]] : 1im*Ψ[tmap[2,i]]
    end
end
function shrnk(p, tmap)
    Nw = tmap[2,end] # work around
    outv = fill(0.0im,Nw)
    for i in 1:size(tmap,2)
        (tmap[1,i] != 0) && (outv[tmap[1,i]] += p[i]);
        (tmap[2,i] != 0) && (outv[tmap[2,i]] += 1.0im*p[i]);
    end
    outv
end

function cohsq(X, pblocks)
    sum(abs2, (sum(@view(X[bl]))) for bl in pblocks)
end
function cohts(X, pblocks)
    Np = length(X)
    vblocks = [let v = fill(0.0,Np)
        v[bl] .= 1.0
        v
    end for bl in pblocks]
    sum(sum(X[bl])*pl for (pl,bl) in zip(vblocks, pblocks))
end
function cohts!(X, pblocks)
    # Np = length(X)
    # vblocks = [
    #     let v = fill(0.0,Np)
    #         v[bl] .= 1.0
    #         v
    #     end for bl in pblocks]
    # sum(sum(X[bl])*pl for (pl,bl) in zip(vblocks, pblocks))
    for bl in pblocks
        @inbounds X[bl] .= sum(@view(X[bl]))
    end
end

function contract_to_intensity(ΨΨstar, block_inds)
    Tmap = get_parameter_map(block_inds, size(ΨΨstar,1))
    Np = size(Tmap,2)
    pblocks = make_pblock_inds(block_inds)
    COH = fill(0.0, Np, Np)
    for bl in pblocks
        COH[bl,bl] .= 1.0
    end
    get(v,w) = (v==0 || w==0) ? 0.0im : ΨΨstar[v,w]
    [COH[i,j]*( # meaning is ([1]-i[2]) * ([1]+i[2])
            get(Tmap[1,i],Tmap[1,j]) +
            get(Tmap[2,i],Tmap[2,j]) -
        1im*get(Tmap[2,i],Tmap[1,j]) +
        1im*get(Tmap[1,i],Tmap[2,j])
            ) for j in 1:Np, i in 1:Np]
end


function get_intesity(pars, form, block_inds)
    pblocks = make_pblock_inds(block_inds)
    BM = real.(contract_to_intensity(form,block_inds))
    Np = length(pars)
    sum(pars[i]*BM[i,j]*pars[j] for i=1:Np, j=1:Np)
end

function get_parameter_ranges(form, block_inds)
    Np = get_npars(block_inds)
    [let p = fill(0.0, Np)
        p[i] = 1.0
        1/sqrt(get_intesity(p, form, block_inds))
     end for i in 1:Np]
end

function normalize_pars!(pars, form, block_inds)
    Intens = get_intesity(pars, form, block_inds)
    pars ./=  sqrt(Intens)
    pars
end


end
