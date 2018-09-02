module PWAHelper

using amplitudes_compass

# using JLD
using Juno

export precalculate_compass_basis, read_precalc_basis
# export precalculate_compass_basis_txt, read_precalc_basis_txt
export get_parameter_map, make_pblock_masks
export extnd, shrnk, cohsq, cohts
export contract_to_intensity


# function precalculate_compass_basis(fin,fout)
#     mm = readdlm(fin); Nd = size(mm,1)
#     m2 = fill(0.0im, Nd, 88)
#     @progress for e in 1:Nd
#         for i in 1:88
#             m2[e,i] = COMPASS_wave(i,mm[e,:]...);
#         end
#     end
#     save(fout,"real",real(m2),"imag",imag(m2))
# end
# function read_precalc_basis(fname)
#     ld = load(fname)
#     ld["real"]+1im*ld["imag"]
# end

function precalculate_compass_basis(fin,fout)
    mm = readdlm(fin); Nd = size(mm,1)
    m2 = fill(0.0im, Nd, 88)
    @progress for e in 1:Nd
        for i in 1:88
            m2[e,i] = COMPASS_wave(i,mm[e,:]...);
        end
    end
    writedlm(fout,[real(m2) imag(m2)])
end
function read_precalc_basis(fname)
    ld = readdlm(fname)
    Nh = div(size(ld,2),2)
    ld[:,1:Nh]+1im*ld[:,(Nh+1):end]
end


function get_parameter_map(block_masks)
    temp = []; numb = []
    for (i,bl) in enumerate(block_masks)
        # push false for the first parameter, true for others
        push!(temp,false,[true for i=1:(sum(bl)-1)]...);
        push!(numb,collect(1:88)[bl]...)
    end
    Tmap = fill(0,2,sum(temp+1))
    count=1
    for (i,b) in enumerate(temp)
        Tmap[1,count] = numb[i]
        count+=1
        if b
            Tmap[2,count] = numb[i];
            count+=1
        end
    end
    Tmap
end
# get_parameter_map([noϵ,posϵ,negϵ,negϵ])

function make_pblock_masks(block_masks)
    just_true = [fill(true,2sum(bl)-1) for bl in block_masks]
    tmp = fill(false,sum(sum(jt) for jt in just_true))
    # to be filed
    pbls = copy(block_masks)
    counter = 1
    for (i,jt) in enumerate(just_true)
        pbls[i] = copy(tmp)
        jtlen = length(jt)
        pbls[i][counter:(counter+jtlen-1)] .= true
        counter += jtlen
    end
    pbls
end


extnd(Ψ, tmap) = [((tmap[2,i]==0) ? Ψ[tmap[1,i]] : 1im*Ψ[tmap[2,i]]) for i in 1:size(tmap,2)]
function shrnk(p, tmap)
    outv = fill(0.0im,88)
    for i in 1:size(tmap,2)
        (tmap[1,i] != 0) && (outv[tmap[1,i]] += p[i]);
        (tmap[2,i] != 0) && (outv[tmap[2,i]] += 1.0im*p[i]);
    end
    outv
end

function cohsq(X, pblocks)
    sum(abs2(sum(X[bl])) for bl in pblocks)
end
function cohts(X, pblocks)
    sum(sum(X[bl])*bl for bl in pblocks)
end

function contract_to_intensity(ΨΨstar, block_masks)
    Tmap = get_parameter_map(block_masks)
    Np = size(Tmap,2)
    pblocks = make_pblock_masks(block_masks)
    COH = sum(bl*bl' for bl in pblocks)
    get(v,w) = (v==0 || w==0) ? 0.0im : ΨΨstar[v,w]
    [COH[i,j]*( # meaning is ([1]-i[2]) * ([1]+i[2])
            get(Tmap[1,i],Tmap[1,j]) +
            get(Tmap[2,i],Tmap[2,j]) -
        1im*get(Tmap[2,i],Tmap[1,j]) +
        1im*get(Tmap[1,i],Tmap[2,j])
            ) for j in 1:Np, i in 1:Np]
end

end
