module FittingPWALikelihood

using LinearAlgebra
using PWAHelper

export createLLHandGRAD, minimize

using NLopt
using ForwardDiff

function createLLHandGRAD(PsiDT, form, block_masks)

    Nd, Nwaves = size(PsiDT);
    transposedPsiDT = transpose(PsiDT)

    Tmap = get_parameter_map(block_masks,Nwaves)

    Np = size(Tmap,2)::Int

    pblocks = make_pblock_inds(block_masks)
    vblocks = [let v = fill(0.0,Np)
            v[bl] .= 1.0
            v
        end for bl in pblocks]

    BM = (real.(contract_to_intensity(form,block_masks)))::Array{Float64,2};

    COHTS(X) = cohts(X,pblocks)
    COHSQ(X) = cohsq(X,pblocks)
    EXTND(Ψ) = extnd(Ψ,Tmap)

    function LLH(pars)
        @inbounds res = sum(log(COHSQ(EXTND(@view(transposedPsiDT[:,e])).*pars)) for e in 1:Nd)
        bilin = Nd * (pars' * BM * pars)
        return -res + bilin
    end

    function GETDV(psi, pars)
        ExtΨ = EXTND(psi)
        cExtΨ = conj.(ExtΨ)
        cv = cExtΨ.*COHTS(ExtΨ.*pars)
        return real.(cv)
    end

    function GETHS(psi, pars)
        ExtΨ = EXTND(psi)
        return sum(let v = ExtΨ.*bl
            real.(2*v*v')
        end for bl in vblocks)
    end

    function HESSIAN(pars)
        Np = size(pars,1)
        hes = fill(0.0, Np, Np)
        for e in 1:Nd
            v = GETDV(@view(transposedPsiDT[:,e]),pars)
            deriv = 2*v
            vale = pars'*v
            num =  GETHS(@view(transposedPsiDT[:,e]),pars)*vale - deriv*deriv'
            hes .-= num ./ vale^2
        end
        hes .+= BM* (2Nd);
        return hes;
    end

    function LLH_and_GRAD!(pars, grad)
        val = 0.0; grad .= 0.0
        for e in 1:Nd
            @inbounds v = GETDV(@view(transposedPsiDT[:,e]),pars)
            vale = pars'*v
            grad .-= v / vale
            val -= log(vale);
        end
        grad .*= 2.0
        BB = BM*pars;
        val += Nd * (pars ⋅ BB);
        grad .+= BB * (2Nd);
        return val;
    end

    function GRAD(pars)
        grad = fill(0.0,length(pars))
        for e in 1:Nd
            @inbounds v = GETDV(@view(transposedPsiDT[:,e]),pars)
            vale = pars'*v
            grad .-= v / vale
        end
        grad .*= 2.0
        BB = BM*pars;
        grad .+= BB* (2Nd);
        return grad;
    end

    return LLH, GRAD, LLH_and_GRAD!, HESSIAN
end


function minimize(minusLogLikelihood, andDerive!;
                  algorithm::Symbol = :LD_LBFGS,
                  verbose::Int = 0,
                  maxeval::Int = 50000,
                  parsprecision::Float64 = 1e-4,
                  llhtolerance::Float64 = 1e-6,
                  starting_pars::Vector{Float64} = error("I need the starting parameters! Say starting_pars = [...]"))
    function to_minimize(x::Vector, grad::Vector)
        if length(grad) > 0
            v = andDerive!(x,grad)
        else
            v = minusLogLikelihood(x)
        end
        verbose==1 && @show v
        verbose==2 && @show v,x
        verbose==3 && @show v,grad
        return v;
    end
    opt = Opt(algorithm, length(starting_pars)) # try LD_LBFGS || LD_MMA || LD_SLSQP
    xtol_rel!(opt,parsprecision)
    maxeval!(opt,maxeval)
    ftol_abs!(opt,llhtolerance);

    min_objective!(opt, to_minimize)

    (minf,pars,ret) = optimize(opt, starting_pars)#rand(size(TT,2)))
    println("got $minf at $pars after some iterations (returned $ret)")

    pars
end


end
