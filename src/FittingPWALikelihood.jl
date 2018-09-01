module FittingPWALikelihood

using PWAHelper

export createLLHandGRAD, minimize

using NLopt
using ForwardDiff

function createLLHandGRAD(PsiDT, form, block_masks)

    Tmap = get_parameter_map(block_masks)
    Np = size(Tmap,2)

    pblocks = make_pblock_masks(block_masks)

    BM = real.(contract_to_intensity(form,block_masks));

    Nd = size(PsiDT,1);

    COHTS(X) = cohts(X,pblocks)
    COHSQ(X) = cohsq(X,pblocks)
    EXTND(Ψ) = extnd(Ψ,Tmap)

    function LLH(pars)
        @inbounds res = sum(log(COHSQ(EXTND(PsiDT[e,:]).*pars)) for e in 1:Nd)
        - res + real(sum(pars[i]*BM[i,j]*pars[j] for i=1:Np, j=1:Np)) * Nd
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
        end for bl in pblocks)
    end

    function HESSIAN(pars)
        Np = size(pars,1)
        hes = fill(0.0, Np, Np)
        for e in 1:Nd
            v = GETDV(PsiDT[e,:],pars)
            deriv = 2*v
            vale = pars'*v
            num =  GETHS(PsiDT[e,:],pars)*vale - deriv*deriv'
            hes .-= num ./ vale^2
        end
        hes .+= BM* (2Nd);
        return hes;
    end

    function LLH_and_GRAD!(pars, grad)
        val = 0.0; grad .= 0.0
        for e in 1:Nd
            @inbounds v = GETDV(PsiDT[e,:],pars)
            vale = pars'*v
            grad .-= v / vale
            val -= log(vale);
        end
        grad .*= 2.0

        BB = BM*pars;
        val += pars'*BB * Nd;
        grad .+= BB* (2Nd);
        return val;
    end

    function GRAD(pars)
        grad = fill(0.0,length(pars))
        for e in 1:Nd
            @inbounds v = GETDV(PsiDT[e,:],pars)
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
                  verbose::Int=0, precisn::Float64 = 1e-6,
                  starting_pars::Vector{Float64} = error("I need the starting parameters!"))
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
    xtol_rel!(opt,precisn)
    maxeval!(opt,500000)

    min_objective!(opt, to_minimize)

    (minf,pars,ret) = optimize(opt, starting_pars)#rand(size(TT,2)))
    println("got $minf at $pars = [m, Γ] after some iterations (returned $ret)")

    pars
end


end
