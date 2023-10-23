module FittingPWALikelihood

using LinearAlgebra
using PWAHelper

export createLLHandGRAD, minimize

using NLopt
using ForwardDiff

function createLLHandGRAD(PsiDT, form, block_masks)

    Nd, Nwaves = size(PsiDT)
    Tmap = get_parameter_map(block_masks, Nwaves)
    Np = size(Tmap, 2)::Int

    pblocks = make_pblock_inds(block_masks)
    vblocks = [
        let v = fill(0.0, Np)
            v[bl] .= 1.0
            v
        end for bl in pblocks
    ]

    BM = (real.(contract_to_intensity(form, block_masks)))::Array{Float64,2}

    COHTS!(X) = cohts!(X, pblocks)
    COHSQ(X) = cohsq(X, pblocks)
    EXTND(Ψ) = extnd(Ψ, Tmap)
    EXTND!(X, Ψ) = extnd!(X, Ψ, Tmap)

    transposedExtendedPsiDT = Array{Complex{Float64}}(undef, Np, Nd)
    for e in 1:Nd
        EXTND!(@view(transposedExtendedPsiDT[:, e]), @view(PsiDT[e, :]))
    end

    # function LLH(pars)
    #     v = EXTND(@view(transposedPsiDT[:,1]))
    #     res = sum(log,
    #         let
    #             EXTND!(v,@view(transposedPsiDT[:,e]))
    #             v .*= pars
    #             COHSQ(v)
    #         end for e in 1:Nd)
    #     bilin = Nd * (pars' * BM * pars)
    #     return -res + bilin
    # end
    # function LLH(pars)
    #     res = sum(log,
    #         let
    #             v = transposedExtendedPsiDT[:,e]
    #             v .*= pars
    #             COHSQ(v)
    #         end for e in 1:Nd)
    #     bilin = Nd * (pars' * BM * pars)
    #     return -res + bilin
    # end
    function LLH(pars)
        res = sum(log,
            @views sum(abs2,
                transposedExtendedPsiDT[bl, e] ⋅ pars[bl]
                for bl in pblocks)
            for e in 1:Nd)
        bilin = Nd * (pars' * BM * pars)
        return -res + bilin
    end

    # function GETDV!(X, psi, pars)
    #     ExtΨ = Array{Complex{Float64}}(undef, length(pars))
    #     EXTND!(ExtΨ, psi)
    #     Y = ExtΨ.*pars
    #     COHTS!(Y)
    #     Y .*= conj.(ExtΨ)
    #     X .= real.(Y)
    #     return
    # end
    # function GETDV!(X, ExtΨ, pars)
    #     Y = ExtΨ.*pars
    #     COHTS!(Y)
    #     Y .*= conj.(ExtΨ)
    #     X .= real.(Y)
    # end
    function GETDV!(X, ExtΨ, Ypars)
        @inbounds Ypars .*= ExtΨ
        @views for bl in pblocks
            @inbounds Ypars[bl] .= sum(Ypars[bl])
        end
        Ypars .*= conj.(ExtΨ)
        X .= real.(Ypars)
    end
    # function GETDV!(X, ExtΨ, pars)
    #     @views for bl in pblocks
    #         X[bl] .= real.(
    #             conj.(ExtΨ[bl]) .*
    #             (pars[bl]' * ExtΨ[bl])
    #         )
    #     end
    # end

    function GETHS(ExtΨ, pars)
        return sum(
            let v = ExtΨ .* bl
                real.(2 * v * v')
            end for bl in vblocks)
    end

    function HESSIAN(pars)
        Np = size(pars, 1)
        hes = fill(0.0, Np, Np)
        v = Array{Float64}(undef, Np)
        Yp = Array{Complex{Float64}}(undef, Np)
        for e in 1:Nd
            Yp .= pars
            GETDV!(v, @view(transposedExtendedPsiDT[:, e]), Yp)
            deriv = 2 * v
            vale = pars' * v
            num = GETHS(@view(transposedExtendedPsiDT[:, e]), pars) * vale - deriv * deriv'
            hes .-= num ./ vale^2
        end
        hes .+= BM * (2Nd)
        return hes
    end

    function LLH_and_GRAD!(pars, grad)
        val = 0.0
        grad .= 0.0
        temp = Array{Float64}(undef, Np)
        Yp = Array{Complex{Float64}}(undef, Np)
        for e in 1:Nd
            Yp .= pars
            GETDV!(temp, @view(transposedExtendedPsiDT[:, e]), Yp)
            vale = (pars' * temp)::Float64
            temp .*= (1.0 / vale)
            grad .-= temp
            val -= log(vale)
        end
        grad .*= 2.0
        BB = BM * pars
        val += Nd * (pars' * BB)
        grad .+= BB .* (2Nd)
        return val
    end

    function GRAD(pars)
        grad = fill(0.0, Np)
        temp = Array{Float64}(undef, Np)
        Yp = Array{Complex{Float64}}(undef, Np)
        for e in 1:Nd
            Yp .= pars
            GETDV!(temp, @view(transposedExtendedPsiDT[:, e]), Yp)
            vale = (pars' * temp)::Float64
            temp .*= 1 / vale
            grad .-= temp
        end
        grad .*= 2.0
        BB = BM * pars
        grad .+= BB * (2Nd)
        return grad
    end

    return LLH, GRAD, LLH_and_GRAD!, HESSIAN
end


function minimize(minusLogLikelihood, andDerive!;
    algorithm::Symbol=:LD_LBFGS,
    verbose::Int=0,
    maxeval::Int=50000,
    # parsprecision::Float64 = 1e-4,
    llhtolerance::Float64=1e-4,
    starting_pars::Vector{Float64}=error("I need the starting parameters! Say starting_pars = [...]"))
    i_iter = 1
    function to_minimize(x::Vector, grad::Vector)
        if length(grad) > 0
            v = andDerive!(x, grad)
        else
            v = minusLogLikelihood(x)
        end
        i_iter += 1
        verbose == 1 && println(i_iter, ": ", v)
        verbose == 2 && @show v, x
        verbose == 3 && @show v, grad
        return v
    end
    opt = Opt(algorithm, length(starting_pars)) # try LD_LBFGS || LD_MMA || LD_SLSQP
    # xtol_rel!(opt,parsprecision)
    maxeval!(opt, maxeval)
    ftol_abs!(opt, llhtolerance)

    min_objective!(opt, to_minimize)

    (minf, pars, ret) = optimize(opt, starting_pars)#rand(size(TT,2)))
    println("got $minf at $pars after some iterations (returned $ret)")

    pars
end


end
