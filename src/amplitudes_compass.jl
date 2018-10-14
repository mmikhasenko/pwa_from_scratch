module amplitudes_compass
# using masses: mπ, mπ2
using QuadGK  # needed for waves notmalization
using DelimitedFiles
using DalitzPlotAnalysis: λ, change_basis, Z, WignerDϵ, WignerD, Wignerd, ClebschGordon

export swap_kin_parameters
export fρ, ff2, fρ3, fσ, ff0_980, ff0_1500, BlttWskpf
export compass_jmels_basis_psi
# export COMPASS_wave
#, COMPASS_waves, COMPASS_wave_short
# export wavenames, wavesfile, Nwaves, thresholds_filter
export get_wavelist, get_wavenames, get_wavebasis

# export buildbasis

const mπ = 0.13956755; const mπ2 = mπ^2;
const mK = 0.493677; const mK2 = mK^2;
const mK0 = 0.497614; const mK02 = mK0^2;

export mπ, mπ2

const BlttWskpf = [z->z/(1+z),
                   z->z^2/(9+3z+z^2),
                   z->z*z*z/(z*z*z+6z*z+45z+225),
                   z->z*z*z*z/(z*z*z*z+10z*z*z+135z*z+1575z+11025),
                   z->z*z*z*z*z/(z*z*z*z*z+15z*z*z*z+315z*z*z+6300z*z+99225z+893025),
                   z->z*z*z*z*z*z/(z*z*z*z*z*z+21z*z*z*z*z+630z*z*z*z+18900z*z*z+496125z*z+9823275z+108056025)];

################### SIGMA #################
# AMP Table 1, M solution: f_2^2
# AMP Table 1, M solution: f_1^1 and f_2^1
# Last item is Katchaev modification
const _a = [0.1131, 0.0]

# AMP Table 1, M solution: c_11^0
# last item - Katchaev modification
const _c = [0.0337, -0.3185, -0.0942, -0.5927, 0.0]
const _sP = [-0.0074, 0.9828]

function _fσ(s::Number)
    (s < 4mπ2) && return 0.0;

    rho00 = sqrt(λ(s, mπ2, mπ2))/s;
    scale = (s / (4 * ((mK+mK0)/2.)^2)) - 1.;

    M00 = 0.0;
    for i = 1:length(_sP)
        M00 += _a[i] / (s - _sP[i]);
    end
    for i = 1:length(_c)
        M00 += scale^(i-1) *_c[i];
    end

  return 1.0/(M00 - 1im*rho00);
end


function _ff0_980(s::Number)
    m = 0.965;
    gPi = 0.165;
    rK = 4.21;
    qsq = λ(s,mπ2,mπ2)/(4*s);
    qsqK = λ(s,mK2,mK2)/(4*s);
    return 1.0/(m*m-s-2.0im*(sqrt(qsq)*gPi+sqrt(convert(Complex{Float64},qsqK))*gPi*rK)/sqrt(s));
end

function _ff0_1500(s::Number)
    m = 1.507;
    G = 0.109;
    return m*G/(m*m-s-1im*m*G);
end

function _fρ(s::Number)
    const mρ = 0.7685; const Γρ = 0.1507;
    # break up momentum
    p  = sqrt(λ(mπ2,mπ2,s)/(4*s));
    p0 = sqrt(λ(mπ2,mπ2,mρ^2))/(2*mρ);
    # ph.sp.
    ρ  = 1/(8*π)*2*p/sqrt(s)
    ρ0 = 1/(8*π)*2*p0/mρ
    # extra factor due to the spin of ρ
    R = 4.94 # it was 5
    ff = p^2/p0^2*(1.0/R^2+p0^2)/(1.0/R^2+p^2)
    mΓ = mρ*Γρ*p/p0*ff  #changed from ρ/ρ0 to p/p0
    return sqrt(ff)/(mρ^2-s -1.0im*mΓ) #
end

function _ff2(s::Number)
    mπ = 0.13956755;
    mπ2 = mπ^2;
    m = 1.274;  # it was 1.2754;
    G = 0.185;  # it was 0.1852;
    qsq_R = 1.0/4.94^2; # it was 5
    qsq = λ(s,mπ2, mπ2)/(4*s);
    qsq0 = λ(m^2,mπ2, mπ2)/(4*m^2);
    ff = sqrt(qsq/qsq0)*m/sqrt(s)*BlttWskpf[2](qsq/qsq_R)/BlttWskpf[2](qsq0/qsq_R);
    return m*G/(m*m-s-1im*m*G*ff) * sqrt(BlttWskpf[2](qsq/qsq_R));
end

#
function _fρ3(s::Number)
  mπ = 0.13956755; mπ2 = mπ^2;
  m = 1.690;
  G = 0.190;
  qsq_R = 1/4.94^2;
  qsq = λ(s,mπ2, mπ2)/(4*s);
  return sqrt(m*G*sqrt(s))/(m*m-s-1im*m*G) * sqrt(BlttWskpf[3](qsq/qsq_R));
end

for name in ["fρ3", "ff2", "fρ", "ff0_1500", "fσ", "ff0_980"]
    @eval const $(Symbol("norm_"*name)) = sqrt(quadgk(
        x->abs2($(Symbol("_"*name))(x))*sqrt((1-4mπ2/x)),4mπ2, Inf)[1]/(16*π^2))
    @eval $(Symbol(name))(σ) = $(Symbol("_"*name))(σ)/$(Symbol("norm_"*name))
end

function ZϵDD(J::Int64,M::Int64,ϵ::Bool,L::Int64,l::Int64, Dϵ, D)
    rng_lm = min(l,J);
    # println("$J,$M,$ϵ")
    sum(ClebschGordon(L,0,l,lm,J,lm)*sqrt((2*L+1)*(2*l+1))*
                Dϵ[J+1,M+1,1+ϵ,4+lm]*D[1+l,4+lm]
        for lm=-min(l,J):min(l,J))
end

isobarsV = [fσ,fρ,ff2,fρ3]
isobarsS = [fσ,ff0_980,ff0_1500]

# pwd()

# basis = []
# wavenames = []
# Nwaves = 0
# wavesfile = []

# constract COMPASS basis

function get_wavelist(path_to_wavelist; path_to_thresholds="", M3pi=3.0)
    (! isfile(path_to_wavelist)) && error("Cannot find $(path_to_wavelist)!")

    wavesInFile = readdlm(path_to_wavelist)

    thresholds = fill(0.0,size(wavesInFile,1))
    if isfile(path_to_thresholds)
        thf = readdlm(path_to_thresholds)
        for v in zip(thf[:,1],thf[:,2])
            thresholds[Int64(v[1])] = v[2]
        end
    else
        warn("Do not consider thresholds!")
    end
    thresholds_filter = thresholds .< M3pi
    wavesfile = wavesInFile[thresholds_filter,:]

    return wavesfile
end

function get_wavenames(wavelist::Array)
    return vcat(wavelist[:,2]...)::Vector{SubString{String}}
end

function get_wavenames(path_to_wavelist; path_to_thresholds="", M3pi=3.0)
    wavelist = get_wavelist(path_to_wavelist; path_to_thresholds=path_to_thresholds, M3pi=M3pi);
    get_wavenames(wavelist)
end

function get_wavebasis(wavelist::Array)
    basis = []
    for i in 1:size(wavelist,1)
        wl = wavelist[i,:]
        # special case
        if wl[2] == "FLAT"
            flat(σ1,cosθ1,ϕ1,cosθ23,ϕ23,s) = 1.0+0.0im
            push!(basis,flat)
            continue
        end
        # special case
        wn, name, J, P, M, ϵ, S, L = wl
        fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
        S = (S≥0 ? S : 0)
        push!(basis, get_compass_jmels_basis_psi(J,M,1*(P==ϵ),L,S))
    end
    return basis
end

function get_compass_jmels_basis_psi(J,M,Pϵ,L,S)
    function wave(σ1,cosθ1,ϕ1,cosθ23,ϕ23,s)
        compass_jmels_basis_psi(QNs=(J,M,Pϵ,L,S),
                                τ1=(σ1,cosθ1,ϕ1,cosθ23,ϕ23),
                                s=s)
    end
    return wave
end

function compass_jmels_basis_psi(;QNs::NTuple{5,Int}=error("give quantum numbers for the wave (J,M,Pϵ,L,S)"),
                                   τ1::NTuple{5,Float64}=error("give kinematic variables (σ1,cosθ1,ϕ1,cosθ23,ϕ23)"),
                                   s::Float64=error("give s=M3π^2"))
    J,M,Pϵ,L,S = QNs
    m1sq, m2sq, m3sq = fill(mπ2,3)
    ###############
    fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
    S = (S ≥ 0 ? S : 0)
    ###############
    τ3 = change_basis(τ1...,mπ2,mπ2,mπ2,s)
    (abs(τ3[2]) ≈ 1.0) && (τ3[2] = sign(τ3[1])*1.0);
    (abs(τ3[4]) ≈ 1.0) && (τ3[4] = sign(τ3[3])*1.0);
    (abs(τ3[2]) > 1.0 || abs(τ3[4]) > 1.0) && error("Something is wrong! (abs($(τ3[1])) > 1.0 || abs($(τ3[3]) > 1.0)")
    #############
    σ1 = τ1[1]
    σ3 = τ3[1]
    R = 1/0.2024;
    bw1 = (L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ1,m1sq)/(4s)*R^2)
    bw3 = (L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ3,m3sq)/(4s)*R^2)
    #############
    # τ1_rev = collect(τ1); τ1_rev[4] *= -1; τ1_rev[5] += π  # compass `convension` (inconsistent, btw)
    return Z(QNs...,τ1[2],τ1[3],-τ1[4],τ1[5]+π)*fi(σ1)*sqrt(bw1) +
           Z(QNs...,τ3[2],τ3[3], τ3[4],τ3[5])*fi(σ3)*sqrt(bw3)
end

#############################################################################
#############################################################################
#############################################################################

# short_basis = []
# let flat(lm,σ1,cosθ23,m1sq,m2sq,m3sq,s) = 1.0+0.0im
#     push!(short_basis,flat)
# end
# for i in 2:Nwaves
#     wn, name, J, P, M, ϵ, S, L = wavesfile[i,:]
#     fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
#     S = (S≥0 ? S : 0)
#     @eval function $(Symbol("sort_wave_$(wn)"))(lm,σ1,cosθ23,m1sq,m2sq,m3sq,s)
# #         println("\nwave ",$J," ",$P," ",$M," ",$ϵ," ",$L," ",$S)
#         # if (sqrt(σ1) < (sqrt(m2sq)+sqrt(m3sq))) || ((sqrt(σ1) > (sqrt(s)-sqrt(m1sq))))
#         #     error("Check your masses. There is inconsitency, probably.")
#         # end
#         pars = change_basis(σ1,cosθ23,m1sq,m2sq,m3sq,s)
#         σ3 = pars[1]; cosθ3 = pars[2]; cosθ12 = pars[3];
#         θ23 = acos(cosθ23); θ3 = acos(cosθ3); θ21 = acos(-cosθ12);
#
#         R = 1/0.2024; #it was 5
#         bw1 = ($L == 0) ? 1.0 : BlttWskpf[$L]($λ(s,σ1,m1sq)/(4s)*R^2)
#         bw3 = ($L == 0) ? 1.0 : BlttWskpf[$L]($λ(s,σ3,m3sq)/(4s)*R^2)
#
#         # lm is fixed!
#         rng_lm = min($S,$J);
#         res_λ1 = Wignerd($(S),lm,0,θ23)*ClebschGordon(2*$(L),0,2*$(S),2*lm,2*$(J),2*lm)
#         res_λ3 = 0.0im;
#         for ν=-rng_lm:1:rng_lm
#             res_λ3 += (ν%2==0 ? 1.0 : -1.0)*
#                 Wignerd($(J), lm, ν, θ3 )*
#                 Wignerd($(S),  ν, 0, θ21)*ClebschGordon(2*$(L),0,2*$(S),2*ν,2*$(J),2*ν)
#         end
#         res = sqrt((2*$(L)+1)*(2*$(S)+1)/(2*$(J)+1))*( #*WignerDϵ(($P==$ϵ),$J,$M,lm,0.0,0.0,0.0)
#             res_λ1*$(fi)(σ1)*sqrt(bw1)+
#             res_λ3*$(fi)(σ3)*sqrt(bw3)
#         )
#         return res
#     end
#     @eval push!(short_basis, $(Symbol("wave_$(wn)")))
# end

# function COMPASS_wave(i,s::Float64,σ1::Float64,
#     cosθ1::Float64,ϕ1::Float64,cosθ23::Float64,ϕ23::Float64)
#     (i <  1) && return 0;
#     (i > Nwaves) && return 0;
#     basis[i](σ1,cosθ1,ϕ1,cosθ23,ϕ23,mπ2,mπ2,mπ2,s)
# end

# function COMPASS_wave_short(i,lm,s::Float64,σ1::Float64,cosθ23::Float64)
#     (i <  1) && return 0;
#     (i > Nwaves) && return 0;
#     short_basis[i](lm,σ1,cosθ23,mπ2,mπ2,mπ2,s)
# end


# function COMPASS_waves(s,σ1,cosθ1,ϕ1,cosθ23,ϕ23)
#     m1sq = mπ2; m2sq = mπ2; m3sq = mπ2;
#     # system 1 <-> 23
#     const Dϵ1 = [(M > J || λ > J) ? 0.0im : WignerDϵ(ϵ==1,J,M,λ,ϕ1,acos(cosθ1),0) for J=0:6, M=0:2, ϵ=0:1, λ=-3:3]
#     const D1  = [(λ > S) ? 0.0im : WignerD(S,λ,0,ϕ23,acos(cosθ23),0) for S=0:6, λ=-3:3]
#     # Z functions
#     const Zϵf1 = Vector{Complex{Float64}}(Nwaves)
#     Zϵf1[1] = 0.5+0.0im;
#     const iV = [f(σ1) for f in isobarsV]
#     const iS = [f(σ1) for f in isobarsS]
#     const R = 1/0.2024; #it was 5
#     const bw1 = [(L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ1,m1sq)/(4s)*R^2) for L=0:6];
#     for i in 2:Nwaves
#         wn, name, J, P, M, ϵ, S, L = wavesfile[i,:]
#         fi = (S ≥ 0) ? iV[S+1] : iS[1-S]
#         S = (S≥0 ? S : 0)
#         Zϵf1[i] = sum(ClebschGordon(L,0,S,lm,J,lm)*sqrt((2*L+1)*(2*S+1))*
#                         Dϵ1[J+1,M+1,1+(ϵ==P),4+lm]*D1[1+S,4+lm]
#                          for lm=-min(S,J):min(S,J))*
#                 fi*sqrt(bw1[L+1]);
#     end
#     # system 3 <-> 12
#     τ3 = Vector{Float64}(4)
#     σ3,τ3[1],τ3[2],τ3[3],τ3[4] = change_basis(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s)
#     τ3[3] *= -1; τ3[4] += π
#     # functions
#     const Dϵ3 = [(M > J || λ > J) ? 0.0im : WignerDϵ(ϵ==1,J,M,λ,τ3[2],acos(τ3[1]),0) for J=0:6, M=0:2, ϵ=0:1, λ=-3:3]
#     const D3  = [(λ > S) ? 0.0im : WignerD(S,λ,0,τ3[4],acos(τ3[3]),0) for S=0:6, λ=-3:3]
#     # Z functions
#     const Zϵf3 = Vector{Complex{Float64}}(Nwaves)
#     Zϵf3[1] = 0.5+0.0im;
#     iV .= [f(σ3) for f in isobarsV]
#     iS .= [f(σ3) for f in isobarsS]
#     const bw3 = [(L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ3,m3sq)/(4s)*R^2) for L=0:6];
#     for i in 2:Nwaves
#         wn, name, J, P, M, ϵ, S, L = wavesfile[i,:]
#         fi = (S ≥ 0) ? iV[S+1] : iS[1-S]
#         S = (S≥0 ? S : 0)
#         Zϵf3[i] = sum(ClebschGordon(L,0,S,lm,J,lm)*sqrt((2*L+1)*(2*S+1))*
#                             Dϵ3[J+1,M+1,1+(ϵ==P),4+lm]*D3[1+S,4+lm] for lm=-min(S,J):min(S,J))*
#                         fi*sqrt(bw3[L+1]);
#     end
#     (Zϵf1 + Zϵf3)
#     # Zϵf1
# end

function swap_kin_parameters(s,τ1...)
    τ3 = collect(change_basis(τ1...,mπ2,mπ2,mπ2,s))
    τ3[4] *= -1 # angle of the oposite particle
    τ3[5] += π # angle of the oposite particle
    (τ3[5] > π) && (τ3[5] -= 2π)
    vcat([s],τ3)
end

# swap_kin_parameters(1.3,0.7,0.1,0.1,0.1,0.1)

# COMPASS_waves(1.4,0.6,0.1,0.1,0.1,0.1) - [COMPASS_wave(i,1.4,0.6,0.1,0.1,0.1,0.1) for i in 1:88]
# # WignerDϵ(true,1,1,1,0.1,0.1,0.0)
# @time for k = 1:1000
#     COMPASS_waves(1.4,0.6,0.1,0.1,0.1,0.1)
# end
# # # # # #
# @time for k = 1:1000
#     [COMPASS_wave(i,1.4,0.6,0.1,0.1,0.1,0.1) for i in 1:88]
# end
# #
# using Traceur
# @trace COMPASS_waves(1.4,0.6,0.1,0.1,0.1,0.1)
#
# using BenchmarkTools
# @btime [COMPASS_wave(i,1.4,0.6,0.1,0.1,0.1,0.1) for i in 1:88]
# @btime COMPASS_waves(1.4,0.6,0.1,0.1,0.1,0.1)

end
