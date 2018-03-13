module amplitudes_compass
# using masses: mπ, mπ2
using DalitzPlotAnalysis: change_basis, Z, WignerDϵ, WignerD, ClebschGordon

export λ, fρ, ff2, fρ3, fσ, ff0_980, ff0_1500, BlttWskpf, COMPASS_wave, COMPASS_waves

"""
    λ(x,y,z)

Kaellen triangle function defined by
```
λ(x,y,z) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x
```
"""
λ(x::Number,y::Number,z::Number) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x

const mπ = 0.139570; const mπ2 = mπ^2;
const mK = 0.493677; const mK2 = mK^2;
const mK0 = 0.497614; const mK02 = mK0^2;
const mρ = 0.77526; # from PDG2016
const Γρ = 0.1491; #  from PDG2016

################### SIGMA #################
# AMP Table 1, M solution: f_2^2
# AMP Table 1, M solution: f_1^1 and f_2^1
# Last item is Katchaev modification
const _a = [0.1131, 0.0]

# AMP Table 1, M solution: c_11^0
# last item - Katchaev modification
const _c = [0.0337, -0.3185, -0.0942, -0.5927, 0.0]
const _sP = [-0.0074, 0.9828]

function fσ(s::Number)
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

  return 1./(M00 - 1im*rho00);
end


function ff0_980(s::Number)
    m = 0.965;
    gPi = 0.165;
    rK = 4.21;
    qsq = λ(s,mπ2,mπ2)/(4*s);
    qsqK = λ(s,mK2,mK2)/(4*s);
    return 1./(m*m-s-2.0im*(sqrt(qsq)*gPi+sqrt(convert(Complex{Float64},qsqK))*gPi*rK)/sqrt(s));
end

function ff0_1500(s::Number)
    m = 1.507;
    G = 0.109;
    return m*G/(m*m-s-1im*m*G);
end

function fρ(s::Number)
    # break up momentum
    p  = sqrt(λ(mπ2,mπ2,s)/(4*s));
    p0 = sqrt(λ(mπ2,mπ2,mρ^2))/(2*mρ);
    # ph.sp.
    ρ  = 1/(8*π)*2*p/sqrt(s)
    ρ0 = 1/(8*π)*2*p0/mρ
    # extra factor due to the spin of ρ
    R = 5
    ff = p^2/p0^2*(1./R^2+p0^2)/(1./R^2+p^2)
    mΓ = mρ*Γρ*ρ/ρ0*ff
    return sqrt(ff)/(mρ^2-s -1.0im*mΓ) #
end

const BlttWskpf = [z->z/(1+z),
                   z->z^2/(9+3z+z^2),
                   z->z*z*z/(z*z*z+6.*z*z+45.*z+225.),
                   z->z*z*z*z/(z*z*z*z+10.*z*z*z+135.*z*z+1575.*z+11025.),
                   z->z*z*z*z*z/(z*z*z*z*z+15.*z*z*z*z+315.*z*z*z+6300.*z*z+99225.*z+893025.),
                   z->z*z*z*z*z*z/(z*z*z*z*z*z+21.*z*z*z*z*z+630.*z*z*z*z+18900.*z*z*z+496125.*z*z+9823275.*z+108056025.)];

function ff2(s::Number)
    mπ = 0.13956755;
    mπ2 = mπ^2;
    m = 1.2754;
    G = 0.1852;
    qsq_R = 1.0/5^2;
    qsq = λ(s,mπ2, mπ2)/(4*s);
    qsq0 = λ(m^2,mπ2, mπ2)/(4*m^2);
    ff = sqrt(qsq/qsq0)*m/sqrt(s)*BlttWskpf[2](qsq/qsq_R)/BlttWskpf[2](qsq0/qsq_R);
    return m*G/(m*m-s-1im*m*G*ff) * sqrt(BlttWskpf[2](qsq/qsq_R));
end

#
function fρ3(s::Number)
  mπ = 0.13956755; mπ2 = mπ^2;
  m = 1.690;
  G = 0.190;
  qsq_R = 1/4.94^2;
  qsq = λ(s,mπ2, mπ2)/(4*s);
  return sqrt(m*G*sqrt(s))/(m*m-s-1im*m*G) * sqrt(BlttWskpf[3](qsq/qsq_R));
end

function ZϵDD(J::Int64,M::Int64,ϵ::Bool,L::Int64,l::Int64, Dϵ, D)
    rng_lm = min(l,J);
    # println("$J,$M,$ϵ")
    sum(ClebschGordon(L,0,l,lm,J,lm)*sqrt((2*L+1)*(2*l+1))*
                Dϵ[J+1,M+1,1+ϵ,4+lm]*D[1+l,4+lm]
        for lm=-min(l,J):min(l,J))
end


# constract COMPASS basis
# pwd()
wavesload = readdlm(pwd()*"/src/wavelist_formated.txt")
isobarsV = [fσ,fρ,ff2,fρ3]
isobarsS = [fσ,ff0_980,ff0_1500]

basis = []
wavenames = []
let flat(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s) = 1.0+0.0im
    push!(basis,flat)
    push!(wavenames,"flat")
end
for i in 2:size(wavesload,1)
    wn, name, J, P, M, ϵ, S, L = wavesload[i,:]
    fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
    S = (S≥0 ? S : 0)
    @eval function $(Symbol("wave_$(wn)"))(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s)
#         println("\nwave ",$J," ",$P," ",$M," ",$ϵ," ",$L," ",$S)
        τ3 = Vector{Float64}(4)
        σ3,τ3[1],τ3[2],τ3[3],τ3[4] = change_basis(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s)
        τ3[3] *= -1; τ3[4] += π
        R = 5;
        bw1 = ($L == 0) ? 1.0 : BlttWskpf[$L]($λ(s,σ1,m1sq)/(4s)*R^2)
        bw3 = ($L == 0) ? 1.0 : BlttWskpf[$L]($λ(s,σ3,m3sq)/(4s)*R^2)
        return Z($J,$M,($P==$ϵ),$L,$S,cosθ1,ϕ1,cosθ23,ϕ23)*$(fi)(σ1)*bw1 +
               Z($J,$M,($P==$ϵ),$L,$S,τ3...)              *$(fi)(σ3)*bw3
    end
    push!(wavenames,name)
    @eval push!(basis, $(Symbol("wave_$(wn)")))
end

function COMPASS_wave(i,s::Float64,σ1::Float64,
    cosθ1::Float64,ϕ1::Float64,cosθ23::Float64,ϕ23::Float64)
    (i <  1) && return 0;
    (i > 88) && return 0;
    basis[i](σ1,cosθ1,ϕ1,cosθ23,ϕ23,mπ2,mπ2,mπ2,s)
end

function COMPASS_waves(s,σ1,cosθ1,ϕ1,cosθ23,ϕ23)
    m1sq = mπ2; m2sq = mπ2; m3sq = mπ2;
    # system 1 <-> 23
    const Dϵ1 = [(M > J || λ > J) ? 0.0im : WignerDϵ(ϵ==1,J,M,λ,ϕ1,acos(cosθ1),0) for J=0:6, M=0:2, ϵ=0:1, λ=-3:3]
    const D1  = [(λ > S) ? 0.0im : WignerD(S,λ,0,ϕ23,acos(cosθ23),0) for S=0:6, λ=-3:3]
    # Z functions
    const Zϵf1 = Vector{Complex{Float64}}(size(wavesload,1))
    Zϵf1[1] = 0.5+0.0im;
    const iV = [f(σ1) for f in isobarsV]
    const iS = [f(σ1) for f in isobarsS]
    const R = 5.0;
    const bw1 = [(L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ1,m1sq)/(4s)*R^2) for L=0:6];
    for i in 2:size(wavesload,1)
        wn, name, J, P, M, ϵ, S, L = wavesload[i,:]
        fi = (S ≥ 0) ? iV[S+1] : iS[1-S]
        S = (S≥0 ? S : 0)
        Zϵf1[i] = sum(ClebschGordon(L,0,S,lm,J,lm)*sqrt((2*L+1)*(2*S+1))*
                        Dϵ1[J+1,M+1,1+(ϵ==P),4+lm]*D1[1+S,4+lm]
                         for lm=-min(S,J):min(S,J))*
                fi*bw1[L+1];
    end
    # system 3 <-> 12
    τ3 = Vector{Float64}(4)
    σ3,τ3[1],τ3[2],τ3[3],τ3[4] = change_basis(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s)
    τ3[3] *= -1; τ3[4] += π
    # functions
    const Dϵ3 = [(M > J || λ > J) ? 0.0im : WignerDϵ(ϵ==1,J,M,λ,τ3[2],acos(τ3[1]),0) for J=0:6, M=0:2, ϵ=0:1, λ=-3:3]
    const D3  = [(λ > S) ? 0.0im : WignerD(S,λ,0,τ3[4],acos(τ3[3]),0) for S=0:6, λ=-3:3]
    # Z functions
    const Zϵf3 = Vector{Complex{Float64}}(size(wavesload,1))
    Zϵf3[1] = 0.5+0.0im;
    iV .= [f(σ3) for f in isobarsV]
    iS .= [f(σ3) for f in isobarsS]
    const bw3 = [(L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ3,m3sq)/(4s)*R^2) for L=0:6];
    for i in 2:size(wavesload,1)
        wn, name, J, P, M, ϵ, S, L = wavesload[i,:]
        fi = (S ≥ 0) ? iV[S+1] : iS[1-S]
        S = (S≥0 ? S : 0)
        Zϵf3[i] = sum(ClebschGordon(L,0,S,lm,J,lm)*sqrt((2*L+1)*(2*S+1))*
                            Dϵ3[J+1,M+1,1+(ϵ==P),4+lm]*D3[1+S,4+lm] for lm=-min(S,J):min(S,J))*
                        fi*bw3[L+1];
    end
    (Zϵf1 + Zϵf3)
    # Zϵf1
end

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
