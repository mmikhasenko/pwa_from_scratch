# module masses
# global const mπ=0.139; global const mπ2=mπ^2;
# global const mρ=0.7755; global const mρ2=mρ^2;
# global const mτ=1.776; global const mτ2=mτ^2;
# export mπ, mπ2, mρ, mρ2, mτ, mτ2
# end
#
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
const _a = [0.1131, 0.1968^2]

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
let flat(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s) = 1
    push!(basis,flat)
    push!(wavenames,"flat")
end
for i in 2:size(wavesload,1)
    wn, name, J, P, M, ϵ, S, L = wavesload[i,:]
    fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
    S = (S≥0 ? S : 0)
    @eval function $(Symbol("wave_$(wn)"))(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s)
#         println("\nwave ",$J," ",$P," ",$M," ",$ϵ," ",$L," ",$S)
        τ3 = [0.0,0.0,0.0,0.0]
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

function COMPASS_wave(i,s,σ1,cosθ1,ϕ1,cosθ23,ϕ23)
    (i <  1) && return 0;
    (i > 88) && return 0;
    basis[i](σ1,cosθ1,ϕ1,cosθ23,ϕ23,mπ2,mπ2,mπ2,s)
end

function COMPASS_waves(s,σ1,cosθ1,ϕ1,cosθ23,ϕ23)
    m1sq = mπ2; m2sq = mπ2; m3sq = mπ2;
    # system 1 <-> 23
    Dϵ1 = [(M > J || λ > J) ? 0.0 : WignerDϵ(ϵ==1,J,M,λ,ϕ1,acos(cosθ1),0) for J=0:6, M=0:2, ϵ=0:1, λ=-3:3]
    D1  = [(λ > S) ? 0.0 : WignerD(S,λ,0,ϕ23,acos(cosθ23),0) for S=0:6, λ=-3:3]
    # Z functions
    Zϵf1 = [0.5+0.0im];
    for i in 2:size(wavesload,1)
        wn, name, J, P, M, ϵ, S, L = wavesload[i,:]
        fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
        S = (S≥0 ? S : 0)
        R = 5;
        bw1 = (L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ1,m1sq)/(4s)*R^2);
        push!(Zϵf1, sum(ClebschGordon(L,0,S,lm,J,lm)*sqrt((2*L+1)*(2*S+1))*
                        Dϵ1[J+1,M+1,1+(ϵ==P),4+lm]*D1[1+S,4+lm]
                         for lm=-min(S,J):min(S,J))*
                fi(σ1)*bw1);
    end
    # system 3 <-> 12
    τ3 = [0.0,0.0,0.0,0.0]
    σ3,τ3[1],τ3[2],τ3[3],τ3[4] = change_basis(σ1,cosθ1,ϕ1,cosθ23,ϕ23,m1sq,m2sq,m3sq,s)
    τ3[3] *= -1; τ3[4] += π
    # functions
    Dϵ3 = [(M > J || λ > J) ? 0.0 : WignerDϵ(ϵ==1,J,M,λ,τ3[2],acos(τ3[1]),0) for J=0:6, M=0:2, ϵ=0:1, λ=-3:3]
    D3  = [(λ > S) ? 0.0 : WignerD(S,λ,0,τ3[4],acos(τ3[3]),0) for S=0:6, λ=-3:3]
    # Z functions
    Zϵf3 = [0.5+0.0im];
    for i in 2:size(wavesload,1)
        wn, name, J, P, M, ϵ, S, L = wavesload[i,:]
        fi = (S ≥ 0) ? isobarsV[S+1] : isobarsS[1-S]
        S = (S≥0 ? S : 0)
        R = 5;
        bw3 = (L == 0) ? 1.0 : BlttWskpf[L](λ(s,σ3,m3sq)/(4s)*R^2);
        push!(Zϵf3, sum(ClebschGordon(L,0,S,lm,J,lm)*sqrt((2*L+1)*(2*S+1))*
                            Dϵ3[J+1,M+1,1+(ϵ==P),4+lm]*D3[1+S,4+lm] for lm=-min(S,J):min(S,J))*
                        fi(σ3)*bw3);
    end
    (Zϵf1 + Zϵf3)
    # Zϵf1
end

# WignerDϵ(true,1,1,1,0.1,0.1,0.0)o
# COMPASS_waves(1.4,0.6,0.1,0.1,0.1,0.1)

end


# For three years I was deeply involved in the data reconstruction and analysis at VES experiment. My work was focused on the selection and analysis of the pi-pi0 system produced via peripheral scattering of the pion beam off a nucleon in the Beryllium target.
# It was required select good event candidates from 10^7 data sample collected in 2011-2012. Neutral pion was reconstructed by the fit of energy clusters in the electromagnetic calorimeter. The shower profile was measured in with electron beam, calibration
# The 1C-kinematical constraint for the pi0 mass was applied.



# 1) Physics simulations
#   I started my research with the physics simulation project which was suggested as a bachelor topic. The studies were dedicated to the idea of the active target detector for VES experiment. A bundle of the scintillating fibers was supposed to replace the target for the scattering experiment aiming to detect the recoil system.
#   Later, I worked on Geant4-simulations of the lead converters to improve acceptance for the neutral particles at VES experiment.
#   During my PhD, I had a chance to run Geant4 based COMPASS simulations. I implemented new event generator and advised a bachelor project for acceptance studies for the reaction pi- p -> pi- eta p at COMPASS.
#
# 2) Data reduction and numerical analysis
#   During my master studies, I was deeply involved in the data reconstruction and analysis at VES experiment. My main project was the selection and analysis of the pi-pi0 system produced via peripheral scattering of the pion beam off a nucleon in the Beryllium target. It was required to select good event candidates from 10^7 data sample collected in 2011-2012.
#
#
# 3)
#
# Neutral pion was reconstructed from showers information in electromagnetic calorimeter
# The 1C-kinematical constraint for the pi0 mass was applied.
#
# The largest computation
#
# List of reactions I worked on, referenced below.
# (1) pi-N->pi-pi0 N',
# (2) pi-p->pi-eta p', eta->gg, eta->3pi,
# (3a/b) pi-p->3pi/pi-2pi0 p',
# (4a/b) tau -> 3p/pi-2pi0i nu,
# (5) B->J/Psi pi K
# (6) e+e-->J/Psi pi+pi-
# 
# 1) Physics simulations
#   = my bachelor project was dedicated to the studies of the active target detector. I worked on the Geant4 model and reconstruction algorithm for the bundle of the scintillating fibers used as the recoil detection system.
#   = I worked on Geant4-simulations of the lead converters to improve acceptance for the neutral particles at VES experiment.
#   = Acceptance studies for (1) was a part of my master project.
#   = Geant4 based COMPASS simulations. I implemented new event generator and advised a bachelor project for acceptance studies for the reaction (2) at COMPASS.
#
# 2) Data reduction and numerical analysis
#   = I was involved in the reconstruction and analyses at VES experiment. I performed the event selection and the partial wave analysis of the pi-pi0 in reaction (1). The significant background from the reaction (3) was simulated and subtracted.
#   = My work at COMPASS was focused on the numerical analyses related to the reaction (3a). I tested various theoretical models to the dynamics in partial waves on the COMPASS data. Pole positions for a1 and pi2 mesons were obtained. The exotic meson candidate was described as the rescattering phenomenon.
#
# 3) Theoretical calculations:
#   = Rescattering corrections in the system of three pions for the reaction (3a).
# Effective Lagrangian and dispersive theory.
#   = Tensor and helicity formalisms. The amplitude model for the reaction (5).
#
# 4) Teaching and advising
#   = Lecturer at the Reaction Theory School
#   = Advising 3 master- and 1 bachelor- students
#
#
#   "Nature of the a1(1420)", Phys. Rev. D91.9, p. 094015, 3 authors, inspire id 1341619.
#   personal contribution: all calculations, text
#
#   "What is the right formalism to search for resonances?" arXiv: 1712.02815 [hep-ph], accepted by EPJC, 10 authors, inspire id 1642229
#   personal contribution: all calculations, text
#
#   "New analysis of etapi tensor resonances measured at the COMPASS
#   experiment", arXiv: 1707.02848 [hep-ph], 235 authors, inspire id 1609260
#   personal contribution: cross-check of calculation and fit results, corrections to the text
#
#   "Amplitude analysis and the nature of the Zc(3900)", Phys. Lett.
#   B772, pp. 200–209, 7 authors, inspire id 1505197
#   personal contribution: cross-check of calculation, corrections to the text
#
#   "Light isovector resonances in pi-pp>3pi p at 190 GeV/c", arXiv:1802.05913 [hep-ex], sent to EPJC, 214 authors, inspire id 1655631
#   personal contribution: appendix A, corrections to text
#
#   "Unitarity approach to the mass-dependent fit of 3pi resonance production data from the COMPASS experiment", EPJ Web Conf. 137, p. 05017, 4 authors, inspire id 1519590
#   personal contributions: formalism, calculations, fits, text.
