# module masses
# global const mπ=0.139; global const mπ2=mπ^2;
# global const mρ=0.7755; global const mρ2=mρ^2;
# global const mτ=1.776; global const mτ2=mτ^2;
# export mπ, mπ2, mρ, mρ2, mτ, mτ2
# end
#
module isobars
# using masses: mπ, mπ2
export fρ, ff2, fρ3, fσ, ff0_980, ff0_1500, BlttWskpf

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

end
