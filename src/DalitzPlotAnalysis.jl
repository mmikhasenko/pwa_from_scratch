module DalitzPlotAnalysis
export ClebschGordon, Wignerd, WignerDϵ, Z, λ
export change_basis, cross_basis_cosθ3, cross_basis_cosθ12, cross_basis_s3

import Base: Math.atan2

import GSL: sf_coupling_3j

"""
    atan2(y::Complex{Float64}, x::Complex{Float64})
The function calculates `atan` for a complex argument. The standard atan(y/x)
    is either corrected by +pi or returned.

    θ=2+0.1im
    println("sin: ", cos(θ))
    println("rec: ", cos(atan(tan(θ))))
    println("+pi: ", cos(π+atan(tan(θ))))
    println("at2: ", cos(atan2(sin(θ),cos(θ))))
    println("sin: ", sin(θ))
    println("rec: ", sin(atan(tan(θ))))
    println("+pi: ", sin(π+atan(tan(θ))))
    println("at2: ", sin(atan2(sin(θ),cos(θ))))
"""
function atan2(y::Complex{Float64}, x::Complex{Float64})
    val = atan(y/x);
    cosval = cos(val)
    return (real(cosval)*real(x) > 0) ? val : val+π;
end

λ(x,y,z) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x


for tp in [:Float64, :(Complex{Float64})]
    @eval function change_basis(s1::$(tp), cosθ23::Float64, m1sq::Float64, m2sq::Float64, m3sq::Float64, s::$(tp))
            # calculate s3 in (23) frame
            s3 = m1sq + m2sq +
            # 2*(  E2 * E1 -
                2*( (s1 + m2sq - m3sq)/(2*sqrt(s1)) * (s-m1sq-s1)/(2*sqrt(s1)) -
                # |p2|*cos(theta23) * (-|p1|)
                sqrt(λ(s1, m2sq, m3sq)/(4*s1))*cosθ23 * (-sqrt(λ(s, m1sq, s1)/(4*s1))) );
            # caluclate p3 in (23) frame
            p23bu = sqrt(λ(s1, m2sq, m3sq)/(4*s1));
            p3_in23_E = (s1 + m3sq - m2sq)/sqrt(4*s1);
            p3_in23_z = -p23bu*cosθ23;
            # boost to lab frame from (23) frame
            γ1 = (s + s1 - m1sq)/sqrt(4*s*s1);
            β1 = sqrt(1.-1./γ1^2);
            p3_b_z = γ1*(β1*p3_in23_E+p3_in23_z);

            p3_cms = sqrt(λ(s,m3sq,s3))/sqrt(4s)
            cosθ3 = -p3_b_z/p3_cms;

            cosθ12_n = m2sq + m3sq + 2* (s3+m2sq-m1sq)/(2*sqrt(s3)) * (s-m3sq-s3)/(2*sqrt(s3)) - s1;
            cosθ12_d = (2 * sqrt(λ(s3, m1sq, m2sq)/(4*(s3))) * sqrt(λ(s, m3sq, s3)/(4*(s3))) );
            cosθ12 = cosθ12_d ≈ 0.0+0.0im ? 2.0*rand()-1.0*one(s1) : cosθ12_n/cosθ12_d;

            return s3, cosθ3, cosθ12
    end
end

function cross_basis_cosθ3(s1, cosθ23, m1sq, m2sq, m3sq, s)

     s3 = cross_basis_s3(s1, cosθ23, m1sq, m2sq, m3sq, s);
     s2 = s + m1sq + m2sq + m3sq - s1 - s3;
     # check that there will not be the Domain Error.
     λλ = λ(s, m1sq, s1)*λ(s, m3sq, s3)
     (λλ < 0) && error("λλ < 0, check masses")
     # directly get the angle
     return (2s*(m1sq+m3sq-s2)+(s+m1sq-s1)*(s+m3sq-s3))/sqrt(λλ)
end

function cross_basis_s3(s1, cosθ23, m1sq, m2sq, m3sq, s)
     # calculate s3 in (23) frame
     return m1sq + m2sq +
     # 2*(  E2 * E1 -
     2*( (s1 + m2sq - m3sq)/(2*sqrt(s1)) * (s-m1sq-s1)/(2*sqrt(s1)) -
     # |p2|*cos(theta23) * (-|p1|)
     sqrt(λ(s1, m2sq, m3sq)/(4*s1))*cosθ23 * (-sqrt(λ(s, m1sq, s1)/(4*s1))) )
end

function cross_basis_cosθ12(s1, cosθ23, m1sq, m2sq, m3sq, s)
     s3 = m1sq + m2sq +
     # 2*(  E2 * E1 -
         2*( (s1 + m2sq - m3sq)/(2*sqrt(s1)) * (s-m1sq-s1)/(2*sqrt(s1)) -
         # |p2|*cos(theta23) * (-|p1|)
         sqrt(λ(s1, m2sq, m3sq)/(4*s1))*cosθ23 * (-sqrt(λ(s, m1sq, s1)/(4*s1))) );

     cosθ12_n = m2sq + m3sq + 2* (s3+m2sq-m1sq)/(2*sqrt(s3)) * (s-m3sq-s3)/(2*sqrt(s3)) - s1;
     cosθ12_d = (2 * sqrt(λ(s3, m1sq, m2sq)/(4*(s3))) * sqrt(λ(s, m3sq, s3)/(4*(s3))) );
     # cosθ12 = cosθ12_d ≈ 0.0+0.0im ? 2.0*rand()-1.0*one(s1) : cosθ12_n/cosθ12_d;
     # cosθ12 = 2.0*rand()-1.0*one(s1)
     (cosθ12_d ≈ zero(s)) && error("λλ == 0: check your masses")
     cosθ12 = cosθ12_n/cosθ12_d

     return cosθ12
end


# recoupling function
# two options for the type:
#   - either all arguments are Float64 or energy variables are complex
for tp in [:Float64, :(Complex{Float64})]
    @eval function change_basis(s1::$(tp), cosθ1::Float64, ϕ1::Float64, cosθ23::Float64, ϕ23::Float64,
        m1sq::Float64, m2sq::Float64, m3sq::Float64, s::$(tp))
        # calculate s3 in (23) frame
        s3 = m1sq + m2sq +
        # 2*(  E2 * E1 -
            2*( (s1 + m2sq - m3sq)/(2*sqrt(s1)) * (s-m1sq-s1)/(2*sqrt(s1)) -
            # |p2|*cos(theta23) * (-|p1|)
            sqrt(λ(s1, m2sq, m3sq)/(4*s1))*cosθ23 * (-sqrt(λ(s, m1sq, s1)/(4*s1))) );
        # caluclate p3 in (23) frame
        p23bu = sqrt(λ(s1, m2sq, m3sq)/(4*s1));
        p3_in23 = [(s1 + m3sq - m2sq)/sqrt(4*s1),
                   -p23bu*sqrt(1-cosθ23^2)*cos(ϕ23),
                   -p23bu*sqrt(1-cosθ23^2)*sin(ϕ23),
                   -p23bu*cosθ23];
        # boost to lab frame from (23) frame
        γ1 = (s + s1 - m1sq)/sqrt(4*s*s1);
        β1 = sqrt(1.-1./γ1^2);
        p3_b = [γ1*(p3_in23[1]+β1*p3_in23[4]),
                p3_in23[2],
                p3_in23[3],
                γ1*(β1*p3_in23[1]+p3_in23[4])];

        ct1 = cosθ1; st1 = sqrt(1.-cosθ1^2); cp1 = cos(ϕ1); sp1 = sin(ϕ1);
        # Rz(phi1) * Ry(theta1) * p3_boost
        p3_rot = [p3_b[1],
                  cp1*ct1* p3_b[2] + (-sp1)* p3_b[3] + cp1*st1* p3_b[4],
                  sp1*ct1* p3_b[2] +   cp1 * p3_b[3] + sp1*st1* p3_b[4],
                     -st1* p3_b[2] +     0 * p3_b[3] +     ct1* p3_b[4]];

        cosθ3 = -p3_rot[4]/sqrt(p3_b[2]^2+p3_b[3]^2+p3_b[4]^2);
        ϕ3 = (p3_rot[2] != zero(p3_rot[2])) ? atan2(-p3_rot[3], -p3_rot[2]) : rand()*one(p3_rot[2]);

        cosθ12_n = m2sq + m3sq + 2* (s3+m2sq-m1sq)/(2*sqrt(s3)) * (s-m3sq-s3)/(2*sqrt(s3)) - s1;
        cosθ12_d = (2 * sqrt(λ(s3, m1sq, m2sq)/(4*(s3))) * sqrt(λ(s, m3sq, s3)/(4*(s3))) );
        cosθ12 = cosθ12_d ≈ 0.0+0.0im ? 2.0*rand()-1.0*one(s1) : cosθ12_n/cosθ12_d;

        n1 = [0.,
             -sqrt(1.-cosθ1^2)*cos(ϕ1),
             -sqrt(1.-cosθ1^2)*sin(ϕ1),
             -cosθ1];
        ct3 = cosθ3; st3 = sqrt(1.-cosθ3^2); cp3 = cos(ϕ3); sp3 = sin(ϕ3);
        # Rz(phi23) * Ry(theta23) * p3_boost
        n1_rot = [n1[1],
                  cp3*ct3* n1[2] + sp3*ct3* n1[3] + (-st3)* n1[4],
                   (-sp3)* n1[2] +    cp3 * n1[3] +   0.0 * n1[4],
                  cp3*st3* n1[2] + sp3*st3* n1[3] +   ct3 * n1[4]];
        # println(n1_rot[3], " ", n1_rot[2])
        ϕ12 = (n1_rot[2] != zero(n1_rot[2])) ? atan2(n1_rot[3], n1_rot[2]) : rand()*one(n1_rot[2]);
        return s3,cosθ3,ϕ3,cosθ12,ϕ12
    end
end

# using BenchmarkTools
# @btime change_basis(0.5834185707670233 + 0.08828205993763047im, -0.8952547092523563, 3.141592653589793, 0.1, 1.5707963267948966,0.139^2,0.139^2,0.139^2,1.4 + 0.2im)
# @btime change_basis(0.5834185707670233 + 0.08828205993763047im, -0.8952547092523563, 3.141592653589793, 0.0,
#     1.5707963267948966,0.139^2,0.139^2,0.139^2,1.4 + 0.2im)

const logfact = [sum(log(n) for n in 1:i) for i in 1:50]
lf(i) = (i>0) ? logfact[i] :  0.0;

function Wignerd(aj, am, an, β)
    (β == zero(β)) && return (am==an) ? one(β) : zero(β)
    (β == zero(β)) && return (am==-an) ? ((aj-am)%2==1 ? -one(β) : one(β))  : zero(β)

    jpm = aj+am;
    jpn = aj+an;
    jmm = aj-am;
    jmn = aj-an;
    mpn = am+an;

    # common prefactor
    pref = (lf(jpm)+lf(jmm)+lf(jpn)+lf(jmn))/2.0

    #
    sb2 = sin(β/2)
    s  = log(sb2)          # log(sin(β/2))
    c  = log(1-sb2^2)/2.0  # log(cos(β/2))


    res = zero(β)
    let k0 = max(0,mpn), sgn = ((k0+jpm)%2==1 ? -1 : 1)
        for k in k0:min(jpm,jpn)
            res += sgn*exp(pref-
                    lf(jpm-k)-lf(k)-lf(k-mpn)-lf(jpn-k)+
                    (2k-mpn)*c + (jpm+jpn-2k)*s)
            sgn *= -1;
        end
    end
    res
end

# Integer argument
# Wignerd(aj::Int64, am::Int64, an::Int64, β) = Wignerd(aj//1,am//1,an//1,β)

# WignerD
WignerD(aj, am, an, α, β, γ) = Wignerd(aj,am,an,β)*cis(-(am*α+an*γ))
# WignerDϵ
function WignerDϵ(ϵ::Bool,aj::Int64, am::Int64, an::Int64, α, β, γ)
    (am < zero(am)) && return 0.0im;
    WDMp = WignerD(aj,am,an,α,β,γ)
    WDMm = (am == 0) ? WDMp : WignerD(aj,-am,an,α, β, γ);
    factor = (2ϵ-1) * (2((aj-am)%2==0)-1);
    return (WDMp - factor * WDMm) / (am == zero(am) ? 2.0 : sqrt(2.0));
end

# # Clebsches and d-function
# function ClebschGordon(j1::Rational{Int64},m1::Rational{Int64},
#                        j2::Rational{Int64},m2::Rational{Int64},
#                         j::Rational{Int64}, m::Rational{Int64})
#     factor = convert(Int64,j1+j2-m)%2==1 ? -1 : 1
#     factor*sqrt(2*j+1)*sf_coupling_3j(convert(Int64,2*j1),convert(Int64,2*j2), convert(Int64,2*j),
#                                       convert(Int64,2*m1),convert(Int64,2*m2),-convert(Int64,2*m))
# end

function ClebschGordon(two_j1::Int64,two_m1::Int64,
                       two_j2::Int64,two_m2::Int64,
                        two_j::Int64, two_m::Int64)
    factor = ((two_j1÷2)+(two_j2÷2-(two_m÷2) % 2 == 0 ? 1.0 : -1.0;
    factor*sqrt(two_j+1)*sf_coupling_3j(two_j1,two_j2,two_j,two_m1,two_m2,two_m)
end

# function Wignerd(j,m1,m2,cosθ)
#     M1=m1; M2=m2; factor=1;
#     if (m1==0 && m2!=0)
#         factor *= (abs(m2)%2==1 ? -1 : 1)
#         M1=m2;M2=m1;
#     end
#     factor *= (M1 < 0) ? (abs(M1)%2==1 ? -1 : 1) : 1
#     (M2 == 0) && return factor*sqrt(4*π/(2*j+1))*sf_legendre_sphPlm(j, abs(M1), cosθ)
# end
# function Wignerd(j,m1,m2,cosθ,ϕ)
#     Wignerd(j,m1,m2,cosθ)*cis(-m1*ϕ)
# end

function Z(J::Int64,M::Int64,L::Int64,l::Int64,cosθ1,ϕ1,cosθ23,ϕ23)
    rng_lm = min(l,J);
    θ1 = acos(cosθ1)
    θ23 = acos(cosθ23)
    sum(ClebschGordon(2*L,0,2*l,2*lm,2*J,2*lm)*sqrt((2*L+1)*(2*l+1))*
        WignerD(J,M,lm,ϕ1,θ1,0.0)*WignerD(l,lm,0,ϕ23,θ23,0.0) for lm=-rng_lm:1:rng_lm)
end


function Z(J::Int64,M::Int64,ϵ::Bool,L::Int64,l::Int64,cosθ1,ϕ1,cosθ23,ϕ23)
    rng_lm = min(l,J);
    θ1 = acos(cosθ1)
    θ23 = acos(cosθ23)
    sum(ClebschGordon(2*L,0,2*l,2*lm,2*J,2*lm)*sqrt((2*L+1)*(2*l+1))*
        WignerDϵ(ϵ,J,M,lm,ϕ1,θ1,0.0)*WignerD(l,lm,0,ϕ23,θ23,0.0) for lm=-rng_lm:1:rng_lm)
end

# Z(1,0,false,0,1,0.1,0.0,0.1,0.0)

end
