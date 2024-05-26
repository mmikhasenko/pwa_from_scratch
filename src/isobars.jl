module Masses
const global mπ = 0.139
const global mπ2 = mπ^2
const global mρ = 0.7755
const global mρ2 = mρ^2
const global mτ = 1.776
const global mτ2 = mτ^2
export mπ, mπ2, mρ, mρ2, mτ, mτ2
end

module isobars
using ..Masses: mπ, mπ2
export λ, U0, U1, f1_I, f1_II, r1

"""
    λ(x,y,z)

Kaellen triangle function defined by
```
λ(x,y,z) = x^2+y^2+z^2-2*x*y-2*y*z-2*z*x
```
"""
λ(x::Number, y::Number, z::Number) = x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x

const mρ = 0.77526 # from PDG2016
const Γρ = 0.1491 #  from PDG2016

for f in [["I", "-"], ["II", "+"]]
    @eval begin
        function $(Symbol("f1r_" * f[1]))(s::Number)
            # break up momentum
            p = sqrt(λ(mπ2, mπ2, s) / (4 * s))
            p0 = sqrt(λ(mπ2, mπ2, mρ^2)) / (2 * mρ)
            # ph.sp.
            ρ = 1 / (8 * π) * 2 * p / sqrt(s)
            ρ0 = 1 / (8 * π) * 2 * p0 / mρ
            # extra factor due to the spin of ρ
            R = 5
            ff = p^2 / p0^2 * (1 ./ R^2 + p0^2) / (1 ./ R^2 + p^2)
            mΓ = mρ * Γρ * ρ / ρ0 * ff
            return sqrt(ff) / (mρ^2 - s + $(Symbol(f[2]))(1.0im * mΓ)) #
        end
    end
end

function r1r(s::Number)
    # break up momentum
    p = sqrt(λ(mπ2, mπ2, s) / (4 * s))
    p0 = sqrt(λ(mπ2, mπ2, mρ^2)) / (2 * mρ)
    # ph.sp.
    ρ0 = 1 / (8 * π) * 2 * p0 / mρ
    # extra factor due to the spin of ρ
    R = 5
    ff = p^2 / p0^2 * (1 ./ R^2 + p0^2) / (1 ./ R^2 + p^2)
    return 2mρ * Γρ / ρ0 * sqrt(ff)
end

"""
    U1r(s)
Parameterization of 'rho(770)'-meson via the Breit-Wigner function.
The function has arbitrary normalization. Use `f1(s)` for the
`quasi-two-body` normalization
"""
U1r(s::Number) = f1r_I(s) * f1r_II(s) * sqrt(λ(mπ2, mπ2, s)) / (8π * s)


### sigma
for f in [["I", "-"], ["II", "+"]]
    @eval begin
        function $(Symbol("f0r_" * f[1]))(s::Number)
            mσ = 0.86 # consistent with CLEO[2000]
            Γσ = 0.88 # from Tornquist[???]
            ρ = 1 / (8 * π) * sqrt(λ(mπ^2, mπ^2, s)) / s
            ρ0 = 1 / (8 * π) * sqrt(λ(mπ^2, mπ^2, mσ^2)) / mσ^2
            gsq = 2Γσ * mσ / ρ0
            mΓ = gsq * ρ / 2
            return 1.0 / (mσ^2 - s + $(Symbol(f[2]))(1.0im * mΓ))
        end
    end
end

"""
    U0r(s)
Parameterization of 'sigma(770)'-meson via the Breit-Wigner function.
The function has arbitrary normalization. Use `f1(s)` for the
`quasi-two-body` normalization
"""
U0r(s::Number) = f0r_I(s) * f0r_II(s) * sqrt(λ(mπ2, mπ2, s)) / (8π * s)

using QuadGK: quadgk

# normalization
for i in ["0", "1"]
    @eval begin
        const $(Symbol("ns" * i)) = quadgk(t -> begin
                s = 4mπ2 + tan(t)
                real($(Symbol("U" * i * "r"))(s) / cos(t)^2 / (2π))
            end, 0, π / 2)[1]
    end
end

# normalized functions
"""
    U1(s)
Parameterization of 'sigma(800)'-meson via the Breit-Wigner function.
It can be referenced to Tornquist.
The function has the `quasi-two-body` normalization.
"""
U1(s::Number) = U1r(s) / ns1
f1_I(s::Number) = f1r_I(s) / sqrt(ns1)
f1_II(s::Number) = f1r_II(s) / sqrt(ns1)
r1(s::Number) = r1r(s) * sqrt(ns1)

"""
    U0(s)
Parameterization of 'sigma(800)'-meson via the Breit-Wigner function.
It can be referenced to Tornquist.
The function has the `quasi-two-body` normalization.
"""
U0(s::Number) = U0r(s) / ns0
f0_I(s::Number) = f0r_I(s) / sqrt(ns0)
f0_II(s::Number) = f0r_II(s) / sqrt(ns0)

end
