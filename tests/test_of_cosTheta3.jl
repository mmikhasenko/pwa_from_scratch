using Base.Test

push!(LOAD_PATH, "src")

using DalitzPlotAnalysis

fewtests = [
    let s = 1.4, mπ = 0.14, mπ2 = mπ^2, z = 2rand() - 1
        σ1 = 4mπ2 + ((√s - mπ)^2 - 4mπ2) * rand()
        τ = change_basis(σ1, 1.0, 0.0, z, 0.0, mπ2, mπ2, mπ2, s)
        cosθ3 = cross_basis_cosθ3(σ1, z, mπ2, mπ2, mπ2, s)
        # that is the test
        τ[2] ≈ cosθ3 #
    end
    for i in 1:1000]

@test sum(fewtests) == size(fewtests, 1)
