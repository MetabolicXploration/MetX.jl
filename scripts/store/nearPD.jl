using LinearAlgebra: eigen, Diagonal, I, det, rank
using Plots

include("sigma.jl")

## ----------------------------------------------
let
    
    δ = 1e-10
    B0 = copy(B)
    B0 = 0.5 * (B0 + B0')
    e, U = eigen(B0)
    τ = max.(δ, real.(e))
    # B0 .= U * Diagonal(τ) * U'
    # B0 .= 0.5 * (A + A')s
    U
end

## ----------------------------------------------
@show isposdef(A)
X = round.(0.5 * (A * A'); digits = 5)
@time P, _ = prox(IndPSD(), X)
# prox!(P, IndPSD(), A) # inline version
@show isposdef(P)

λA, _ = eigen(A)
λX, _ = eigen(X)
λP, _ = eigen(P)


# p = plot()
# T(x) = log10(abs(x))
# plot!(p, T.(sort(λA)); label = "A", lw = 3)
# plot!(p, T.(sort(λX)); label = "X", lw = 3)
# plot!(p, T.(sort(λP)); label = "P", lw = 3)

nothing