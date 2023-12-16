@time begin
    using Plots
    using MetXNetHub
    import MetXBase
    using LinearAlgebra
end

# ------------------------------------------------------------------
# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function basis_rxns(X::AbstractMatrix; 
        tol::Float64 = 1e-10,
        include = nothing # TODO: implement this (try to include the given cols in the independent set)
    )

    sum(abs2,X) == 0 && return Int[]
    _,R,E = qr(X, ColumnNorm())
    diagr = diag(R)
    diagr1 = abs(first(diagr))
    r = findlast(abs.(diagr) .>= tol * diagr1)
    isnothing(r) && return Int[]
    idx = sort(view(E, 1:r))
    return idx

end

function echelonize1(X::T, v; 
        tol::Real = 1e-10
    ) where T <:DenseArray

    _, N = size(X)
    c0 = zero(eltype(X))
    c1 = one(eltype(X))

    idxrow = basis_rxns(Matrix(X'))
    Mred = length(idxrow)

    idxd = basis_rxns(X; tol)
    idxf = setdiff(1:N, idxd)
    idxmap = vcat(idxd, idxf)
    Tv = view(X, idxrow, idxmap)
    iTv = inv(Tv[1:Mred, 1:Mred])
    IG = iTv * Tv
    # trimming zeros
    @inbounds for i in eachindex(IG)
        abs(IG[i]) < tol && (IG[i] = c0)
    end
    bnew = iTv * v[idxrow]
    # trimming ones
    @inbounds for i in 1:Mred
        abs(1.0 - IG[i,i]) < tol && (IG[i,i] = c1)
    end
    
    G = IG[:, (length(idxd) + 1):end]
    
    Nd, Nf = length(idxd), length(idxf)

    # TODO: check this
    basis = zeros(N, Nf)
    basis[idxd, :] .= -IG[1:Nd, idxf]
    basis[idxf, :] = Matrix(I, Nf, Nf)

    return idxf, idxd, idxmap, G, bnew, basis

end

# ------------------------------------------------------------------
function echelonize2(X, v; tol::Real = 1e-10)

    # @info "echelonize2"

    M, N = size(X)
    @assert M <= N
    Ab = hcat(X, v)
    A, idxd = _rref(Ab; tol)
    Nd = length(idxd)
    A = view(A, 1:Nd, :)
    idxf = setdiff(1:N, idxd)
    Nf = length(idxf)
    G = A[:, idxf]
    be = A[:, end]
    basis = zeros(N, Nf)
    basis[idxd, :] = -A[:, idxf]
    basis[idxf, :] = Matrix(I, Nf, Nf)
    idxmap = vcat(idxd, idxf)

    return idxf, idxd, idxmap, G, be, basis

end

# Derived from https://github.com/blegat/RowEchelon.jl
function _rref(A::AbstractMatrix; tol::Float64=1e-8)

    T = eltype(A)
    nr, nc = size(A)
    idxd = Int[]
    i = j = 1
    @inbounds while i <= nr && j <= nc
        (m, mi) = findmax(abs, view(A, i:nr, j))
        mi = mi+i - 1
        if m <= tol
            if tol > 0
                A[i:nr,j] .= zero(T)
            end
            j += 1
        else
            push!(idxd, j) 
            for k=j:nc
                A[i, k], A[mi, k] = A[mi, k], A[i, k]
            end
            d = A[i,j]
            for k = j:nc
                A[i,k] /= d
            end
            for k = 1:nr
                if k != i
                    d = A[k,j]
                    for l = j:nc
                        A[k,l] -= d*A[i,l]
                    end
                end
            end
            i += 1
            j += 1
        end
    end
    return A, idxd
end

## ------------------------------------------------------------------
let

    # model_id = "iJO1366"
    # model_id = "ecoli_core"
    # model_id = "ECC2"
    model_id = "iJR904"
    net = pull_net(model_id)
    @show size(net)

    tol = 1e-10
    S = Matrix(net.S)
    b = Vector(net.b)

    @info "echelonize1"
    @time begin
        idxf1, idxd1, idxmap1, G1, be1, basis1 = echelonize1(S, b; tol)
    end
    @info "echelonize2"
    @time begin
        idxf2, idxd2, idxmap2, G2, be2, basis2 = echelonize2(S, b; tol)
    end
    @info "MetXBase.echelonize"
    @time begin
        idxf3, idxd3, idxmap3, G3, be3, basis3 = MetXBase.echelonize(S, b; tol)
    end

    @assert G3 == G2
    # plot(spy(basis1), spy(basis2))
    # plot(spy(G1), spy(G2))
    # sum(sort(G1[:]) .- sort(G2[:]))
    # p = plot()
    # plot!(p, sort(G1[:]))
    # plot!(p, sort(G2[:]))
    # plot!(p, sort(basis1[:]))
    # plot!(p, sort(basis2[:]))
    # p
    return nothing
end