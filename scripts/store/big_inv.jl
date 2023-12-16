## ------------------------------------------------------------------
@time begin
    using LinearAlgebra
    # using GenericLinearAlgebra
end

## ------------------------------------------------------------------
function _inplaceinverse!(dest::AbstractArray, source::AbstractArray; δ = 1e-10)
    copyto!(dest, source)
    try
        inv!(cholesky!(Hermitian(dest)))
    catch err
        if err isa PosDefException
            nearPD!(dest, δ)
            # @show isposdef(dest)
            inv!(cholesky!(Hermitian(dest)))
            return 
        end
        rethrow(err)
    end
end
## ------------------------------------------------------------------
let 
    S = rand(BigFloat, 10, 10)
    # S = rand(Float64, 10, 10)
    S = S'S


    # S = Hermitian(S'S)
    LinearAlgebra.inv!(cholesky!(Hermitian(S)))
end