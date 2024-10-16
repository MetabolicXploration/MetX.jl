@time begin
    using MetX
    using Gurobi
end

## ------------------------------------------------------------------
# ISSUE: Echelon Pivot
## ------------------------------------------------------------------

## ------------------------------------------------------------------
# Implement pivoting. Right now we can not control the set of independent 
# colums which is selected. It can even be the last column, 
# which in a net extended matrix (Sb) is not a flux column.
let
    global net0 = pull_net("ECC2")
    global net = fva_strip(net0, Gurobi.Optimizer, verbose = true)

    # EP
    global enet = EchelonLEPModel(net)
    @show size(net.S)
    @assert maximum(enet.idxd) > size(net.S, 2) # This should fail
    @assert maximum(enet.idxi) > size(net.S, 2)
end





