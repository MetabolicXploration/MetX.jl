@time begin
    using MetX
    using Gurobi
end

## ------------------------------------------------------------------
# ISSUE: PosDefException
## ------------------------------------------------------------------

## ------------------------------------------------------------------
# Either because the real ill-condition state of the metabolic network, 
# because 'fva_strip' is not reducing accordingly the network, 
# or EP is accumulating erros, the Cov matrix fail the isposdef test beyong repare.
let
    net0 = pull_net("iJR904")
    net = fva_strip(net0, Gurobi.Optimizer, verbose = true)

    # EP
    epm = FluxEPModelT0(net)
    config!(epm; 
        epsconv = 1e-6, 
        maxiter = Int(1e4)
    )
    try; converge!(epm)
        catch err; @show err
    end
    @show maximum(epm.Σi)
    @show isposdef(epm.Σi)
end    

## ------------------------------------------------------------------
# Artifitially reducing the ill conditioning avoid the error
let
    net0 = pull_net("iJR904")
    clampbounds!(net0, -100.0, 100.0)
    for i in eachindex(net0.S)
        abs(net0.S[i]) < 1e-4 || continue
        net0.S[i] = 0.0
    end
    net = fva_strip(net0, Gurobi.Optimizer, verbose = true)

    # EP
    epm = FluxEPModelT0(net)
    config!(epm; 
        epsconv = 1e-6, 
        maxiter = Int(1e3)
    )
    try; converge!(epm)
        catch err; @show err
    end
    @show maximum(epm.Σi)
    @show isposdef(epm.Σi)
end    

