@time begin
    import MetNets
    import MetEP
    import MetLP
    import MetXNetHub
    import MetXBase
    import MetXOptim
    import MetXOptim: GLPK
end

## ---------------------------------------------
function to_MetNets(netX::MetXBase.MetNet)
    MetNets.MetNet(;
        S = netX.S, 
        b = netX.b, 
        c = netX.c, 
        lb = netX.lb, 
        ub = netX.ub, 
        mets = netX.mets, 
        rxns = netX.rxns, 
    )
end

function to_MetX(net::MetNets.MetNet)
    MetXBase.MetNet(;
        S = net.S, 
        b = net.b, 
        c = net.c, 
        lb = net.lb, 
        ub = net.ub, 
        mets = net.mets, 
        rxns = net.rxns, 
    )
end

## ---------------------------------------------
let
    global model_id = "ECC2"
    global netX = MetXNetHub.pull_net(model_id)
    global biom_id = MetXBase.extras(netX, "BIOM")
    global glc_id = MetXBase.extras(netX, "EX_GLC")
    global netX = MetXOptim.box(netX, GLPK.Optimizer; protect_obj = true, verbose = true)

    global opm = MetXOptim.fba(netX, GLPK.Optimizer)
    @show MetXOptim.solution(opm, biom_id)
    
    # global net = to_MetNets(netX)
    # global net = MetLP.fva_preprocess(net)
    
    # global netX = to_MetX(net)
    # global opm = MetXOptim.fba(netX, GLPK.Optimizer)
    # @show MetXOptim.solution(opm, biom_id)

end