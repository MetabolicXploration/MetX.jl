## - - - - - - - - - - - - - - - - - - - - - -
# GUROBI
using Gurobi
GRB_ENV = Gurobi.Env()
OPTIMIZER = () -> Gurobi.Optimizer(GRB_ENV)

## - - - - - - - - - - - - - - - - - - - - - -
function _num_av(ys, dx; w = 5)
    dydx = zeros(eltype(ys), length(ys))
    for i in eachindex(ys)
        i0 = max(firstindex(ys), i - w)
        i1 = min(lastindex(ys), i + w)
        dydx[i] = mean(diff(ys[i0:i1])) ./ dx
    end
    return dydx
end

## - - - - - - - - - - - - - - - - - - - - - -
function _fixxed_net(net, ider, b)
    net = deepcopy(net)
    ider = rxnindex(net, ider)
    lb!(net, ider, b)
    ub!(net, ider, b)
    empty_fixxed!(net)
    empty_void_iders!(net)
    return emptyless_model(net)
end

## - - - - - - - - - - - - - - - - - - - - - -
function random_net(; ignore = nothing, margin = 1e-1)

    model_id = "ecoli_core"
    net0 = pull_net(model_id)
    m, n = size(net0)

    # randomlly change S columns
    if isnothing(ignore)
        ignore = (rxn) -> begin
            contains(rxn, "BIOMASS_Ecoli_core_w_GAM") && return true
            contains(rxn, "EX_") && return true
            return false
        end
    end
    for (rxni, rxn) in enumerate(colids(net0))
        ignore(rxn) === true && continue
        net0.S[:, rxni] .*= rand(m) .* (1.0 - margin) .+ margin
    end

    # posdef
    net1 = posdef(net0)
    
    return net1
end

## - - - - - - - - - - - - - - - - - - - - - -
function EP_Urange(f::Function, net::MetNet, ider, Urange)
    u_bk = ub(net, ider)
    for u in Urange
        ub!(net, ider, u)
        _net = box(net, OPTIMIZER)
        epm = FluxEPModelT0(_net)
        converge!(epm)
        f(epm) === true && break
    end
    ub!(net, ider, u_bk)
    return nothing
end

## - - - - - - - - - - - - - - - - - - - - - -
function LP_Urange(f::Function, net::MetNet, ider, Urange)
    lpm = FBAOpModel(net, OPTIMIZER)
    for u in Urange
        ub!(lpm, ider, u)
        optimize!(lpm)
        f(lpm) === true && break
    end
    return nothing
end

## - - - - - - - - - - - - - - - - - - - - - -
function sample_random_net(; 
        ignore = nothing, 
        margin = 1e-1, 
        biom_th = 1e-1, 
        dzdU_th = 0.1, 
        Uider = "BWD_EX_glc__D_e", 
        Urange = range(9.0, 10.0; length = 5)
    )

    
    for it in 1:100
        try
            # net
            net0 = random_net(; ignore, margin)
            
            # z filter
            z = _objective_value(net0)
            z > biom_th || continue

            # dzdU filter
            zs = Float64[]
            LP_Urange(net0, Uider, Urange) do lpm
                z = objective_value(lpm)
                push!(zs, z)
            end
            dzdU = mean(_num_av(zs, step(Urange); w = 2))
            dzdU > dzdU_th || continue

            return net0
        catch e; end
    end
    error("Sampling failed!")
end

## - - - - - - - - - - - - - - - - - - - - - -
function _objective_value(net)
    lpm = FBAOpModel(net, OPTIMIZER)
    optimize!(lpm)
    return objective_value(lpm)
end

## - - - - - - - - - - - - - - - - - - - - - -

nothing
