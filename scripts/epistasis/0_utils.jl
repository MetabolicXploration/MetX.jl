## - - - - - - - - - - - - - - - - - - - - - -
# GUROBI
using Gurobi
GRB_ENV = Gurobi.Env()
OPTIMIZER = () -> Gurobi.Optimizer(GRB_ENV)

## - - - - - - - - - - - - - - - - - - - - - -
function EP_Urange(f::Function, net::MetNet, ider, Urange)
    u_bk = ub(net, ider)
    for u in Urange
        ub!(net, ider, u)
        _net = box(net, OPTIMIZER)
        epm = FluxEPModelT0(_net)
        extras!(epm, :EP_Urange, (;u))
        converge!(epm)
        f(epm) === true && break
    end
    ub!(net, ider, u_bk)
    return nothing
end

## - - - - - - - - - - - - - - - - - - - - - -
function EP_Urange(
        f::Function, net::MetNet, 
        ider1, Urange1, ider2, Urange2;
        box_kwargs = (;), 
        # ep_kwargs = (;), 
    )
    u1_bk = ub(net, ider1)
    u2_bk = ub(net, ider2)
    for u1 in Urange1
        for u2 in Urange2
            ub!(net, ider1, u1)
            ub!(net, ider2, u2)
            _net = box(net, OPTIMIZER; box_kwargs...)
            epm = FluxEPModelT0(_net)
            extras!(epm, :EP_Urange, (;u1, u2))
            converge!(epm)
            f(epm) === true && break
        end
    end
    ub!(net, ider1, u1_bk)
    ub!(net, ider2, u2_bk)
    return nothing
end

## - - - - - - - - - - - - - - - - - - - - - -
function _flatten_product(f::Function, i0, i1)
    _vec = zeros(length(i0) * length(i1))
    li = 1
    for t in Iterators.product(i0, i1)
        _vec[li] = f(li, t)
        li += 1   
    end
    return _vec
end

function _flatten_product(i0, i1; dim = 1)
    return _flatten_product(i0, i1) do _, t
        return t[dim]
    end
end


## - - - - - - - - - - - - - - - - - - - - - -
nothing