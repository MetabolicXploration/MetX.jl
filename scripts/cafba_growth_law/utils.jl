## - - - - - - - - - - - - - - - - - - - - - -
# GUROBI
using Gurobi
GRB_ENV = Gurobi.Env()
OPTIMIZER = () -> Gurobi.Optimizer(GRB_ENV)

## - - - - - - - - - - - - - - - - - - - - - -
# Adding cost to net
function _cafba_model(net0;
    ϕmax = 0.5, 
    ignore = (rxn) -> false
)
# posdef
net1 = posdef(net0; ignore)

# Add cost struct
# Add nutrient uptake reations
net2 = resize(net1; 
    nmets = size(net1, 1) + 1,
    nrxns = size(net1, 2) + 1,
)

cost_rxni = findfirst_empty_spot(net2, colids)
@assert !isnothing(cost_rxni)
cost_meti = findfirst_empty_spot(net2, rowids)
@assert !isnothing(cost_meti)

# cost vector
net2.S[cost_meti, :] .= 0
net2.S[cost_meti, cost_rxni] = -1.0
net2.mets[cost_meti] = "COST"
net2.metNames[cost_meti] = "COST DUMMY MET"
net2.lb[cost_rxni] = 0.0
net2.ub[cost_rxni] = ϕmax
net2.rxns[cost_rxni] = "EX_COST"
net2.rxnNames[cost_rxni] = "COST CONTROL RXN"
net2.subSystems[cost_rxni] = "COST"
net2.c[cost_rxni] = 0.0

# Add costs sectors
# for rxn in colids(net2)
#     # IGNORE
#     rxn == "EX_COST" && continue

#     # GROWTH (WR) = 0.169
#     if rxn == "BIOMASS_Ecoli_core_w_GAM" 
#         stoi!(net2, cost_meti, rxn, wr)
#         continue
#     end

#     # INTAKES (WC)
#     # if startswith(rxn, "BACK_EX_")
#     if rxn in wc_iders
#         stoi!(net2, cost_meti, rxn, wc)
#         continue
#     end

#     # INTERNALS (WE) 
#     stoi!(net2, cost_meti, rxn, wE)
# end
return net2
end

## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
function _cafba_w!(net0, ider, w)
    cost_meti = rowindex(net0, "COST")
    stoi!(net0, cost_meti, ider, w)
end

## - - - - - - - - - - - - - - - - - - - - - -
function compute_pred_glcgln_z(glc_z, gln_z; z_c = 1.16)
    pred_glcgln_z = gln_z + glc_z - (2 * (glc_z * gln_z) / z_c)
    pred_glcgln_z = pred_glcgln_z / (1 - ((glc_z * gln_z) / z_c^2))
    return pred_glcgln_z
end

## - - - - - - - - - - - - - - - - - - - - - -
function _plot_heat_map(;
        xs, ys, zs, 
        title = "",
        xlabel = "",
        ylabel = "",
        limits = (nothing, nothing, nothing, nothing),
        dim1_T = (x1) -> identity.(x1),
        dim2_T = (x2) -> identity.(x2),
        cs_T = (cs) -> log10.(cs),
        label = "",
    )

    # Plot
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,1:4]; 
        title, limits, xlabel, ylabel
    )
    x1 = dim1_T(collect(xs)) # koma len
    x2 = dim2_T(collect(ys)) # rxn idx
    w = collect(zs)
    sidx = sortperm(w; rev = true)
    # cs = log10.(w[sidx] ./ maximum(w[sidx]))
    cs = cs_T(w[sidx] ./ maximum(w[sidx]))
    scatter!(ax, x1[sidx], x2[sidx]; 
        colormap = :viridis, 
        markersize = 20, 
        # marker = :square,
        marker = '◼',
        color = cs, 
        alpha = 1.0,
    )
    Colorbar(g[1:3, 5]; 
        label,
        colormap = :viridis, 
        limits = extrema(cs), 
    )
    f
end