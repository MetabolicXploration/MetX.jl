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
function _cafba_w!(net0, ider, w)
cost_meti = rowindex(net0, "COST")
stoi!(net0, cost_meti, ider, w)
end
