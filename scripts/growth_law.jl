## - - - - - - - - - - - - - - - - - - - - - -
@time begin 
    using GLPK
    using MetXEP
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using CairoMakie
    using Statistics
    using Distributions: mvnormal_c0
    using GLPK
end

## - - - - - - - - - - - - - - - - - - - - - -
# Adding cost to net
function _cofba_model(net0; 
        wr = 0.169, 
        wc = 0.1,
        wE = 8.3e-4,
        ϕmax = 0.5
    )
    # posdef
    net1 = posdef(net0)

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
    for rxn in colids(net2)
        # IGNORE
        rxn == "EX_COST" && continue

        # GROWTH (WR) = 0.169
        if rxn == "BIOMASS_Ecoli_core_w_GAM" 
            stoi!(net2, cost_meti, rxn, wr)
            continue
        end

        # INTAKES (WC)
        if startswith(rxn, "BACK_EX_")
            stoi!(net2, cost_meti, rxn, wc)
            continue
        end

        # INTERNALS (WE) 
        stoi!(net2, cost_meti, rxn, wE)
    end

    # FBA TEST
    # biom = "BIOMASS_Ecoli_core_w_GAM"
    # lpm = FBAOpModel(net2, GLPK.Optimizer)
    # optimize!(lpm)
    # solution(lpm, biom)

    return net2

end

## - - - - - - - - - - - - - - - - - - - - - -
let

    global _wcs = range(0.0, 0.10; length = 100)

    global glc_zs = Float64[]
    global gln_zs = Float64[]
    global glcgln_zs = Float64[]

    # nets
    model_id = "ecoli_core" 
    net0 = pull_net(model_id)
    lb!(net0, "EX_glc__D_e", -1000.0)
    ub!(net0, "EX_glc__D_e", 0.0)
    lb!(net0, "EX_gln__L_e", -1000.0)
    ub!(net0, "EX_gln__L_e", 0.0)

    for wc in _wcs

         # REV_EX_glc__D_e
        global net2 = _cofba_model(net0; 
            wr = 0.169, 
            wc = wc,
            # wE = 8.3e-4,
            wE = 6e-4,
            ϕmax = 0.485
        )

        println("-"^40)
        @show wc

        # glc
        ub!(net2, "REV_EX_glc__D_e", 1000.0)
        ub!(net2, "REV_EX_gln__L_e", 0.0)
        lpm = FBAOpModel(net2, GLPK.Optimizer)
        optimize!(lpm)
        glc_z = solution(lpm, "BIOMASS_Ecoli_core_w_GAM")
        @show glc_z
        push!(glc_zs, glc_z)

        # gln
        ub!(net2, "REV_EX_glc__D_e", 0.0)
        ub!(net2, "REV_EX_gln__L_e", 1000.0)
        lpm = FBAOpModel(net2, GLPK.Optimizer)
        optimize!(lpm)
        gln_z = solution(lpm, "BIOMASS_Ecoli_core_w_GAM")
        @show gln_z
        push!(gln_zs, gln_z)

        # glc gln
        ub!(net2, "REV_EX_glc__D_e", 1000.0)
        ub!(net2, "REV_EX_gln__L_e", 1000.0)
        lpm = FBAOpModel(net2, GLPK.Optimizer)
        optimize!(lpm)
        glcgln_z = solution(lpm, "BIOMASS_Ecoli_core_w_GAM")
        @show glcgln_z
        push!(glcgln_zs, glcgln_z)
    end

    # FBA TEST
    
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "wc [gh/mmol]", ylabel = "λ [1/h]")

    z_c = 1.16
    pred_glcgln_zs = gln_zs .+ glc_zs .- (2 .* (glc_zs .* gln_zs) / z_c)
    pred_glcgln_zs = pred_glcgln_zs ./ (1 .- ((glc_zs .* gln_zs) ./ z_c^2))

    lines!(ax, _wcs, gln_zs; label = "λgln cafba",
        linewidth = 4, linestyle = :dot
    )
    lines!(ax, _wcs, glc_zs; label = "λglc cafba", 
        linewidth = 4, linestyle = :dot
    )
    lines!(ax, _wcs, glcgln_zs; label = "λ12 cafba", 
        linewidth = 4, linestyle = :solid
    )
    lines!(ax, _wcs, pred_glcgln_zs; label = "λ12 phen", 
        linewidth = 4, linestyle = :solid
    )
    axislegend(ax; position = :rt)
    f
end