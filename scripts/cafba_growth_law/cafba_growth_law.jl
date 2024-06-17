## - - - - - - - - - - - - - - - - - - - - - -
@time begin 
    using Makie
    using CairoMakie
    using Gurobi
    using MetXEP
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Statistics
    using Distributions: mvnormal_c0
end

## - - - - - - - - - - - - - - - - - - - - - -
include("utils.jl")

## - - - - - - - - - - - - - - - - - - - - - -
# CAFBA GROWTH
# TODO: FIX THIS
let
    global _wc_glc_s = collect(range(0.0, 3.2; length = 10))
    global _wc_gln_s = collect(range(0.0, 0.9; length = 2))
    
    global wcs = []
    global wc_glc_s = Float64[]
    global wc_gln_s = Float64[]

    global glc_zs = Dict()
    global gln_zs = Dict()
    global glcgln_zs = Dict()

    # nets
    model_id = "ecoli_core" 
    net0 = pull_net(model_id)
    bounds!(net0, "EX_glc__D_e", -1000.0, 0.0)
    bounds!(net0, "EX_gln__L_e", -1000.0, 0.0)
    
    # _cafba_model
    global net2 = _cafba_model(net0; ϕmax = 0.485)

    # Internals
    for rxn in reactions(net2)
        rxn == "EX_COST" && continue
        contains(rxn, "_EX_") && continue
        # @show rxn
        # _cafba_w!(net2, rxn, 8.3e-4)
        _cafba_w!(net2, rxn, 8.0e-4)
    end
    # Biomass
    _cafba_w!(net2, "FWD_BIOMASS_Ecoli_core_w_GAM", 0.169)
    
    for (wc1, wc2) in Iterators.product(_wc_glc_s, _wc_gln_s)

        println("-"^40)
        @show (wc1, wc2)

        # glc
        _cafba_w!(net2, "BWD_EX_glc__D_e", wc1)
        _cafba_w!(net2, "BWD_EX_gln__L_e", 0)
        lpm = FBAOpModel(net2, OPTIMIZER)
        optimize!(lpm)
        glc_z = solution(lpm, "FWD_BIOMASS_Ecoli_core_w_GAM")
        @show glc_z

        # filter
        global z_filter = (z) -> z < 0 ? 1e-3 : z

        push!(wcs, (wc1, wc2))
        push!(wc_glc_s, wc1)
        push!(wc_gln_s, wc2)
        glc_zs[(wc1, wc2)] = z_filter(glc_z)
        # gln_zs[(wc1, wc2)] = z_filter(gln_z)
        # glcgln_zs[(wc1, wc2)] = z_filter(glcgln_z)
        continue

        # gln
        _cafba_w!(net2, "BWD_EX_glc__D_e", 0)
        _cafba_w!(net2, "BWD_EX_gln__L_e", wc2)
        lpm = FBAOpModel(net2, OPTIMIZER)
        optimize!(lpm)
        gln_z = solution(lpm, "FWD_BIOMASS_Ecoli_core_w_GAM")
        @show gln_z

        # glc gln
        _cafba_w!(net2, "BWD_EX_glc__D_e", wc1)
        _cafba_w!(net2, "BWD_EX_gln__L_e", wc2)
        lpm = FBAOpModel(net2, OPTIMIZER)
        optimize!(lpm)
        glcgln_z = solution(lpm, "FWD_BIOMASS_Ecoli_core_w_GAM")
        @show glcgln_z

        # filter
        global z_filter = (z) -> z < 0 ? 1e-3 : z

        push!(wcs, (wc1, wc2))
        push!(wc_glc_s, wc1)
        push!(wc_gln_s, wc2)
        glc_zs[(wc1, wc2)] = z_filter(glc_z)
        gln_zs[(wc1, wc2)] = z_filter(gln_z)
        glcgln_zs[(wc1, wc2)] = z_filter(glcgln_z)
    end
    
end

## - - - - - - - - - - - - - - - - - - - - - -
# PREDICTED GROWTH
let
    z_c = 0.9
    # z_c = 2.2
    global pred_glcgln_zs = Dict()
    for wp in wcs
        glc_z = glc_zs[wp]
        gln_z = gln_zs[wp]
        pred_glcgln_z = compute_pred_glcgln_z(glc_z, gln_z; z_c)
        pred_glcgln_zs[wp] = z_filter(pred_glcgln_z)
    end
end

## - - - - - - - - - - - - - - - - - - - - - -
# PLOTS
## - - - - - - - - - - - - - - - - - - - - - -
# Correlation
let
    
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,1:4]; 
        title = "ecoli_core", 
        limits = (nothing, nothing, nothing, nothing), 
        xlabel = "ws index", 
        ylabel = "λ"
    )

    
    xs = Float64[pred_glcgln_zs[wc] for wc in wcs]
    ridx = sortperm(xs)
    # lines!(ax, xs[ridx]; label = "λ12 heuristic")
    
    xs = Float64[glcgln_zs[wc] for wc in wcs]
    # lines!(ax, xs[ridx]; label = "λ12 cafba")
    
    xs = Float64[glc_zs[wc] for wc in wcs]
    lines!(ax, xs[ridx]; label = "λglc cafba")

    xs = Float64[gln_zs[wc] for wc in wcs]
    lines!(ax, xs[ridx]; label = "λgln cafba")
    
    axislegend(ax; position = :lb)
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
# model growth heat map
# wc vs wc, 
let
    xs = wc_glc_s
    ys = wc_gln_s
    zs = map(wcs) do wp
        z = glcgln_zs[wp]
        return z < 0 ? 1e-3 : z
    end

    _plot_heat_map(;
        xs, ys, zs, 
        title = "ecoli_core",
        xlabel = "wc_glc",
        ylabel = "wc_gln",
        limits = (nothing, nothing, nothing, nothing),
        dim1_T = (x1) -> identity.(x1),
        dim2_T = (x2) -> identity.(x2),
        cs_T = (cs) -> log10.(cs),
        label = "log λ12",
    )
end



## - - - - - - - - - - - - - - - - - - - - - -
let
    xs = wc_glc_s
    ys = wc_gln_s
    zs = map(wcs) do wp
        z = glcgln_zs[wp]
        return z < 0 ? 1e-3 : z
    end

    _plot_heat_map(;
        xs, ys, zs, 
        title = "ecoli_core",
        xlabel = "wc_glc",
        ylabel = "wc_gln",
        limits = (nothing, nothing, nothing, nothing),
        dim1_T = (x1) -> identity.(x1),
        dim2_T = (x2) -> identity.(x2),
        cs_T = (cs) -> log10.(cs),
        label = "log λ12",
    ) 
end
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "wc [gh/mmol]", ylabel = "log λ [1/h]")

    z_c = 1.16
    # z_c = 2.2
    pred_glcgln_zs = gln_zs .+ glc_zs .- (2 .* (glc_zs .* gln_zs) ./ z_c)
    pred_glcgln_zs = pred_glcgln_zs ./ (1 .- ((glc_zs .* gln_zs) ./ z_c^2))

    idx_ = 1:10
    Tx = v -> log.(v .+ 1e-3)
    Ty = v -> v
    lines!(ax, Tx(wcs), Ty(gln_zs); label = "λgln cafba",
        linewidth = 4, linestyle = :dot
    )
    lines!(ax, Tx(wcs), Ty(glc_zs); label = "λglc cafba", 
        linewidth = 4, linestyle = :dot
    )
    lines!(ax, Tx(wcs), Ty(glcgln_zs); label = "λ12 cafba", 
        linewidth = 4, linestyle = :solid
    )
    lines!(ax, Tx(wcs), Ty(pred_glcgln_zs); label = "λ12 phen", 
        linewidth = 4, linestyle = :solid
    )
    axislegend(ax; position = :rt)
    f
end