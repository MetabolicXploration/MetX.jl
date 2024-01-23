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
include("utils.jl")

## - - - - - - - - - - - - - - - - - - - - - -
let

    global _wc_glc_s = collect(range(0.0, 1.5; length = 50))
    global _wc_gln_s = collect(range(0.0, 1.2; length = 50))
    
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
        _cafba_w!(net2, rxn, 8.3e-4)
    end
    # Biomass
    _cafba_w!(net2, "FWD_BIOMASS_Ecoli_core_w_GAM", 0.169)
    
    for (wc1, wc2) in Iterators.product(_wc_glc_s, _wc_gln_s)

        println("-"^40)
        @show (wc1, wc2)

        # glc
        _cafba_w!(net2, "BWD_EX_glc__D_e", wc1)
        _cafba_w!(net2, "BWD_EX_gln__L_e", 0)
        lpm = FBAOpModel(net2, GLPK.Optimizer)
        optimize!(lpm)
        glc_z = solution(lpm, "FWD_BIOMASS_Ecoli_core_w_GAM")
        @show glc_z

        # gln
        _cafba_w!(net2, "BWD_EX_glc__D_e", 0)
        _cafba_w!(net2, "BWD_EX_gln__L_e", wc2)
        lpm = FBAOpModel(net2, GLPK.Optimizer)
        optimize!(lpm)
        gln_z = solution(lpm, "FWD_BIOMASS_Ecoli_core_w_GAM")
        @show gln_z

        # glc gln
        _cafba_w!(net2, "BWD_EX_glc__D_e", wc1)
        _cafba_w!(net2, "BWD_EX_gln__L_e", wc2)
        lpm = FBAOpModel(net2, GLPK.Optimizer)
        optimize!(lpm)
        glcgln_z = solution(lpm, "FWD_BIOMASS_Ecoli_core_w_GAM")
        @show glcgln_z

        push!(wcs, (wc1, wc2))
        push!(wc_glc_s, wc1)
        push!(wc_gln_s, wc2)
        glc_zs[(wc1, wc2)] = glc_z
        gln_zs[(wc1, wc2)] = gln_z
        glcgln_zs[(wc1, wc2)] = glcgln_z
    end
    
end

## - - - - - - - - - - - - - - - - - - - - - -
# wc vs wc
let
    # 
    xs = wc_glc_s
    ys = wc_gln_s
    zs = map(wcs) do wp
        z = glcgln_zs[wp]
        return z < 0 ? 1e-3 : z
    end
    sidxs = sortperm(zs; rev = true)

    # Plot
    title = "ecoli_core"
    xlabel = "wc_glc"
    ylabel = "wc_gln"
    dim1_T = (x1) -> identity.(x1)
    dim2_T = (x2) -> identity.(x2)
    cs_T = (cs) -> log10.(cs)
    limits = (nothing, nothing, nothing, nothing)
    dim1_bar_width = 1.0
    dim2_bar_width = 1.0
    colgap = -20
    rowgap = -5

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
        # label = "log10 λ12",
        label = "log10 λ12",
        colormap = :viridis, limits = extrema(cs), 
    )
    f
end



## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "wc [gh/mmol]", ylabel = "log λ [1/h]")

    z_c = 1.16
    # z_c = 2.2
    pred_glcgln_zs = gln_zs .+ glc_zs .- (2 .* (glc_zs .* gln_zs) / z_c)
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