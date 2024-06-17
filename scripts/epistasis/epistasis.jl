## - - - - - - - - - - - - - - - - - - - - - -
@time begin 
    using Gurobi
    using MetX
    using CairoMakie
    using DataFrames
    using Base.Iterators: product, flatten
end

## - - - - - - - - - - - - - - - - - - - - - -
include("0_utils.jl")

## - - - - - - - - - - - - - - - - - - - - - -
let
    # nets
    global net0 = pull_net("ecoli_core")
    global net1 = posdef(net0)
    global biom_id = "FWD_BIOMASS_Ecoli_core_w_GAM"

    # medium
    global exch1 = "BWD_EX_glc__D_e"
    ub!(net1, exch1, 10.0)
    global exch2 = "BWD_EX_lac__D_e"
    ub!(net1, exch2, 10.0)

    # EP_Urange
    global epm_12pool = Dict();
    global Urange1 = range(9.0, 10.0; length = 25)
    global Urange2 = range(9.0, 10.0; length = 25)
    box_kwargs = (; nths = 1, verbose = false)

    @time EP_Urange(net1, exch1, Urange1, exch2, Urange2; box_kwargs) do _epm
        u1, u2 = extras(_epm, :EP_Urange)
        @show u1, u2
        @show mean(_epm, biom_id)
        epm_12pool[(u1, u2)] = _epm
    end

    return nothing
end

## - - - - - - - - - - - - - - - - - - - - - -
# mean_12mats
let
    global mean_12mats = Dict()
    for rxn in reactions(net1)
        _mat = mean_12mats[rxn] = zeros(length(Urange1), length(Urange2))
        for (i, u1) in enumerate(Urange1)
            for (j, u2) in enumerate(Urange2)
                println("- "^20)
                @show rxn
                @show u1, u2
                epm = epm_12pool[(u1, u2)]
                av = hascolid(epm, rxn) ? mean(epm, rxn) : NaN
                @show av
                _mat[i,j] = av
            end
        end
    end
end

## - - - - - - - - - - - - - - - - - - - - - -
# plots
## - - - - - - - - - - - - - - - - - - - - - -
let
    # rxn = rand(keys(mean_12mats))
    rxn = "FWD_BIOMASS_Ecoli_core_w_GAM"
    mat = mean_12mats[rxn]
    
    xs = _flatten_product(Urange1, Urange2; dim = 1)
    ys = _flatten_product(Urange1, Urange2; dim = 2)
    zs = _flatten_product(eachindex(Urange1), eachindex(Urange2)) do li, t
        return mat[t...]
    end

    # Plot
    title = string("ecoli_core", "\n", rxn)
    xlabel = L"U_1"
    ylabel = L"U_2"
    label = "av"
    dim1_T = (x1) -> identity.(x1)
    dim2_T = (x2) -> identity.(x2)
    cs_T = (cs) -> (cs)
    limits = (nothing, nothing, nothing, nothing)

    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,1:4]; 
        title, limits, xlabel, ylabel
    )
    x1 = dim1_T(collect(xs)) #
    x2 = dim2_T(collect(ys)) #
    w = collect(zs)
    sidx = sortperm(w; rev = true)
    # cs = log10.(w[sidx] ./ maximum(w[sidx]))
    # cs = cs_T(w[sidx] ./ maximum(w[sidx]))
    cs = w[sidx]
    scatter!(ax, x1[sidx], x2[sidx]; 
        colormap = :viridis, 
        markersize = 20, 
        marker = '◼',
        color = cs, 
        alpha = 1.0,
    )
    Colorbar(g[1:3, 5]; 
        colormap = :viridis, limits = extrema(cs), label
    )
    f
end


## - - - - - - - - - - - - - - - - - - - - - -
# derivatives
## - - - - - - - - - - - - - - - - - - - - - -
function _mixed_derivatives(mat, i1, i2, i1step, i2step)
    ((i1 == 1) || (i2 == 1)) && return 0.0
    v12 = mat[i1 - 1, i2 - 1]
    v2 = mat[i1, i2 - 1]
    v1 = mat[i1 - 1, i2]
    v0 = mat[i1, i2]
    dv12 = ((v12 - v2) - (v1 - v0)) / (i1step * i2step)
    return dv12
end

## - - - - - - - - - - - - - - - - - - - - - -
# ddv_dU1U2_mats
let
    global ddv_dU1U2_mats = Dict()
    for rxn in reactions(net1)
        _mat = ddv_dU1U2_mats[rxn] = zeros(length(Urange1), length(Urange2))
        for (i1, u1) in enumerate(Urange1)
            for (i2, u2) in enumerate(Urange2)
                println("- "^20)
                @show rxn
                @show u1, u2
                ddv_dU1U2 = _mixed_derivatives(
                    mean_12mats[rxn], i1, i2, 
                    step(Urange1), step(Urange2)
                )
                @show ddv_dU1U2
                _mat[i1,i2] = ddv_dU1U2
            end
        end
    end
end

## - - - - - - - - - - - - - - - - - - - - - -
# plots
## - - - - - - - - - - - - - - - - - - - - - -
let
    # rxn = rand(keys(mean_12mats))
    rxn = "FWD_BIOMASS_Ecoli_core_w_GAM"
    @show rxn
    mat = ddv_dU1U2_mats[rxn]
    
    xs = _flatten_product(Urange1, Urange2; dim = 1)
    ys = _flatten_product(Urange1, Urange2; dim = 2)
    zs = _flatten_product(eachindex(Urange1), eachindex(Urange2)) do li, t
        return mat[t...]
    end

    # Plot
    title = string("ecoli_core", "\n", rxn)
    xlabel = L"U_1"
    ylabel = L"U_2"
    label = "av"
    dim1_T = (x1) -> identity.(x1)
    dim2_T = (x2) -> identity.(x2)
    cs_T = (cs) -> (cs)
    limits = (nothing, nothing, nothing, nothing)

    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,1:4]; 
        title, limits, xlabel, ylabel
    )
    x1 = dim1_T(collect(xs)) #
    x2 = dim2_T(collect(ys)) #
    w = collect(zs)
    sidx = sortperm(w; rev = true)
    # cs = log10.(w[sidx] ./ maximum(w[sidx]))
    # cs = cs_T(w[sidx] ./ maximum(w[sidx]))
    cs = w[sidx]
    scatter!(ax, x1[sidx], x2[sidx]; 
        colormap = :viridis, 
        markersize = 20, 
        marker = '◼',
        color = cs, 
        alpha = 1.0,
    )
    Colorbar(g[1:3, 5]; 
        colormap = :viridis, limits = extrema(cs), label
    )
    f
end
