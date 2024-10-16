## - - - - - - - - - - - - - - - - - - - - - -
@time begin 
    using Gurobi
    using MetXEP
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using CairoMakie
    using Statistics
    using SparseArrays
    using Base.Threads
    using Distributions: mvnormal_c0
end

## - - - - - - - - - - - - - - - - - - - - - -
include("utils.jl")

## - - - - - - - - - - - - - - - - - - - - - -
# FBA TEST
let
    biom_th = 1e-1
    N = 100
    z1s = zeros(N)
    z2s = zeros(N)
    @threads :static for it in 1:N
        net1 = sample_random_net()
        net2 = _fixxed_net(net1, "BWD_EX_glc__D_e", 10.0)
        z1 = _objective_value(net1)
        z2 = _objective_value(net2)
        z1s[it] = z1
        z2s[it] = z2
        # @show objective_value(lpm)
    end

    # Plots
    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = "net1", ylabel = "net2"
    )
    scatter!(ax, z1s, z2s; label = "z")
    lines!(ax, z1s, z1s; label = "y=x", linestyle = :dash)
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
# Z derivatives
let
    # Random.seed!(123)

    N = 100

    global dZdUs = zeros(BigFloat, N)
    global Z0s = zeros(BigFloat, N)
    global Z1s = zeros(BigFloat, N)
    global davdUs = []
    global av0s = []
    global av0s_Urange = []
    global av1s = []

    global Urange = range(10.0, 9.0; length = 5)

    @threads :static for it in 1:N
        try
            @show it
            # nets
            net0 = sample_random_net(;)
            net1 = _fixxed_net(net0, "BWD_EX_glc__D_e", 10.0)
            ub!(net0, "BWD_EX_glc__D_e", 10.0)
            
            # net0
            # numerical derivatives
            
            _Zs = Float64[]
            _avs_Urange = []
            EP_Urange(net0, "BWD_EX_glc__D_e", Urange) do epm
                # Z = mean(epm, "FWD_BIOMASS_Ecoli_core_w_GAM")
                F = exp(-big(free_energy(epm)[1]))
                push!(_Zs, F)
                _avs = Dict(rxn => mean(epm, rxn) for rxn in colids(epm))
                push!(_avs_Urange, _avs)
                return false
            end
            dZdU = mean(_num_av(_Zs, step(Urange); w = 2))
            _davdUs_dict = Dict()
            for rxn in colids(net0)
                _rxn_avs = get.(_avs_Urange, rxn, NaN)
                _davdUs_dict[rxn] = mean(_num_av(_rxn_avs, step(Urange); w = 2))
            end
            push!(davdUs, _davdUs_dict)
            push!(av0s_Urange, _avs_Urange)
            @show dZdU

            # partition functions and means
            _net0 = fva_strip(net0, OPTIMIZER)
            epm0 = FluxEPModelT0(_net0)
            converge!(epm0)
            Z0 = exp(-big(free_energy(epm0)[1]))
            _av0s = Dict(rxn => mean(epm0, rxn) for rxn in colids(_net0))
            @show Z0
            
            _net1 = fva_strip(net1, OPTIMIZER)
            epm1 = FluxEPModelT0(_net1)
            converge!(epm1)
            Z1 = exp(-big(free_energy(epm1)[1]))
            _av1s = Dict(rxn => mean(epm1, rxn) for rxn in colids(_net1))
            @show Z1
            
            Z0s[it] = Z0
            Z1s[it] = Z1
            dZdUs[it] = dZdU
            push!(av0s, _av0s)
            push!(av1s, _av1s)

        catch e
            @show e
            # rethrow(e)
        end
    end
    # net1
    
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure(;
        size = (1000, 1000)
    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - 
    # dZdU
    _fil = (x) -> !isnan(x) && !isinf(x)
    idxs = intersect(
        findall(_fil, Z0s),
        findall(_fil, Z1s),
        findall(_fil, dZdUs),
    )

    _Z0s = Z0s[idxs]
    _Z1s = Z1s[idxs]
    _dZdUs = dZdUs[idxs]

    ax = Axis(f[1,1:3]; 
        title = "Z", 
        xlabel = "dZdU numerical", 
        ylabel = "dZdU analytic"
    )
    scatter!(ax, log.(abs.(_dZdUs)), log.(abs.(_Z0s)); label = "Z0")
    scatter!(ax, log.(abs.(_dZdUs)), log.(abs.(_Z1s)); label = "Z1")
    xs = log.(abs.(_dZdUs))
    lines!(ax, sort(xs), sort(xs); label = "y=x")
    _corr0 = cor(log.(abs.(_dZdUs)), log.(abs.(_Z0s))) |> Float64
    _corr1 = cor(log.(abs.(_dZdUs)), log.(abs.(_Z1s))) |> Float64
    @show _corr0
    @show _corr1
    axislegend(ax; position = :lt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - 
    # dav_dU
    _rxns = intersect(keys.(davdUs)..., keys.(av0s)..., keys.(av1s)...) |> collect
    rxn = rand(_rxns)
    # rxn = "BWD_FUM"
    ax = Axis(f[1,4:6]; 
        title = rxn, 
        xlabel = "davdU numerical", 
        ylabel = "davdU analytic"
    )
    
    _rxn_davdUs = get.(davdUs, rxn, NaN)
    _rxn_av0s = get.(av0s, rxn, NaN)
    _rxn_av1s = get.(av1s, rxn, NaN)
    # _rxn_davdU_anal = (_Z1s ./ _Z0s) .* (_rxn_av1s .- _rxn_av0s)
    _rxn_davdU_anal = (_rxn_av1s .- _rxn_av0s)
    # _rxn_davdU_anal = (_Z1s ./ _Z0s) 
    # _rxn_davdU_anal = exp.(F0s .- F1s)
    scatter!(ax, _rxn_davdUs, _rxn_davdU_anal)

    # - - - - - - - - - - - - - - - - - - - - - - - - - 
    # avUrange
    ax = Axis(f[1,7:9]; 
        title = rxn, 
        xlabel = "UB", 
        ylabel = "mean"
    )
    for _av_dat in av0s_Urange
        _avs = get.(_av_dat, rxn, NaN)
        lines!(ax, Urange, _avs)
    end

    f
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    _rxns = intersect(keys.(davdUs)..., keys.(av0s)..., keys.(av1s)...) |> collect
    # rxn = _rxns[2]
    # _rxn_davdUs = get.(davdUs, rxn, NaN)
end

## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - -