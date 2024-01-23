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
end

## - - - - - - - - - - - - - - - - - - - - - -
include("utils.jl")

## - - - - - - - - - - - - - - - - - - - - - -
let 
    # net0
    global model_id = "ecoli_core"
    global net0 = pull_net(model_id)
    lb!(net0, "EX_glc__D_e", 0)
    lb!(net0, "EX_gln__L_e", -10)
    
    # posdef
    global net1 = posdef(net0)

    lim_id = "BACK_EX_gln__L_e"
    # lim_id = "BACK_EX_glc__D_e"

    global ubs = 1.0:0.1:10.0

    global F0s = Float64[]
    global S0s = Float64[] 
    global av0s = []
    for u in ubs
        lb!(net1, lim_id, 0)
        ub!(net1, lim_id, u)
        global lep = box(net1, GLPK.Optimizer; verbose = false)
        global epm0 = FluxEPModelT0(lep)
        config!(epm0; verbose = false)
        converge!(epm0)
        av = Dict(rxn => mean(epm0, rxn) for rxn in colids(epm0))
        F = free_energy(epm0)[1]
        S = entropy(epm0)
        push!(S0s, S)
        push!(F0s, F)
        push!(av0s, av)
    end

    global F1s = Float64[]
    global S1s = Float64[]
    global av1s = []
    for u in ubs
        net2 = _fixxed_net(net1, lim_id, u)
        global lep = box(net2, GLPK.Optimizer; verbose = false)
        global epm0 = FluxEPModelT0(lep)
        config!(epm0; verbose = false)
        converge!(epm0)
        av = Dict(rxn => mean(epm0, rxn) for rxn in colids(epm0))
        F = free_energy(epm0)[1]
        S = entropy(epm0)
        push!(S1s, S)
        push!(F1s, F)
        push!(av1s, av)
    end
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "UB", title = "Z")
    scatter!(ax, ubs, F0s; label = "F0")
    scatter!(ax, ubs, F1s; label = "F1")
    # scatter!(ax, ubs, S0s; label = "S0", color = :red)
    # scatter!(ax, ubs, S1s; label = "S1", color = :blue)
    axislegend(ax; position = :lt)
    f
end


## - - - - - - - - - - - - - - - - - - - - - -
# dZ/dU
let
    global Z0s = exp.(-big.(F0s))
    global Z1s = exp.(-big.(F1s))
    # dZ0s = diff(Z0s) ./ step(ubs)
    # push!(dZ0s, last(dZ0s))
    dZ0s = _num_av(Z0s, step(ubs); w = 1)
    
    f = Figure()

    ax = Axis(f[1,1]; 
        xlabel = L"\log~~\frac{\partial Z}{\partial U_i}", 
        ylabel = L"\log~~Z", 
    )
    lines!(ax, log.(abs.(dZ0s)), log.(abs.(dZ0s)); label = L"y=x", color = :black)
    scatter!(ax, log.(abs.(dZ0s)), -F0s; label = L"Z", color = :red)
    scatter!(ax, log.(abs.(dZ0s)), -F1s; label = L"Z(S_{-i}, b(U_{i}), c_{-i})", color = :blue)
    axislegend(ax; position = :lt)

    ax = Axis(f[1,2]; 
        xlabel = L"UB", 
        ylabel = L"\log~~Z", 
    )
    # scatter!(ax, ubs, log.(Z0s); label = L"Z", color = :red)
    scatter!(ax, ubs, -F0s; label = L"Z", color = :red)
    # scatter!(ax, ubs, log.(Z1s); label = L"Z(S_{-i}, b(U_{i}), c_{-i})", color = :blue)
    scatter!(ax, ubs, -F1s; label = L"Z(S_{-i}, b(U_{i}), c_{-i})", color = :blue)
    # axislegend(ax; position = :lt)

    # ax = Axis(f[1,3]; 
    #     xlabel = "Ub", 
    #     ylabel = "ratio", 
    # )
    # scatter!(ax, ubs, log.(dZ0s./Z1s); label = L"Z", color = :red)
    # lines!(ax, ubs, fill(log(1), length(ubs)); label = L"Z", color = :red)
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    
    # biomass
    ider = "BIOMASS_Ecoli_core_w_GAM"
    # ider = colids(net1, rand(1:132))
    _av0s = getindex.(av0s, ider)
    _av1s = getindex.(av1s, ider)

    f = Figure(;)

    ax = Axis(f[1,2]; 
        title = ider,
        xlabel = L"UB", 
        ylabel = L"\bar{v}", 
    )
    scatter!(ax, ubs, _av0s; label = L"Z0", color = :red)
    scatter!(ax, ubs, _av1s; label = L"Z1", color = :blue)
    axislegend(ax; position = :lt)
    
    dav_num = _num_av(_av0s, step(ubs); w = 1)
    # Z0s = exp.(big.(S0s))
    # Z1s = exp.(big.(S1s))
    # dav_anal = (Z1s ./ Z0s) .* (_av1s .- _av0s)
    dav_anal = exp.(F0s .- F1s) .* (_av1s .- _av0s)
    # dav_anal = (_av1s .- _av0s)
    # dav_anal = (Z1s ./ Z0s) 
    # dav_anal = exp.(F0s .- F1s)

    ax = Axis(f[1,1]; 
        xlabel = "numeric", 
        ylabel = "analytic", 
    )
    scatter!(ax, log.(dav_num), log.(dav_anal), color = :red)

    ax = Axis(f[1,3]; 
        xlabel = "UB", 
        ylabel = L"\frac{\partial \bar{v}}{\partial U_{i}}", 
    )
    # scatter!(ax, ubs, dav_anal, color = :red)
    scatter!(ax, ubs, dav_num, color = :red)

    
    return f
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "LB", title = "averages")
    for idx in 1:132
        _av0s = getindex.(av0s, idx)
        # _av0s .-= minimum(_av0s)
        _av1s = getindex.(av1s, idx)
        # _av1s .-= minimum(_av1s)
        # lines!(ax, ubs, _av0s; label = "av0", color = :red)
        # lines!(ax, ubs, _av1s; label = "av1", color = :blue)
        lines!(ax, ubs, _av1s - _av0s; label = "av1", color = :blue)
    end
    # axislegend(ax; position = :lt)
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "numerical", ylabel = "analytic")

    global corrs = Float64[]
    for idx in 1:132
        _av0s = getindex.(av0s, idx)
        _av1s = getindex.(av1s, idx)
        dav0s = diff(_av0s) ./ step(ubs)
        dav1s = exp.(F0s .- F1s) .* ( _av1s .- _av0s )
        dav1s = dav1s[2:end]
        # scatter!(ax, log10.(_av0s), log10.(_av1s))
        scatter!(ax, dav0s, dav1s)
        push!(corrs, cor(dav1s, dav0s))
    end
    # axislegend(ax; position = :lt)
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "rxn index (sorted)", ylabel = "corr index")
    scatter!(ax, sort(corrs))
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
let

end