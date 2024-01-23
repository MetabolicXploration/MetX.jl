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
function _build_box(m, n)
    return MetNet(;
        S = zeros(m, n),
        b = zeros(m),
        lb = zeros(n),
        ub = ones(n),
        rxns = ["rxn$i" for i in 1:n],
        mets = ["met$i" for i in 1:m]
    )
end

## - - - - - - - - - - - - - - - - - - - - - -
# ecoli
let
    global Fs = Float64[]
    global Ss = Float64[]
    global S2s = Float64[]
    global lbs = -10.0:0.1:-1.0
    global net0 = pull_net("ecoli_core")
    for l in lbs
        lb!(net0, "EX_glc__D_e", l)
        net1 = box(net0, GLPK.Optimizer)
        epm0 = FluxEPModelT0(net1)
        converge!(epm0)
        S2 = entropy(MultivariateNormal(epm0))
        S = entropy(epm0)
        F = free_energy(epm0)[1]
        # VF = log(-F) / size(net, 2)
        # VS = log(S)
        push!(Fs, F)
        push!(Ss, S)
        push!(S2s, S)
    end
end

## - - - - - - - - - - - - - - - - - - - - - -
# Plots
let
    # global VFs = log.(-Fs)
    # global VSs = log.(Ss)
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "bound", ylabel = "vol")
    # scatter!(ax, lbs, VFs; label = "VF")
    # scatter!(ax, ns, VSs; label = "VS")
    scatter!(ax, lbs, Ss; label = "S")
    scatter!(ax, lbs, S2s; label = "S1")
    axislegend(ax; position = :lb)
    ax = Axis(f[1,2]; xlabel = "bound", ylabel = "vol")
    scatter!(ax, lbs, -Fs; label = "F")
    f
end

## - - - - - - - - - - - - - - - - - - - - - -
# TODO: make it work
let
    global Fs = Float64[]
    global Ss = Float64[]
    global ns = 1:1:100
    for n in ns
        @show n
        net = pull_net("ecoli_core")
        net1 = box(net, GLPK.Optimizer)
        # net = _build_box(n, n)
        epm0 = FluxEPModelT0(net)
        converge!(epm0)
        S = entropy(epm0)
        F = free_energy(epm0)[1]
        # VF = log(-F) / size(net, 2)
        # VS = log(S)
        push!(Fs, F)
        push!(Ss, S)
    end
end