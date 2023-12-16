@time begin
    using MetX
    using MetXEP
    import MetXEP: erf
    import MetXGrids
    using Plots
    using Gurobi
    using Distributions
end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
let 
    setprecision(BigFloat, 2*16)

    # Net
    global model_id = "ecoli_core"
    global net0 = pull_net("ecoli_core")
    global net = box(net0, Gurobi.Optimizer)

    # EP
    global epm = FluxEPModelT0(net)
    biom_id = extras(net, "BIOM")

    βs = [
        0.0
        10.0.^range(0.0, 6.0; length = 60);
    ]

    Fs, Ss = Float64[], Float64[]
    for β in βs

        println("-"^50)
        @show β

        beta!(epm, biom_id, β)
        converge!(epm)
        @show convergence_status(epm)

        S_Q_EP = entropy(epm)
        @show S_Q_EP
        
        F_EP = free_energy(epm)
        @show F_EP

        isinf(F_EP) && break

        push!(Fs, F_EP)
        push!(Ss, S_Q_EP)
    end
    βs = βs[eachindex(Ss)]

    ϵ = 1e-3
    T = (x) -> (x .- minimum(x) .+ ϵ) / (maximum(x) .- minimum(x))
    scatter(
        Fs, Ss; 
        xlabel = "\$F_{EP}\$",
        ylabel = "\$S_{Q_{EP}}\$",
    )


end





## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
# ## ------------------------------------------------------------------
# ## ------------------------------------------------------------------
# # Normal Commulative
# Ncdf(x) = 0.5*(1.0+erf(big(x)/sqrt(2.0)))
# Ncdf(x, av, sd) = Ncdf((x - av) / sd)

# # (2π)^{n/2}|Σ|^{1/2}
# ZQ(Σ) = big(2*π)^(size(Σ, 1)/2)*sqrt(det(big.(Σ)))
# ZQ(epm::FluxEPModelT0) = ZQ(MetXEP._Σi(epm))

# # s -> var
# ZQn(Σ, μ, s, l, u) = ZQ(Σ)*(Ncdf(u, μ, sqrt(s)) - Ncdf(l, μ, sqrt(s)))
# function ZQn(epm::FluxEPModelT0, rxn)
#     rxn = rxnindex(epm, rxn)
#     μ = MetXEP._μ(epm, rxn)
#     s = MetXEP._s(epm, rxn)
#     l = lb(epm, rxn)
#     u = ub(epm, rxn)
#     Σ = MetXEP._Σi(epm)
#     return ZQn(Σ, μ, s, l, u)
# end

# function _F_EP(epm)
#     # F_EP = (N-1)logZ_Q - ∑logZ_Qn
#     N = length(epm.idxi)
#     _F = (N - 1)*log(ZQ(epm))
#     for i in epm.idxi
#         _F += -log(ZQn(epm, i))
#     end
#     return _F
# end



# ## ------------------------------------------------------------------
# let 
#     global net0 = pull_net("ecoli_core")
#     global net = box(net0, Gurobi.Optimizer)
#     global enet = EchelonMetNet(net)
#     global nbins = 200
#     global fbox = BoxGrid(enet, nbins)

#     global epm = FluxEPModelT0(net)
#     biom_id = extras(net, "BIOM")
#     beta!(epm, biom_id, 1e5)
#     converge!(epm)
#     @show convergence_status(epm)

#     global Q = MultivariateNormal(epm)

#     # @show MetXEP._normal_entropy(MetXEP._Σi(epm))
#     # @show entropy(Q)
#     # @show entropy(epm)
#     S_Q_EP = entropy(epm)
#     @show S_Q_EP
#     # S_grid = entropy(enet, fbox)
#     # @show S_grid
#     _F_EP = Float64(-F_EP(epm))
#     @show _F_EP

#     return nothing
# end

# ## ------------------------------------------------------------------
# # Grid entropy
# let 
#     # entropy(f::Function, enet::EchelonMetNet, fbox::BoxGrid)
    
# end