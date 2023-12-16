# ------------------------------------------------------------------
@time begin
    using MetXBase
    using MetXOptim
    using Gurobi
    using MetXEP
    using MetXMC
    using MetXNetHub
    using Plots
    using ProjAssistant
    using Distributions
end

## ------------------------------------------------------------------
# TODO: Add to tutorials

# /Users/Pereiro/.julia/dev/Metabolic-EP/data

## ------------------------------------------------------------------
function _trunc_marginal(epm::FluxEPModelT0, rxn)
    rxn = rxnindex(epm, rxn)
    
    μ = [epm.μd; epm.μi]
    μ = μ[epm.idxmap_inv][rxn] * epm.scalefact

    s = [epm.sd; epm.si]
    s = s[epm.idxmap_inv][rxn] * epm.scalefact^2

    l, u = lb(epm, rxn), ub(epm, rxn)

    return truncated(Normal(μ, sqrt(s)), l, u)
end

## ------------------------------------------------------------------
# Test marginals
let 
    global model_id = "toy_net4D"
    global net0 = pull_net(model_id)
    global net = box(net0, Gurobi.Optimizer)
    
    global epm = FluxEPModelT0(net)
    converge!(epm)
    @show convergence_status(epm)
    global bmn = EPBoxedMvNormal(epm)
    global hrm = HRModel(net, Gurobi.Optimizer)
    sample!(hrm, Int(1e5)) # warmup

    
    # Marginals
    ps = Plots.Plot[]
    
    for rxn in net.rxns
        
        p = plot(; title = rxn, xlabel = "norm. flux")

        rxni = rxnindex(net, rxn)

        nsamples = Int(5e6)
        nbins = 100
        l, u = bounds(net, rxn)
        bins = range(l, u; length = nbins)
        bins_norm = (bins .- minimum(bins)) ./ (maximum(bins) .- minimum(bins))
        
        # Analytical
        tN = _trunc_marginal(epm, rxn)
        hist = pdf.([tN], bins)
        hist = hist ./ sum(hist)
        plot!(p, bins_norm, hist; 
            lw = 3, c = :red,
            label = string("EP Analytical")
        )
        
        # EP Sampled
        @time bins, hist = sample_histogram!(bmn, rxni, bins; 
            nsamples, rw = one)
        hist = hist ./ sum(hist)
        plot!(p, bins_norm, hist; 
            lw = 3, c = :blue,
            label = string("EP Sampled"),
        )

        # HR Sampled
        nsamples = Int(5e6)
        @time bins, hist = sample_histogram!(hrm, rxni, bins; 
            nsamples, rw = one)
        hist = hist ./ sum(hist)
        plot!(p, bins_norm, hist; 
            lw = 3, c = :black,
            label = string("HR Sampled"),
        )

        push!(ps, p)
    end

    fn = sfig(ps, [@__DIR__], 
        basename(@__FILE__), model_id, "marginals_test", ".png";
    )
    @show fn
end
