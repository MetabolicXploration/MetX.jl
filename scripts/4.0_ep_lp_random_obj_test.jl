## ------------------------------------------------------------------
@time begin
    using MetX
    using ProjFlows
    using Plots
    using Distributions
    import Gurobi
    import Clp
end

## ------------------------------------------------------------------
function _multi_linobj_optimize!(opm::OpModel, obj_idxs::Vector, obj_senses::Vector;
        eps = 0.0
    )
    _iter = zip(obj_idxs, obj_senses)
    _n = length(_iter)
    for (i, (obj_idx, obj_sense)) in enumerate(zip(obj_idxs, obj_senses))
        # optimize
        set_linear_obj!(opm, obj_idx, obj_sense)
        optimize!(opm)
        # let the last unbounded
        i == _n && break
        # bound
        sol = solution(opm, obj_idx)
        lb!(opm, obj_idx, sol - eps)
        ub!(opm, obj_idx, sol + eps)
    end
    return opm
end

## ------------------------------------------------------------------
# TODO: move to EP tests/tutorials
let

    # lep
    global net0 = pull_net("ecoli_core") 
    bounds!(net0, "SUCDi", 0.0, 0.0) # fix degeneracy
    global lep0 = lepmodel(net0)
    global lep1 = fva_strip(lep0, Clp.Optimizer)


    # objective setup
    for ops in [
            (;
                obj_idxs = ["BIOMASS_Ecoli_core_w_GAM"],
                obj_senses = [1.0],
                beta_factor = [5.0],
                beta_max_order = 7,
            ),
            (;
                obj_idxs = ["EX_o2_e", "ATPM"],
                obj_senses = [-1.0, 1.0],
                beta_factor = [1.0, 6.0],
                beta_max_order = 5,
            ),
            (;
                obj_idxs = ["EX_glc__D_e", "EX_o2_e", "ATPM"],
                obj_senses = [-1.0, -1.0, 1.0],
                beta_factor =  [3.0, 1.0, 6.0],
                beta_max_order = 6,
            ),
            (;
                obj_idxs = ["EX_glc__D_e", "EX_lac__D_e"],
                obj_senses = [1.0, 1.0],
                beta_factor = [6.0, 1.0],
                beta_max_order = 6,
            )
        ]

        # Extract
        obj_idxs = ops.obj_idxs
        obj_senses = ops.obj_senses
        beta_factor = ops.beta_factor
        beta_max_order = ops.beta_max_order

        println("\n", "="^40)

        # LP
        @show obj_idxs
        global opm = FBAOpModel(lep1, Clp.Optimizer)
        _multi_linobj_optimize!(opm, obj_idxs, obj_senses)
        lp_sol = solution(opm)
        @show solution(opm, obj_idxs)

        # # Check degeneracy
        # lep2 = deepcopy(lep1)
        # lb1, ub1 = fva(lep1, Gurobi.Optimizer)
        # lb!(lep2, obj_idxs, solution(opm, obj_idxs))
        # ub!(lep2, obj_idxs, solution(opm, obj_idxs))
        # lb2, ub2 = fva(lep2, Gurobi.Optimizer)
        # @show Float64(log10(prod(big, ub1 - lb1)))
        # @show Float64(log10(prod(big, ub2 - lb2)))

        # for (rxn, lb, ub) in zip(colids(lep2), lb2, ub2)
        #     if abs(ub - lb) > 1.0
        #         println(rxn)
        #     end
        # end

        # EP
        global epm = FluxEPModelT0(lep1)
        damp = 0.95
        verbose = false
        epsconv = 1e-6
        maxiter = 1000
        config!(epm; damp, verbose, epsconv, maxiter)
        converge!(epm) # beta zero

        
        beta_pool = 10.0.^range(-3, beta_max_order; length = 60)
        mse_vec = Float64[]
        max_err_vec = Float64[]
        ep_sol_dict = Dict()
        for _beta in beta_pool
            beta!(epm, obj_idxs, beta_factor .* obj_senses .* _beta)
            converge!(epm)
            ep_sol = mean(epm)
            mse = mean(x -> x^2, ep_sol - lp_sol)
            max_err = maximum(abs, (ep_sol - lp_sol) ./ max.(lp_sol, 1.0))

            for rxn in colids(epm)
                ep_sols = get!(ep_sol_dict, rxn, Float64[])
                push!(ep_sols, mean(epm, rxn))
            end

            push!(mse_vec, mse)
            push!(max_err_vec, max_err)
        end

        # test
        @assert minimum(mse_vec) < 1e-4
        @assert minimum(max_err_vec) < 1e-2
        @show minimum(mse_vec)
        @show minimum(max_err_vec)
        
        # Plots
        ps = Plots.Plot[]
        
        p = plot(;xlabel = "log10(beta)", ylabel = "log10(mse)")
        plot!(p, log10.(beta_pool), log10.(mse_vec); 
            label = "", c = :black, alpha = 0.7, lw = 2
        )
        push!(ps, p)
        
        p = plot(;xlabel = "log10(beta)", ylabel = "max sq error")
        plot!(p, log10.(beta_pool), max_err_vec; 
            label = "", c = :black, alpha = 0.7, lw = 2
        )
        # push!(ps, p)

        p = plot(;xlabel = "log10(beta)", ylabel = "ep - lp")
        for rxn in colids(epm)
            ep = get!(ep_sol_dict, rxn, Float64[])
            lp = solution(opm, rxn)
            plot!(p, log10.(beta_pool), ep .- lp; 
                label = "", c = :black, alpha = 0.7, lw = 2
            )
        end
        push!(ps, p)
        
        # Save
        objid = join(obj_idxs, "-")
        sfig(ps, 
            [@__DIR__], 
            "ep_to_lp", objid,
            ".png";
            layout = (2, 1)
        ) |> println

    end # for ops 
end