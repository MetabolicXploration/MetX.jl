## ------------------------------------------------------------------
import Pkg

## ------------------------------------------------------------------
include("metx_base.jl")

## ------------------------------------------------------------------
let 
    curr = Base.active_project()
    try
        for pkg in MET_X_PKG_NAMES

            println("-"^60)
            path = joinpath(Pkg.devdir(), pkg)

            # activate
            Pkg.activate(path)
            
            # dev MetX packages
            _foreach_dep() do dep
                dep.name in MET_X_PKG_NAMES || return
                println("dev dep: ", dep.name)
                Pkg.develop(dep.name)
            end

            # resolve
            Pkg.resolve()

            println()
        end
        finally; Pkg.activate(curr)
    end
end