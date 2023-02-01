## ------------------------------------------------------------------
import Pkg
let 
    curr = Base.active_project()
    try
        for pkg in [
                "MetX", 
                "MetXBase", 
                "MetXOptim", 
                "MetXGrids", "MetXMC", "MetXEP", 
                "MetXPlots", 
                "MetXNetHub", "MetXCultureHub", 
                "MetXTutorials"
            ]
            println("-"^60)
            path = joinpath(Pkg.devdir(), pkg)
            Pkg.activate(path)
            Pkg.resolve()
            println()
        end
        finally; Pkg.activate(curr)
    end
end