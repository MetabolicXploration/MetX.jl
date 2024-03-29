using MetX
using Test
import Pkg

@testset "MetX.jl" begin
    
    # TODO: Do this using the Test interface.
    totest =  [
        "MetXBase", "MetXNetHub", "MetXOptim", "MetXEP", 
        "MetXMC", "MetXGrids"
    ]

    for pkg in totest
        
        println("\n"^3, "="^60, "\n")
        println(uppercase(pkg))
        println("\n", "="^60, "\n"^3)

        Pkg.test(pkg)
    end

end
