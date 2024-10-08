# TODO: Use Preferences.jl for configuration
module MetX

    using Reexport

    #! include .
    include("devmode.jl")

    @reexport using MetXBase
    @reexport using MetXGEMs
    @reexport using MetXOptim
    @reexport using MetXEP
    @reexport using MetXMC
    # @reexport using MetXGrids
    @reexport using MetXNetHub
    @reexport using MetXCultureHub

end