module Sigil

# using DifferentialEquations
# using Plots

mutable struct Star
    M::AbstractFloat
    m::AbstractFloat
    mf::AbstractFloat
    X::AbstractFloat
    Y::AbstractFloat
    grid::Array{AbstractFloat}
    r::Array{AbstractFloat}
    l::Array{AbstractFloat}
    P::Array{AbstractFloat}
    T::Array{AbstractFloat}
    filename::String
end

Star(M) = Star(M, 0.01 * Msun, 0.8 * M, 0.7, 0.28,
               zeros(1), zeros(1), zeros(1), zeros(1), zeros(1), "")
Star(M, X, Y) = Star(M, 0.01 * Msun, 0.8 * M, X, Y,
                     zeros(1), zeros(1), zeros(1), zeros(1), zeros(1), "") 
Star(M, filename) = Star(M, 0.01 * Msun, 0.8 * M, 0.7, 0.28,
                     zeros(1), zeros(1), zeros(1), zeros(1), zeros(1), filename) 

export score,
    init_guess,
    deriv!,
    Msun,
    Rsun,
    Star,
    load1,
    load2,
    shootf!,
    plot_star
    

include("constants.jl")
include("opacity.jl")
include("equations.jl")
include("physics.jl")
include("shootf.jl")
include("plot.jl")

end
