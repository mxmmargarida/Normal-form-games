module NormalGames

using Random
using JSON

using JuMP, Gurobi
using IterTools
using Combinatorics

using DataFrames, CSV

include("NormalForm/Instances.jl")
include("NormalForm/Equilibria.jl")

end # module NormalGames
