# Basic usage
# first, do 
#] activate .
# in NormalGames folder (outside src)
# load module
using NormalGames
println("--- generating game with 3 players and maximum 5 strategies ---")
# any function in NormalGames needs this prefix to be used
# we also need to load packages with using f used outside the module
MyGame = NormalGames.GenerateGame(3,7)

println("--- Access the info in MyGame ---")
println("number of players: ", MyGame.n)
println("number of strategies: ", MyGame.strat)
println("Normal form game between player 1 and 2 ", MyGame.polymatrix[(1,2)])
println("Normal form game between player 2 and 1 ",  MyGame.polymatrix[(2,1)])

println("Compute Nash Equilibria")
# without profiling
println("Indicator constraints method")
NE_mixed1, NE_u1, total_time1, _, numb_NE1 = NormalGames.All_NE(MyGame,1)
println("Formulation 1 by Sandholm et al. (2005)")
NE_mixed2, NE_u2,  total_time2, _, numb_NE3 = NormalGames.All_NE(MyGame,2)
# println("Formulation 2 by Sandholm et al. (2005)")
# Regret_val, NE_u3, NE_mixed3 = NashEquilibria3(MyGame)
# println("Formulation 3 by Sandholm et al. (2005)")
# Regret_prob_val, NE_u4, NE_mixed4 = NashEquilibria4(MyGame)
# println("Formulation 4 by Sandholm et al. (2005)")
# Regret_comb_val, NE_u5, NE_mixed5 = NashEquilibria5(MyGame)

println("Indicator constraints method with callback")
NE_mixed_callback, NE_u_callback, total_time_callback, time_call_back, numb_call_back = NormalGames.All_NE(MyGame, "callback 1")
println("Formulation 1 by Sandholm et al. (2005) with callback")
NE_mixed2_callback, NE_u2_callback,  total_time2_callback, numb_NE3_callback = NormalGames.All_NE(MyGame,"callback 2")

# Save results
# Sample data
using DataFrames, CSV
data = DataFrame(Numb_NE = [length(NE_mixed1),length(NE_mixed2),length(NE_mixed_callback),length(NE_mixed2_callback)],Name = ["Indicator Method", "Formulation 1", "Indicator Method Callback", "Formulation 2 Callback"],Time = [total_time1, total_time2, total_time_callback,total_time2_callback])
CSV.write("output.csv", data)

# Save instance
NormalGames.save_game_to_json("Problems/MyGame.json", MyGame)
loaded_game = NormalGames.read_game_from_json("Problems/MyGame.json")

loaded_game.polymatrix == MyGame.polymatrix

# profiling: runs code many times
# using BenchmarkTools
# println("Indicator constraints method")
# b1 = @benchmark NE_mixed, NE_u = All_NE(MyGame,1)
# println("Formulation 1 by Sandholm et al. (2005)")
# b2 = @benchmark NE_mixed, NE_u = All_NE(MyGame,2)
# println("Indicator constraints benchmark ", b1)
# println("Formulation 1 by Sandholm et al. (2005) ", b2)

# OPT, utilities, mixed, supp, m, y_var = NashEquilibria(MyGame)
# println("Social welfare = ", OPT)
# println("Mixed strategies of NE ", mixed)
# @constraint(m, sum(sum(1-y_var[p][j] for j in 1:MyGame.strat[p] if supp[p][j]>=0.5)+sum(y_var[p][j] for j in 1:MyGame.strat[p] if supp[p][j]<0.5) for p in 1:MyGame.n)>=1) 

# With correlated equilibria, we can only compute 1. If we want to compute all, then we need to find a constraint that allow us to have y=1 if and only if tau>0. It is not obvious how to do it.
println("Correlated equilibria with the support of the first NE found")
# restrict CE computation to the support of the first NE found
using IterTools
indices = collect(product(ntuple(i -> 1:MyGame.strat[i], length(MyGame.strat))...)) 
supp = []
for ind in indices
    aux = 1
    for p in 1:MyGame.n
        aux *= NE_mixed2_callback[1][p][ind[p]]
    end
    if aux > 10^(-6)
        push!(supp,ind)
    end
end
# note that supp is optional
SW, CE_mixed, CE_u, CE_time = NormalGames.CorrelatedEquilibria(MyGame,supp)
SW1, CE_mixed1, CE_u1, CE_time1 = NormalGames.CorrelatedEquilibria(MyGame)

println("CE social welfare is ", SW, " and for the NE is ", [sum(ne) for ne in NE_u2_callback])
println("The CE with maximum social welfare is ", SW1)

# Ongoing step: all functions for NE verify if the polymatrices are negative. If they are an error is given. So missing to do a preprocessing to make all entries non-negative (add it to Instances). Double check that this is not needed for CE.

# The next step is to code SGM (or the RL approach from the paper shared with Goolnosh) and think how to use callbacks for computing all correlated equilibria


