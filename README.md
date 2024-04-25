# NormalGames <!-- omit from toc -->

Julia library for the equilibria computation of normal form (also known as strategic form) games.

## Table of Contents <!-- omit from toc -->

- [Installation](#installation)
- [Basic Usage](#basic-usage)
  - [Generation of instances](#generation-of-instances)
  - [Computation of equilibria](#computation-of-equilibria)
- [References](#references)

## Installation

The package is not available on the registry. Please clone the repository and
add it as a [local package](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-a-local-package). You also need to add [JuMP](https://jump.dev/JuMP.jl/stable/) and [Gurobi](https://github.com/jump-dev/Gurobi.jl).

```julia
]add https://github.com/mxmmargarida/Normal-form-games.git
```

## Basic Usage

In the following, we essentially go over test.jl

Start by activating the Julia enviroment:

```julia
] activate .
```

Load the module:

```julia
using NormalGames
```

If the above raises an error, it may mean that you need to run the commands below to install recorded dependencies
```julia
using Pkg
Pkg.instantiate()
```

### Generation of instances

Generate an instance with 3 players and at most 7 pure strategies per player:

```julia
MyGame = NormalGames.GenerateGame(3,7)
```

Access the information of the generated game:

```julia
println("number of players: ", MyGame.n)
println("number of strategies: ", MyGame.strat)
println("Normal form game between player 1 and 2 ", MyGame.polymatrix[(1,2)])
println("Normal form game between player 2 and 1 ",  MyGame.polymatrix[(2,1)])
```

For instance the matrix in ```MyGame.polymatrix[(2,1)]```

Save an instance:

```julia
NormalGames.save_game_to_json("Problems/MyGame.json", MyGame)
```

Load a saved instance:

```julia
loaded_game = NormalGames.read_game_from_json("Problems/MyGame.json")
```

### Computation of equilibria

Compute all Nash equilibria using indicator constraints and a cutting plane method (no use of lazy constraints):

```julia
println("Indicator constraints method")
NE_mixed1, NE_u1, total_time1, _, numb_NE1 = NormalGames.All_NE(MyGame,1)
```

Compute all Nash equilibria using indicator constraints and a branch-and-cut method (use of lazy constraints):

```julia
println("Indicator constraints method with callback")
NE_mixed_callback, NE_u_callback, total_time_callback, time_call_back, numb_call_back = NormalGames.All_NE(MyGame, "callback 1")
```

Compute all Nash equilibria using Formlulation 1 of \[[1](#readme-ref1)\] (Big-M formulation), and a cutting plane method (no use of lazy constraints):

```julia
println("Formulation 1 by Sandholm et al. (2005)")
NE_mixed2, NE_u2,  total_time2, _, numb_NE3 = NormalGames.All_NE(MyGame,2)
```

Compute all Nash equilibria using Formlulation 1 of \[[1](#readme-ref1)\] (Big-M formulation), and a branch-and-cut method (use of lazy constraints):

```julia
println("Formulation 1 by Sandholm et al. (2005) with callback")
NE_mixed2_callback, NE_u2_callback,  total_time2_callback, numb_NE3_callback = NormalGames.All_NE(MyGame,"callback 2")
```

Compute one Nash equilibrium using Formlulation 2 of \[[1](#readme-ref1)\]

```julia
println("Formulation 2 by Sandholm et al. (2005)")
Regret_val, NE_u3, NE_mixed3 = NormalGames.NashEquilibria3(MyGame)
```

Compute one Nash equilibrium using Formlulation 3 of \[[1](#readme-ref1)\]

```julia
println("Formulation 3 by Sandholm et al. (2005)")
Regret_prob_val, NE_u4, NE_mixed4 = NormalGames.NashEquilibria4(MyGame)
```

Compute one Nash equilibrium using Formlulation 3 of \[[1](#readme-ref1)\]

```julia
println("Formulation 4 by Sandholm et al. (2005)")
Regret_comb_val, NE_u5, NE_mixed5 = NormalGames.NashEquilibria5(MyGame)
```

Save results:

```julia
using DataFrames, CSV
data = DataFrame(Numb_NE = [length(NE_mixed1),length(NE_mixed2),length(NE_mixed_callback),length(NE_mixed2_callback)],Name = ["Indicator Method", "Formulation 1", "Indicator Method Callback", "Formulation 2 Callback"],Time = [total_time1, total_time2, total_time_callback,total_time2_callback])
CSV.write("output.csv", data)
```

Compute one correlated equilibria. It is not possible to compute more than one (for that need to find a constraint that allow us to have y=1 if and only if tau>0. It is not obvious how to do it.) 

```julia
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
```

## References

<a id="readme-ref1"></a> \[1\] Sandholm, Tuomas, Andrew Gilpin, Vincent Conitzer. Mixed-integer programming methods for finding Nash equilibria. In AAAI, p. 495-501, 2005. ([link](https://cdn.aaai.org/AAAI/2005/AAAI05-078.pdf))