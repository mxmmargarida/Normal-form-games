struct NormalGame
    n::Int # number of players
    strat::Vector{Int} # number of strategies for each player
    polymatrix::Dict{Tuple{Int, Int}, Matrix{Float64}} # matrices for nomal form representation
end


# generate randin matrix
function GenerateGame(n::Int = 2,max_strat::Int =10)
    # n - number of players; by default, it is 2
    # max_strat - maximum number of strategies by player; by default, it is 10
    if n <= 0 || max_strat <= 0 
        throw(ArgumentError("Inputs must be strictly positive integers"))
    end
    # for each player decide the number of strategies
    strat = rand(1:max_strat,n)
    # for each pair of players create two matrices
    polymatrix = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    for p1 in 1:n
        for p2 in p1:n
            if p1 == p2 # it is assumed that no utility is produced independently from the other players (as this could be included in one of the polymatrices). Still, if we add that in these tables, the code for computing equilibria remains correct.
                polymatrix[(p1, p2)] = zeros(strat[p1], strat[p2])
            else
                # Assign matrices for (p1, p2) pair
                polymatrix[(p1, p2)] = rand(1:100, strat[p1], strat[p2])
                polymatrix[(p2, p1)] = rand(1:100, strat[p2], strat[p1])
            end
        end
    end
    return NormalGame(n,strat,polymatrix)
end

# save instance
function save_game_to_json(filename::AbstractString, my_game::NormalGame)
    json_string = JSON.json(my_game)
    # Write the JSON string to a file
    open(filename, "w") do io
    write(io, json_string)
    end
end

# Read instance
function read_game_from_json(filename::AbstractString)::NormalGame
    # Read the JSON data from the file
    json_data = JSON.parsefile(filename)

    # Extract fields from the JSON data
    n = json_data["n"]
    strat = json_data["strat"]
 
    # Initialize an empty dictionary for the polymatrix
    polymatrix = Dict{Tuple{Int, Int}, Matrix{Float64}}()

    # Populate the polymatrix dictionary
    for (key, value) in json_data["polymatrix"]
        # Convert the string key to a tuple of integers
        tuple_key = eval(Meta.parse(key))
        
        # Convert the array of arrays to a matrix of Float64
        if length(value)>1
            matrix_value = reduce(hcat,value)
        else
            matrix_value = reshape(Float64.(value[1]), (length(value[1]),1))
        end
        
        # Assign the matrix to the corresponding key in the dictionary
        polymatrix[tuple_key] = matrix_value
    end
    
    # Create and return the NormalGame struct
    return NormalGame(n, strat, polymatrix)
end