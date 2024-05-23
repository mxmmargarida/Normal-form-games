# Compute Nash equilibria: MIP approach

function has_negative_entry(PolyMatrix::Dict{Tuple{Int, Int}, Matrix{Float64}})
    for (key, matrix) in PolyMatrix
        if any(x -> x < 0, matrix)
            return true
        end
    end
    return false
end

"""
    NashEquilibria(my_game,model="NaN",supp=[],y=[],v=[],x=[])

    Method using indicator constraints
"""
function NashEquilibria(my_game,model="NaN",supp=[],y=[],v=[],x=[])
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    if model != "NaN"
        @constraint(model, sum(sum(1-y[p][j] for j in 1:my_game.strat[p] if supp[p][j]>=0.5)+sum(y[p][j] for j in 1:my_game.strat[p] if supp[p][j]<0.5) for p in 1:my_game.n)>=1) 
    else
        # start model
        model = Model(Gurobi.Optimizer)

        set_silent(model)
        set_time_limit_sec(model, 300.0)
        MOI.set(model, MOI.NumberOfThreads(), 1)

        # Variables
        # define probability distribution for each player
        x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
        # define indicator variable for each strategy of each player
        y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
        # define expected utility value for each player
        @variable(model, v[i = 1:my_game.n])

        # Constraints
        # if x>0 then y=1
        @constraint(model, y .>= x) # componentwise constraint
        for p in 1:my_game.n
            # add probability distribution constraint
            @constraint(model, sum(x[p]) == 1)
            # if x>0 (then y = 1) then sum(polymatrix[(p,i)] i=1:my_game.n) = max utility of p
            for j in 1:my_game.strat[p]
                profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
                @constraint(model, profitj <= v[p])
                @constraint(model, y[p][j] --> {profitj == v[p]})
            end
        end

        # Objective function
        @objective(model, Max, sum(v))
    end

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        # println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), value.(v), [value.(x[p]) for p in 1:my_game.n], [value.(y[p]) for p in 1:my_game.n], model, y, v, x
    else
        # println("Model is not solved and/or not feasible.")
        return "NaN","NaN","NaN", "NaN", model, y, v, x
    end
end

function NashEquilibria_CallBack(my_game)
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    # start model
    model = Model(Gurobi.Optimizer)
    model[:callback_time] = 0.
    model[:callback_calls] = 0.
    model[:NE_mixed] = []
    model[:NE_u] = []

    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution for each player
    x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
    # define indicator variable for each strategy of each player
    y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
    # define expected utility value for each player
    @variable(model, v[i = 1:my_game.n])

    # Constraints
    # if x>0 then y=1
    @constraint(model, y .>= x) # componentwise constraint
    for p in 1:my_game.n
        # add probability distribution constraint
        @constraint(model, sum(x[p]) == 1)
        # if x>0 (then y = 1) then sum(polymatrix[(p,i)] i=1:my_game.n) = max utility of p
        for j in 1:my_game.strat[p]
            profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
            @constraint(model, profitj <= v[p])
            @constraint(model, y[p][j] --> {profitj == v[p]})
        end
    end

    # Objective function
    @objective(model, Max, sum(v))

    # Tell JuMP to use lazy constraints with the callback function
    MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> no_good_cut_callback(cb_data,x,y,v,model,my_game))

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        # println("Model is solved and feasible and solving time is ",solve_time(model))
        return model[:NE_mixed], model[:NE_u], solve_time(model), model[:callback_time], model[:callback_calls]
    else
        # println("Model is not solved and/or not feasible, and solve time is ",solve_time(model))
        return model[:NE_mixed], model[:NE_u], solve_time(model), model[:callback_time], model[:callback_calls]
    end
end

function no_good_cut_callback(cb_data,x,y,v,model,my_game)
    status = callback_node_status(cb_data, model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        push!(model[:NE_mixed],[callback_value.(cb_data,x[p]) for p in 1:my_game.n])
        push!(model[:NE_u],callback_value.(cb_data,v))
        model[:callback_time] += @elapsed begin # accumulated time to use call back
            supp = [callback_value.(cb_data,y[p]) for p in 1:my_game.n]
            con = @build_constraint(sum(sum(1-y[p][j] for j in 1:my_game.strat[p] if supp[p][j]>=0.5)+sum(y[p][j] for j in 1:my_game.strat[p] if supp[p][j]<0.5) for p in 1:my_game.n)>=1) 
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
        end
        model[:callback_calls] += 1 # count number of callbacks
    end
end

"""
    NashEquilibria2(my_game,model="NaN",supp=[],y=[],v=[],x=[])

    Method using by Sandholm et al. (2005) extended to polymatrix games: Formulation 1
"""
function NashEquilibria2(my_game,model="NaN",supp=[],y=[],v=[],x=[])
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    if model != "NaN"
        @constraint(model, sum(sum(1-y[p][j] for j in 1:my_game.strat[p] if supp[p][j]>=0.5)+sum(y[p][j] for j in 1:my_game.strat[p] if supp[p][j]<0.5) for p in 1:my_game.n)>=1) 
    else
        # start model
        model = Model(Gurobi.Optimizer)

        set_silent(model)
        set_time_limit_sec(model, 300.0)
        MOI.set(model, MOI.NumberOfThreads(), 1)

        # Variables
        # define probability distribution for each player
        x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
        # define indicator variable for each strategy of each player
        y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
        # define expected utility value for each player
        @variable(model, v[i = 1:my_game.n])
        # define expected objective value of strategy sj for player i (same notation as the paper)
        u = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
        # define regret of strategy sj for player i (same notation as the paper)
        r = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]

        # Big-M: maximum difference between two strategies
        V = []
        for p in 1:my_game.n
            aux_v = 0
            for j in 1:my_game.strat[p]
                for i in 1:my_game.strat[p]
                    aux_j_i = sum(maximum(my_game.polymatrix[(p,k)][j,:] - my_game.polymatrix[(p,k)][i,:]) for k in 1:my_game.n)
                    if aux_j_i > aux_v
                        aux_v = aux_j_i
                    end
                end
            end
            push!(V, aux_v)
        end

        # Constraints
        for p in 1:my_game.n
            # add probability distribution constraint
            @constraint(model, sum(x[p]) == 1)
            for j in 1:my_game.strat[p]
                profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
                @constraint(model, u[p][j] == profitj)
                @constraint(model, v[p] >= u[p][j])
                @constraint(model, r[p][j] == v[p]-u[p][j])
                @constraint(model, x[p][j] <= 1- y[p][j])
                @constraint(model, r[p][j] <= V[p]*y[p][j])
            end
        end

        # Objective function
        @objective(model, Max, sum(v))
    end

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        # println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), value.(v), [value.(x[p]) for p in 1:my_game.n], [value.(y[p]) for p in 1:my_game.n], model, y, v, x
    else
        # println("Model is not solved and/or not feasible.")
        return "NaN","NaN","NaN", "NaN", model, y, v, x
    end
end

function NashEquilibria_CallBack2(my_game,model="NaN",supp=[],y=[],v=[],x=[])
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    # start model
    model = Model(Gurobi.Optimizer)
    model[:callback_time] = 0.
    model[:callback_calls] = 0.
    model[:NE_mixed] = []
    model[:NE_u] = []

    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution for each player
    x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
    # define indicator variable for each strategy of each player
    y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
    # define expected utility value for each player
    @variable(model, v[i = 1:my_game.n])
    # define expected objective value of strategy sj for player i (same notation as the paper)
    u = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # define regret of strategy sj for player i (same notation as the paper)
    r = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]

    # Big-M: maximum difference between two strategies
    V = []
    for p in 1:my_game.n
        aux_v = 0
        for j in 1:my_game.strat[p]
            for i in 1:my_game.strat[p]
                aux_j_i = sum(maximum(my_game.polymatrix[(p,k)][j,:] - my_game.polymatrix[(p,k)][i,:]) for k in 1:my_game.n)
                if aux_j_i > aux_v
                    aux_v = aux_j_i
                end
            end
        end
        push!(V, aux_v)
    end

    # Constraints
    for p in 1:my_game.n
        # add probability distribution constraint
        @constraint(model, sum(x[p]) == 1)
        for j in 1:my_game.strat[p]
            profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
            @constraint(model, u[p][j] == profitj)
            @constraint(model, v[p] >= u[p][j])
            @constraint(model, r[p][j] == v[p]-u[p][j])
            @constraint(model, x[p][j] <= 1- y[p][j])
            @constraint(model, r[p][j] <= V[p]*y[p][j])
        end
    end

    # Objective function
    @objective(model, Max, sum(v))

    # Tell JuMP to use lazy constraints with the callback function
    MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> no_good_cut_callback(cb_data,x,y,v,model,my_game))

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        # println("Model is solved and feasible and solving time is ",solve_time(model))
        return model[:NE_mixed], model[:NE_u], solve_time(model), model[:callback_time], model[:callback_calls]
    else
        # println("Model is not solved and/or not feasible, and solve time is ",solve_time(model))
        return model[:NE_mixed], model[:NE_u], solve_time(model), model[:callback_time], model[:callback_calls]
    end
end

function All_NE(my_game, method = 1)
    NE_mixed = []
    NE_u = []
    count_NE = 0
    total_time = 0
    if method == 1
        total_time += @elapsed begin
            OPT, utilities, mixed, supp, m, y_var, v_var, x_var = NashEquilibria(my_game)
        end
    elseif method == "callback 1"
        return NashEquilibria_CallBack(my_game)
    elseif method == "callback 2"
        return NashEquilibria_CallBack2(my_game)
    else
        total_time += @elapsed begin
            OPT, utilities, mixed, supp, m, y_var, v_var, x_var = NashEquilibria2(my_game)
        end
    end
    total_time += @elapsed begin
        while OPT != "NaN"
            count_NE = count_NE+1
            push!(NE_mixed,mixed)
            push!(NE_u,utilities)
            if method == 1
                OPT, utilities, mixed, supp, m, y_var = NashEquilibria(my_game,m,supp,y_var, v_var, x_var)
            # elseif method == 2
            else
                OPT, utilities, mixed, supp, m, y_var = NashEquilibria2(my_game,m,supp,y_var, v_var, x_var)
            end
        end
    end
    return NE_mixed, NE_u, total_time, "NaN", count_NE
end

"""
    NashEquilibriaPNS(my_game)

    Method using the enumeration procedure by Porter et al. (2008)
    This has options to compute all equilibria or stop when one is found, to compute only pure equilibria, and to run until a certain time limit is achieved.
"""
function NashEquilibriaPNS(my_game,all_NE = true, just_pure =false, time_limit = false, time_stop = 60)
    if just_pure
        supp_size = collect(product(ntuple(i -> 1:1, my_game.n)...))
    else
        supp_size = collect(product(ntuple(i -> 1:my_game.strat[i], my_game.n)...))
    end
    ### sort support sizes: flatten supp_size before sorting
    sorted_indices = sort(supp_size[:], by=t -> (sum(t), max_abs_diff(t)))
    ### solve feasibility problem for each assignment of strategies for the support size
    # list of expected utilities for each Nash equilibrium
    NE_u = []
    # list of Nash equilibria
    NE_mixed = []
    stop_enumeration = false
    start_time = time()
    total_time = @elapsed begin
        for s in sorted_indices
            # all supports of size s[i]
            combinations_per_dimension = [Combinatorics.combinations(1:my_game.strat[i], s[i]) for i in 1:my_game.n]
            combined_combinations = collect(Iterators.product(combinations_per_dimension...))
            for supp in combined_combinations
                OPT, U_opt, xOpt = FeasibilityProgram(my_game,supp)
                # if there is an equilibrium it must be added
                if OPT !="NaN"
                    push!(NE_u,U_opt)
                    push!(NE_mixed,xOpt)
                    if all_NE==false
                        stop_enumeration = true
                        break
                    end
                end
                # if we use a time limit than, when we go over it, we must stop
                if time_limit && time()-start_time > time_stop
                    return "TL", NE_u, NE_mixed
                end
            end
            if stop_enumeration
                break
            end
        end
    end
    return total_time, NE_u, NE_mixed
end

# function used on the sorting of the support sizes
max_abs_diff = t -> maximum(abs(t[i] - t[j]) for i in 1:length(t) for j in i+1:length(t))

# auxiliar function for the PNS method
function FeasibilityProgram(my_game,supp)
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution for each player
    x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
    # define expected utility value for each player
    @variable(model, v[i = 1:my_game.n])

    # Constraints
    for p in 1:my_game.n
        # add probability distribution constraint
        @constraint(model, sum(x[p]) == 1)
        for j in 1:my_game.strat[p]
            if j in supp[p]
                @constraint(model, v[p] == sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n))
            else
                @constraint(model, v[p] >= sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n))
                @constraint(model, x[p][j] == 0)
            end
        end
    end
    
    # Objective function
    @objective(model, Max, sum(v))

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        # println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), value.(v), [value.(x[p]) for p in 1:my_game.n]
    else
        # println("Model is not solved and/or not feasible.")
        return "NaN","NaN","NaN"
    end
end

"""
    NashEquilibria3(my_game)

    Method using by Sandholm et al. (2005) extended to polymatrix games: Formulation 2
    This method cannot be used to compute multiple equilibria as its feasible region contains non-equilibrium solutions
"""
function NashEquilibria3(my_game)
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    # start model
    model = Model(Gurobi.Optimizer)

    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution for each player
    x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
    # define indicator variable for each strategy of each player
    y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
    # define expected utility value for each player
    @variable(model, v[i = 1:my_game.n])
    # define expected objective value of strategy sj for player i (same notation as the paper)
    u = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # define regret of strategy sj for player i (same notation as the paper)
    r = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # f_si for player i (same notation as the paper)
    f = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]

    # Big-M: maximum difference between two strategies
    V = []
    for p in 1:my_game.n
        aux_v = 0
        for j in 1:my_game.strat[p]
            for i in 1:my_game.strat[p]
                aux_j_i = sum(maximum(my_game.polymatrix[(p,k)][j,:] - my_game.polymatrix[(p,k)][i,:]) for k in 1:my_game.n)
                if aux_j_i > aux_v
                    aux_v = aux_j_i
                end
            end
        end
        push!(V, aux_v)
    end

    # Constraints
    for p in 1:my_game.n
        # add probability distribution constraint
        @constraint(model, sum(x[p]) == 1)
        for j in 1:my_game.strat[p]
            profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
            @constraint(model, u[p][j] == profitj)
            @constraint(model, v[p] >= u[p][j])
            @constraint(model, r[p][j] == v[p]-u[p][j])
            @constraint(model, x[p][j] <= 1- y[p][j])
            @constraint(model, f[p][j] >= r[p][j])
            @constraint(model, f[p][j] >= V[p]*y[p][j])
        end

        # Objective function: minimize the sum of the regrets of strategies that have positive probability
        @objective(model, Min, sum(sum(f[p][j]-V[p]*y[p][j] for j in 1:my_game.strat[p]) for p in 1:my_game.n))
    end

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), value.(v), [value.(x[p]) for p in 1:my_game.n] #, [value.(y[p]) for p in 1:my_game.n], model, y, v, x
    else
        println("Model is not solved and/or not feasible.")
        return "NaN","NaN","NaN" #, "NaN", model, y, v, x
    end
end

"""
    NashEquilibria4(my_game)

    Method using by Sandholm et al. (2005) extended to polymatrix games: Formulation 3
    This method cannot be used to compute multiple equilibria as its feasible region contains non-equilibrium solutions
"""
function NashEquilibria4(my_game)
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    # start model
    model = Model(Gurobi.Optimizer)

    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution for each player
    x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
    # define indicator variable for each strategy of each player
    y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
    # define expected utility value for each player
    @variable(model, v[i = 1:my_game.n])
    # define expected objective value of strategy sj for player i (same notation as the paper)
    u = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # define regret of strategy sj for player i (same notation as the paper)
    r = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # g_si for player i (same notation as the paper)
    g = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]

    # Big-M: maximum difference between two strategies
    V = []
    for p in 1:my_game.n
        aux_v = 0
        for j in 1:my_game.strat[p]
            for i in 1:my_game.strat[p]
                aux_j_i = sum(maximum(my_game.polymatrix[(p,k)][j,:] - my_game.polymatrix[(p,k)][i,:]) for k in 1:my_game.n)
                if aux_j_i > aux_v
                    aux_v = aux_j_i
                end
            end
        end
        push!(V, aux_v)
    end

    # Constraints
    for p in 1:my_game.n
        # add probability distribution constraint
        @constraint(model, sum(x[p]) == 1)
        for j in 1:my_game.strat[p]
            profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
            @constraint(model, u[p][j] == profitj)
            @constraint(model, v[p] >= u[p][j])
            @constraint(model, r[p][j] == v[p]-u[p][j])
            @constraint(model, g[p][j] >= x[p][j])
            @constraint(model, g[p][j] >= 1-y[p][j])
            @constraint(model, r[p][j] <= V[p]*y[p][j])
        end

        # Objective function: minimize the probabilities of the strategies with positive regret
        @objective(model, Min, sum(sum(g[p][j]-(1-y[p][j]) for j in 1:my_game.strat[p]) for p in 1:my_game.n))
    end

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), value.(v), [value.(x[p]) for p in 1:my_game.n] #, [value.(y[p]) for p in 1:my_game.n], model, y, v, x
    else
        println("Model is not solved and/or not feasible.")
        return "NaN","NaN","NaN" #, "NaN", model, y, v, x
    end
end

"""
    NashEquilibria5(my_game)

    Method using by Sandholm et al. (2005) extended to polymatrix games: Formulation 4
    This method cannot be used to compute multiple equilibria as its feasible region contains non-equilibrium solutions
"""
function NashEquilibria5(my_game)
    if has_negative_entry(my_game.polymatrix)
        throw(ArgumentError("A polymatrix game must have no negative entries. A pre-processing is necessary"))
    end
    # start model
    model = Model(Gurobi.Optimizer)

    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution for each player
    x = [@variable(model, [1:my_game.strat[i]], lower_bound = 0.0, upper_bound = 1.0) for i = 1:my_game.n]
    # define indicator variable for each strategy of each player
    y = [@variable(model, [1:my_game.strat[i]], Bin) for i = 1:my_game.n]
    # define expected utility value for each player
    @variable(model, v[i = 1:my_game.n])
    # define expected objective value of strategy sj for player i (same notation as the paper)
    u = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # define regret of strategy sj for player i (same notation as the paper)
    r = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # f_si for player i (same notation as the paper)
    f = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]
    # g_si for player i (same notation as the paper)
    g = [@variable(model, [1:my_game.strat[i]]) for i = 1:my_game.n]

    # Big-M: maximum difference between two strategies
    V = []
    for p in 1:my_game.n
        aux_v = 0
        for j in 1:my_game.strat[p]
            for i in 1:my_game.strat[p]
                aux_j_i = sum(maximum(my_game.polymatrix[(p,k)][j,:] - my_game.polymatrix[(p,k)][i,:]) for k in 1:my_game.n)
                if aux_j_i > aux_v
                    aux_v = aux_j_i
                end
            end
        end
        push!(V, aux_v)
    end

    # Constraints
    for p in 1:my_game.n
        # add probability distribution constraint
        @constraint(model, sum(x[p]) == 1) # constrait 1 in the paper
        for j in 1:my_game.strat[p]
            profitj = sum(my_game.polymatrix[(p,i)][j,:]'*x[i] for i in 1:my_game.n) 
            @constraint(model, u[p][j] == profitj) # constraint 2 in the paper
            @constraint(model, v[p] >= u[p][j]) # constraint 3 in the paper
            @constraint(model, r[p][j] == v[p]-u[p][j]) # constraint 4 in the paper
            @constraint(model, f[p][j]>=r[p][j]/V[p])
            @constraint(model, f[p][j]>=y[p][j])
            @constraint(model, g[p][j] >= x[p][j])
            @constraint(model, g[p][j] >= 1-y[p][j])
        end

        # Objective function: minimize the probabilities and regret of the strategies with positive regret (combination of the formulations 2 and 3)
        @objective(model, Min, sum(sum(g[p][j]+f[p][j] for j in 1:my_game.strat[p]) for p in 1:my_game.n))
    end

    # Solve model
    optimize!(model)

    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), value.(v), [value.(x[p]) for p in 1:my_game.n] #, [value.(y[p]) for p in 1:my_game.n], model, y, v, x
    else
        println("Model is not solved and/or not feasible.")
        return "NaN","NaN","NaN" #, "NaN", model, y, v, x
    end
end

"""
    CorrelatedEquilibria(my_game)

    Method using to compute correlated equilibria. It just solves a linear program
"""
function CorrelatedEquilibria(my_game,supp = [])
    # supp: find CE restricted to a given support
    # start model
    model = Model(Gurobi.Optimizer)

    set_silent(model)
    set_time_limit_sec(model, 300.0)
    MOI.set(model, MOI.NumberOfThreads(), 1)

    # Variables
    # define probability distribution over all combinations of strategies
    indices = collect(product(ntuple(i -> 1:my_game.strat[i], length(my_game.strat))...)) # Cartesian product
    @variable(model, tau[indices], lower_bound = 0.0, upper_bound = 1.0)

    # add probability distribution constraint
    @constraint(model, sum(tau) == 1)

    M = Dict{Tuple{Int, Int}, Matrix{AffExpr}}()
    v = 0 # social welfare
    for p in 1:my_game.n
        for i in 1:my_game.n
            M[(p,i)] = zeros(AffExpr, my_game.strat[p], my_game.strat[i])
            M[(i,p)] = zeros(AffExpr, my_game.strat[i], my_game.strat[p])
            for ind in indices
                j,k = ind[p],ind[i]
                M[(p,i)][j,k] += tau[ind]
                M[(i,p)][k,j] += tau[ind]
            end
            v += sum(my_game.polymatrix[(p,i)].*M[(p,i)])
        end
        for j in 1:my_game.strat[p]
            for k in 1:my_game.strat[p]
                if j !=k # Calculate the expected profit for each pair of straties (correlated equilibrium constraint)
                    profitj = sum(my_game.polymatrix[(p,i)][j,:]'*M[(p,i)][j,:] for i in 1:my_game.n)
                    profitk = sum(my_game.polymatrix[(p,i)][k,:]'*M[(p,i)][j,:] for i in 1:my_game.n) # deviation
                    @constraint(model, profitj >= profitk)
                end
            end
        end
    end

    # restrict the CE to the support supp
    if supp != []
        for ind in indices
            if !(ind in supp) # ind not in supp
                @constraint(model, tau[ind] == 0)
            end
        end
    end


    # Objective function
    @objective(model, Max, v) 

    # Solve model
    optimize!(model)

    # each player expected utility
    M = Dict{Tuple{Int, Int}, Matrix{Float64}}()
    CE_u = zeros(my_game.n)
    for p in 1:my_game.n
        for i in 1:my_game.n
            M[(p,i)] = zeros(my_game.strat[p], my_game.strat[i])
            M[(i,p)] = zeros(my_game.strat[i], my_game.strat[p])
            for ind in indices
                j,k = ind[p],ind[i]
                M[(p,i)][j,k] += value.(tau[ind])
                M[(i,p)][k,j] += value.(tau[ind])
            end
            CE_u[p] +=  sum(my_game.polymatrix[(p,i)].*M[(p,i)])
        end
    end
    
    # Solution
    if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
        println("Model is solved and feasible and solving time is ",solve_time(model))
        return objective_value(model), [value.(tau[ind]) for ind in indices], CE_u, solve_time(model)
    else
        println("Model is not solved and/or not feasible, and solve time is ",solve_time(model))
        return "NaN", "NaN", "NaN", solve_time(model)
    end
end

# With correlated equilibria, we can only compute 1. If we want to compute all, then we need to find a constraint that allow us to have y=1 if and only if tau>0. It is not obvious how to do it.
# Below is a tentative code (some error in the callback), but it is not correct due to the lack of the condition above.

# function CorrelatedEquilibria_CallBack(my_game)
#     # start model
#     model = Model(Gurobi.Optimizer)
#     model[:callback_time] = 0.
#     model[:callback_calls] = 0.
#     model[:CE_mixed] = [] 
#     model[:CE_u] = []

#     set_silent(model)
#     set_time_limit_sec(model, 300.0)
#     MOI.set(model, MOI.NumberOfThreads(), 1)

#     # Variables
#     # define probability distribution for each player
#     indices = collect(product(ntuple(i -> 1:my_game.strat[i], length(my_game.strat))...)) # Cartesian product
#     @variable(model, tau[indices], lower_bound = 0.0, upper_bound = 1.0)
#     # define indicator variable to obtain support
#     @variable(model, y[indices], Bin) # to compute 1 CE, the y variables are irrelevant. For the rest of the CE we use no-good cuts with y.

#     # Constraints
#     # if tau>0 then y=1
#     @constraint(model, y .>= tau) # componentwise constraint
#     # add probability distribution constraint
#     @constraint(model, sum(tau) == 1)

#     M = Dict{Tuple{Int, Int}, Matrix{AffExpr}}()
#     v = 0 # social welfare
#     for p in 1:my_game.n
#         for i in 1:my_game.n
#             M[(p,i)] = zeros(AffExpr, my_game.strat[p], my_game.strat[i])
#             M[(i,p)] = zeros(AffExpr, my_game.strat[i], my_game.strat[p])
#             for ind in indices
#                 j,k = ind[p],ind[i]
#                 M[(p,i)][j,k] += tau[ind]
#                 M[(i,p)][k,j] += tau[ind]
#             end
#         end
#         for j in 1:my_game.strat[p]
#             for k in 1:my_game.strat[p]
#                 if j !=k # Calculate the expected profit for each pair of straties (correlated equilibrium constraint)
#                     profitj = sum(my_game.polymatrix[(p,i)][j,:]'*M[(p,i)][j,:] for i in 1:my_game.n)
#                     profitk = sum(my_game.polymatrix[(p,i)][k,:]'*M[(p,i)][j,:] for i in 1:my_game.n) # deviation
#                     @constraint(model, profitj >= profitk)
#                     v += profitj
#                 end
#             end
#         end
#     end

#     println("made it here: all constraints done")
#     # Objective function
#     # we can only optimize social welfare if we want to compute 1 equilibrium. Otherwise, we must minimize 
#     # @objective(model, Max, v) 

#     # Tell JuMP to use lazy constraints with the callback function
#     MOI.set(model, MOI.LazyConstraintCallback(), (cb_data) -> no_good_cut_callbackCE(cb_data,tau,y,model,my_game,indices))

#     # Solve model
#     optimize!(model)

#     # MISSINF to compute CE_u: see callback

#     # Solution
#     if termination_status(model) == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT
#         # println("Model is solved and feasible and solving time is ",solve_time(model))
#         return model[:CE_mixed], model[:CE_u], solve_time(model), model[:callback_time], model[:callback_calls]
#     else
#         # println("Model is not solved and/or not feasible, and solve time is ",solve_time(model))
#         return model[:CE_mixed], model[:CE_u], solve_time(model), model[:callback_time], model[:callback_calls]
#     end
# end

# function no_good_cut_callbackCE(cb_data,tau,y,model,my_game,indices)
#     status = callback_node_status(cb_data, model)
#     if status == MOI.CALLBACK_NODE_STATUS_INTEGER
#         push!(model[:CE_mixed],[callback_value.(cb_data,tau[ind]) for ind in indices])
#         # Compute individual expected profits
#         #push!(model[:CE_u],[ my_game.polymatrix*callback_value.(cb_data,tau[ind]) for p in 1:my_game.n]
#         model[:callback_time] += @elapsed begin # accumulated time to use call back
#             supp = [callback_value.(cb_data,tau[ind]) for ind in indices]
#             println("supp=",supp)
#             con = @build_constraint(sum(1-y[ind] for ind in indices if supp[ind]>=10.0^(-6))+sum(y[ind] for ind in indices if supp[ind]<10.0^(-6)) >=1) 
#             MOI.submit(model, MOI.LazyConstraint(cb_data), con)
#         end
#         model[:callback_calls] += 1 # count number of callbacks
#     end
# end

