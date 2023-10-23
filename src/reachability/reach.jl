function reach(X0, DT::AbstractDecisionTree, P::AbstractReachProblem,
               abort; preprocess_fun=identity, preprocess_undo_fun=identity,
               fixpoint_check=FixpointCheck(), state_filter=nothing,
               verbose=false)
    total_post_time = 0.0
    used_termination_criterion = false
    reach_states_step = ApproximatedSet[]  # reachable states at time steps
    reach_states = abort(X0, 1) != no_abort ? [X0] : []  # reachable states
    queue_steps = Tuple{ApproximatedSet, Int}[]
    distribute_X0!(X0, P, reach_states_step, queue_steps)
    while !isempty(queue_steps)
        X_root, i = pop!(queue_steps)

        skip = false
        abort_res = abort(X_root, i)
        if abort_res != no_abort
            verbose && println("reached termination criterion")
            used_termination_criterion = abort_res == time_bound
            skip = true
        elseif i > 1 && check(fixpoint_check, X_root)
            verbose && println("found a fixpoint")
            skip = true
        end
        if skip
            continue
        end

        X_root_preprocessed = preprocess_fun(X_root)
        verbose && println("i = $i, |Q| = $(length(queue_steps))")
        queue_tree = Tuple{ApproximatedSet, AbstractDecisionTree}[]
        push!(queue_tree, (X_root_preprocessed, DT))
        leaves = Tuple{ApproximatedSet, Action}[]
        actions = Set{Action}()

        # propagate set down the decision tree (in-order traversal)
        while !isempty(queue_tree)
            X, dt_node = pop!(queue_tree)
            if isleaf(dt_node)
                # leaf node
                push!(leaves, (preprocess_undo_fun(X), dt_node.payload))
                push!(actions, dt_node.payload)
            else
                # inner node: continue along predicate
                pred = dt_node.predicate
                if X ⊆ pred
                    # traverse to left child
                    push!(queue_tree, (X, dt_node.left))
                elseif isdisjoint(X, pred)
                    # traverse to right child
                    push!(queue_tree, (X, dt_node.right))
                else
                    # split and traverse to both children
                    X1, X2 = split(P.alg, X, pred, fixpoint_check; check_special_cases=false)
                    push_splits!(queue_tree, X2, dt_node.right)
                    push_splits!(queue_tree, X1, dt_node.left)
                end
            end
        end

        leaves = merge_leaves!(leaves, actions, X_root, i, fixpoint_check; verbose=verbose)

        for leaf in leaves
            # process set arrived at the leaf node of the decision tree
            X, action = leaf
            sol, Y, comp_time = post(P, X, action)
            total_post_time += comp_time
            push!(reach_states_step, Y)
            push!(reach_states, sol)
            Z = isnothing(state_filter) ? Y : state_filter(Y)
            if !isempty(Z)
                # continue with the next time step
                push!(queue_steps, (Z, i+1))
            end
        end
    end

    if used_termination_criterion
        println("result does not generalize to unbounded time")
    else
        println("result generalizes to unbounded time")
    end
    println("time spent in post computation of dynamical system: $total_post_time seconds")
    return reach_states_step, reach_states
end


function distribute_X0!(X0::LazySet, P, reach_states_step, queue_steps)
    push!(reach_states_step, ApproximatedSet(X0, 0.0))
    push!(queue_steps, (ApproximatedSet(initialize(P, X0), 0.0), 1))
end


function distribute_X0!(X0::UnionSetArray, P, reach_states_step, queue_steps)
    for Y0 in array(X0)
        distribute_X0!(Y0, P, reach_states_step, queue_steps)
    end
end


function merge_leaves!(leaves, actions, X_root, i, fixpoint_check; verbose)
    if length(actions) == 1
        verbose && println("only one action")
        X_app = _overapproximate(X_root.Y)
        if wants_to_learn(fixpoint_check, X_app, i)
            verbose && println("approximating anyway")
            learn!(fixpoint_check, X_app; verbose=verbose)
            leaves = [(ApproximatedSet(X_app, X_root.time), first(actions))]
        else
            leaves = [(X_root, first(actions))]
        end
    else
        verbose && println("$(length(actions)) actions, exploring all")
        learn!(fixpoint_check, _overapproximate(X_root.Y); verbose=verbose)
    end
    return leaves
end


function push_splits!(queue_tree, X::Nothing, tree)
    # ignore
end


function push_splits!(queue_tree,
                      X::ApproximatedSet{T1, <:Union{<:AbstractHyperrectangle, <:IntervalBox}},
                      tree) where T1
    push!(queue_tree, (X, tree))
end


function push_splits!(queue_tree, X::ApproximatedSet{T1, <:UnionSetArray}, tree) where T1
    for Y in array(X.Y)
        push_splits!(queue_tree, ApproximatedSet(Y, X.time), tree)
    end
end
