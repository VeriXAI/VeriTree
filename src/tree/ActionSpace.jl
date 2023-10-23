struct ActionSpace{ST}
    action2sets::Dict{Action, Vector{ST}}
end


function ActionSpace(DT::AbstractUnidimensionalBinaryDecisionTree,
                     state_space;
                     ST=OpenHyperrectangle)
    a2sets = Dict{Action, Vector{ST}}()

    queue = Tuple{ST, AbstractUnidimensionalBinaryDecisionTree}[]
    push!(queue, (state_space, DT))
    while !isempty(queue)
        X, T = pop!(queue)
        if isleaf(T)
            # leaf node: store the set for this action
            action = T.payload
            sets = get(a2sets, action, nothing)
            if isnothing(sets)
                sets = Vector{ST}()
                a2sets[action] = sets
            end
            push!(sets, X)
        else
            # intermediate node: split along predicate
            P = T.predicate
            if X ⊆ P
                # traverse to left child
                push!(queue, (X, T.left))
            elseif isdisjoint(X, P)
                # traverse to right child
                push!(queue, (X, T.right))
            else
                # split and traverse to both children
                B1, B2 = split(X, P; check_special_cases=false)
                push!(queue, (B2, T.right))
                push!(queue, (B1, T.left))
            end
        end
    end

    return ActionSpace(a2sets)
end


@recipe function _plot(AS::ActionSpace; vars=nothing, label=nothing, action2sets=sort(AS.action2sets))
    for (i, (ai, Xi)) in enumerate(action2sets)
        if isempty(Xi)
            continue
        end
        @series begin
            color --> i
            if isnothing(label)
                label := ai.name
            end
            if isnothing(vars)
                Xi[1]
            else
                project(Xi[1], vars)
            end
        end
        @series begin
            color --> i
            label := ""
            if isnothing(vars)
                Xi[2:end]
            else
                [project(Xij, vars) for Xij in Xi[2:end]]
            end
        end
    end
end
