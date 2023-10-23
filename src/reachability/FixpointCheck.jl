struct FixpointCheck{T, D}
    tracker::T
    dimensions::D
    union::Bool
    learning_predicate
    split::Bool
    skip_redundancy_check::Bool
end

function FixpointCheck(; tracker=BoxTracker(), dimensions=nothing, union::Bool=false,
                         learning_predicate=(X, i) -> false, split::Bool=false,
                         skip_redundancy_check::Bool=false)
    @assert !isnothing(tracker)
    return FixpointCheck(tracker, dimensions, union, learning_predicate, split, skip_redundancy_check)
end


function wants_to_learn(F::FixpointCheck, X, i)
    return F.learning_predicate(X, i)
end


function learn!(F::FixpointCheck, X::ApproximatedSet; verbose::Bool=false)
    learn!(F, X.X; verbose=verbose)
end


function learn!(F::FixpointCheck, X::Union{LazySet, IntervalBox}; verbose::Bool=false)
    Y = isnothing(F.dimensions) ? X : project(X, F.dimensions)
    track!(F.tracker, Y; verbose=verbose)
end


function learn!(F::FixpointCheck, X::Vector)
    # ignore Taylor models
end


function check(F::FixpointCheck, X::ApproximatedSet)
    return check(F, X.X)
end


function check(F::FixpointCheck, X::LazySet)
    if F.skip_redundancy_check
        return false
    end

    X_π = isnothing(F.dimensions) ? X : project(X, F.dimensions)

    if F.union
        tracker = [UnionSetArray(F.tracker)]
    else
        tracker = F.tracker
    end
    for Y in tracker
        if X_π == Y
            continue
        elseif X_π ⊆ Y
            return true
        end
    end
    return false
end


function check(F::FixpointCheck, B::IntervalBox)
    H = convert(Hyperrectangle, B)
    return check(F, H)
end


# Taylor model
function check(F::FixpointCheck, X::Vector)
    Y = _overapproximate(X)
    return check(F, Y)
end


function split(F::FixpointCheck, H::Union{AbstractHyperrectangle, IntervalBox}; verbose::Bool=false)
    if F.skip_redundancy_check
        return H
    end

    if check(F, H)
        verbose && println("saved redundant leaf")
        return nothing
    end
    if !F.split
        return H
    end
    union = partition(H, F.tracker, F.dimensions)
    if !isnothing(F.dimensions)
        union = merge_dimensions(union, H, F.dimensions)
    end
    verbose && println("partition into $(length(union)) smaller sets")
    return UnionSetArray(union)
end
