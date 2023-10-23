struct ApproximatedSet{T1, T2, T3, N<:Real}
    X::T1  # approximated set
    Y::T2  # approximating set
    X_time::T3
    time::N
end


function ApproximatedSet(X, time::Real)
    return ApproximatedSet(X, X, nothing, time)
end


# convenience constructor to avoid filling with `nothing`
function ApproximatedSet(X::Nothing, time::Real=0.0)
    return nothing
end


function X_time(X::ApproximatedSet)
    if !isnothing(X.X_time)
        return X.X_time
    end
    return X.X
end


function dim(X::ApproximatedSet)
    return dim(X.Y)
end


function isempty(X::ApproximatedSet)
    if isempty(X.Y)
        return true
    end
    return isempty(X.X)
end

function isbox(::ApproximatedSet)
    return false
end

function isbox(::ApproximatedSet{<:AbstractHyperrectangle})
    return true
end

function isbox(::ApproximatedSet{<:IntervalBox})
    return true
end
