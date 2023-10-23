# default: ignore the algorithm
function split(alg, X, P, fixpoint_check; check_special_cases::Bool=true)
    return split(X, P; fixpoint_check=fixpoint_check, check_special_cases=check_special_cases)
end


function split(alg::IntervalArith, X, P, fixpoint_check; check_special_cases::Bool=true)
    X1, X2 = split(X, P; fixpoint_check=fixpoint_check, check_special_cases=check_special_cases)
    if !isnothing(X1)
        X1 = convert(IntervalBox, X1)
    end
    if !isnothing(X2)
        X2 = convert(IntervalBox, X2)
    end
    return (X1, X2)
end


function split(alg::IntervalArith, X::ApproximatedSet, P, fixpoint_check;
               check_special_cases::Bool=true)
    X1, X2 = split(alg, X.Y, P, fixpoint_check; check_special_cases=check_special_cases)
    return ApproximatedSet(X1, X.time), ApproximatedSet(X2, X.time)
end


function split(B::OpenHyperrectangle{N}, H::OrthogonalHalfSpace;
               fixpoint_check=NoFixpointCheck(), check_special_cases::Bool=true,
               verbose::Bool=false) where N
    if check_special_cases
        if isdisjoint(B, H)
            return EmptySet{N}(dim(B))
        elseif B ⊆ H
            return B
        end
    end

    B1 = deepcopy(B)
    B2 = deepcopy(B)
    i = H.i
    b = H.b
    B1.hi[i] = b
    B2.lo[i] = b
    if H.lt
        # x <= b
        B1.hi_open[i] = H.open
        B2.lo_open[i] = !H.open
    else
        # x >= b
        B1.hi_open[i] = !H.open
        B2.lo_open[i] = H.open
    end
    B1 = split(fixpoint_check, B1; verbose=verbose)
    B2 = split(fixpoint_check, B2; verbose=verbose)
    return B1, B2
end


function split(B::AbstractHyperrectangle{N}, H::OrthogonalHalfSpace;
               fixpoint_check=NoFixpointCheck(), check_special_cases::Bool=true,
               verbose::Bool=false) where N
    if check_special_cases
        if isdisjoint(B, H)
            return EmptySet{N}(dim(B))
        elseif B ⊆ H
            return B
        end
    end

    i = H.i
    b = H.b
    l1 = _vector(low(B))
    h1 = _vector(high(B))
    l2 = deepcopy(l1)
    h2 = deepcopy(h1)
    h1[i] = b
    l2[i] = b
    B1 = Hyperrectangle(low=l1, high=h1)
    B2 = Hyperrectangle(low=l2, high=h2)
    B1 = split(fixpoint_check, B1; verbose=verbose)
    B2 = split(fixpoint_check, B2; verbose=verbose)
    return B1, B2
end


function split(B::IntervalBox, H::OrthogonalHalfSpace{N};
               fixpoint_check=NoFixpointCheck(), check_special_cases::Bool=true,
               verbose::Bool=false) where N
    if check_special_cases
        if isdisjoint(B, H)
            return EmptySet{N}(dim(B))
        elseif B ⊆ H
            return B
        end
    end

    i = H.i
    b = H.b
    itv1 = inf(B.v[i]) .. b
    itv2 = b .. sup(B.v[i])
    B1 = IntervalBox(vcat(B.v[1:i-1], itv1, B.v[i+1:end]))
    B2 = IntervalBox(vcat(B.v[1:i-1], itv2, B.v[i+1:end]))
    B1 = split(fixpoint_check, B1; verbose=verbose)
    B2 = split(fixpoint_check, B2; verbose=verbose)
    return B1, B2
end


function split(X::ApproximatedSet, P;
               fixpoint_check=NoFixpointCheck(), check_special_cases::Bool=true,
               verbose::Bool=false)
    # perform split for the approximating set
    X1, X2 = split(X.Y, P; fixpoint_check=fixpoint_check,
                   check_special_cases=check_special_cases, verbose=verbose)
    return ApproximatedSet(X1, X.time), ApproximatedSet(X2, X.time)
end


# convert to a modifiable vector
function _vector(v::AbstractVector)
    return v
end


function _vector(v::IntervalArithmetic.StaticArrays.SVector)
    return Vector(v)
end
