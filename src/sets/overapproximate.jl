function _overapproximate(X::LazySet)
    return box_approximation(X)
end

function _overapproximate(X::IntervalBox)
    return convert(Hyperrectangle, X)
end

function _overapproximate(X::Vector{<:TaylorModelN})
    return box_approximation(set(overapproximate(X, Zonotope)))
end

function _overapproximate(X::Vector{<:TaylorModel1})
    return set(overapproximate(TaylorModelReachSet(X, 0..0), Hyperrectangle))
end
