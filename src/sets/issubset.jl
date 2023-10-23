# note: the conditions are exactly inverted for `isdisjoint`
function ⊆(B::OpenHyperrectangle, H::OrthogonalHalfSpace)
    i = H.i
    b = H.b
    if H.lt
        # x <= b
        hi = high(B, i)
        return b > hi || (b == hi && (!H.open || B.hi_open[i]))
    else
        # x >= b
        lo = low(B, i)
        return b < lo || (b == lo && (!H.open || B.lo_open[i]))
    end
end


function ⊆(B::AbstractHyperrectangle, H::OrthogonalHalfSpace)
    i = H.i
    b = H.b
    if H.lt
        # x <= b
        hi = high(B, i)
        return b > hi || (b == hi && !H.open)
    else
        # x >= b
        lo = low(B, i)
        return b < lo || (b == lo && !H.open)
    end
end


function ⊆(B::IntervalBox, H::OrthogonalHalfSpace)
    return convert(Hyperrectangle, B) ⊆ H
end


function ⊆(X::ApproximatedSet, Y; dimensions=nothing)
    # sufficient subset check for the approximating set
    if isnothing(dimensions) || !(Y isa Vector{<:TaylorModelN}) || !(Y isa Vector{<:TaylorModel1})
        return X.Y ⊆ Y
    end
    return ⊆(X.Y, Y; dimensions=dimensions)
end
