# note: the conditions are exactly inverted for `⊆`
@commutative function isdisjoint(B::OpenHyperrectangle, H::OrthogonalHalfSpace)
    i = H.i
    b = H.b
    if H.lt
        # x <= b
        lo = low(B, i)
        return b < lo || (b == lo && (H.open || B.lo_open[i]))
    else
        # x >= b
        hi = high(B, i)
        return b > hi || (b == hi && (H.open || B.hi_open[i]))
    end
end


@commutative function isdisjoint(B::AbstractHyperrectangle, H::OrthogonalHalfSpace)
    i = H.i
    b = H.b
    if H.lt
        # x <= b
        lo = low(B, i)
        return b < lo || (b == lo && H.open)
    else
        # x >= b
        hi = high(B, i)
        return b > hi || (b == hi && H.open)
    end
end

@commutative function isdisjoint(B::IntervalBox, H::OrthogonalHalfSpace)
    return isdisjoint(convert(Hyperrectangle, B), H)
end

@commutative function isdisjoint(X::ApproximatedSet, Y)
    return isdisjoint(X.Y, Y)  # perform disjointness check for the approximating set
end
