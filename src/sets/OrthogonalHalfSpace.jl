struct OrthogonalHalfSpace{N<:Real}
    i::Int
    n::Int
    b::N
    lt::Bool
    open::Bool
end


function dim(H::OrthogonalHalfSpace)
    return H.dim
end


function complement(H::OrthogonalHalfSpace)
    return OrthogonalHalfSpace(H.i, H.n, H.b, !H.lt, !H.open)
end


# conversion to HalfSpace
function HalfSpace(H::OrthogonalHalfSpace{N}) where {N}
    if H.lt
        o = one(N)
        b = H.b
    else
        o = -one(N)
        b = -H.b
    end
    return HalfSpace(ReachabilityAnalysis.LazySets.SingleEntryVector(H.i, H.n, o), b)
end


function ∈(x::AbstractVector, H::OrthogonalHalfSpace)
    xi = x[H.i]
    if H.open && H.lt
        return xi < H.b
    elseif H.open && !H.lt
        return xi > H.b
    elseif !H.open && H.lt
        return xi <= H.b
    else
        return xi >= H.b
    end
end


@recipe function _plot(H::OrthogonalHalfSpace)
    return HalfSpace(H)
end
