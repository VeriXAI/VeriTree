@commutative function intersection(B::OpenHyperrectangle{N},
                                   H::OrthogonalHalfSpace;
                                   check_special_cases::Bool=true) where {N}
    if check_special_cases
        if isdisjoint(B, H)
            return EmptySet{N}(dim(B))
        elseif B ⊆ H
            return B
        end
    end

    B1 = deepcopy(B)
    i = H.i
    b = H.b
    if H.lt
        # x <= b
        B1.hi[i] = b
        B1.hi_open[i] = H.open
    else
        # x >= b
        B1.lo[i] = b
        B1.lo_open[i] = H.open
    end
    return B1
end


function intersection(B1::OpenHyperrectangle{N}, B2::OpenHyperrectangle) where {N}
    B3 = deepcopy(B1)
    @inbounds for i in 1:dim(B1)
        l1 = B1.lo[i]
        l2 = B2.lo[i]
        if l1 < l2 || (l1 == l2 && (B1.lo_open[i] || B2.lo_open[i]))
            B3.lo[i] = l2
            B3.lo_open[i] = B2.lo_open[i]
        else
            B3.lo[i] = l1
            B3.lo_open[i] = B1.lo_open[i]
        end

        h1 = B1.hi[i]
        h2 = B2.hi[i]
        if h1 > h2 || (h1 == h2 && (B1.hi_open[i] || B2.hi_open[i]))
            B3.hi[i] = h2
            B3.hi_open[i] = B2.hi_open[i]
        else
            B3.hi[i] = h1
            B3.hi_open[i] = B1.hi_open[i]
        end

        if B3.lo[i] > B3.hi[i] || (B3.lo[i] == B3.hi[i] && (B3.lo_open[i] || B3.hi_open[i]))
            return EmptySet{N}()
        end
    end
    return B3
end


@commutative function intersection(B1::OpenHyperrectangle, B2::AbstractHyperrectangle)
    return intersection(B1, OpenHyperrectangle(B2))
end
