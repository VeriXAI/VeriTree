struct OpenHyperrectangle{N, NVT<:AbstractVector{N}, BVT<:AbstractVector{Bool}}
    lo::NVT
    hi::NVT
    lo_open::BVT
    hi_open::BVT

    function OpenHyperrectangle(lo::NVT, hi::NVT, lo_open::BVT, hi_open::BVT
                               ) where {N, NVT<:AbstractVector{N}, BVT<:AbstractVector{Bool}}
        @assert length(lo) == length(hi) == length(lo_open) == length(hi_open)
        return new{N, NVT, BVT}(lo, hi, lo_open, hi_open)
    end
end


# conversion from Hyperrectangle
function OpenHyperrectangle(H::AbstractHyperrectangle)
    n = dim(H)
    lo = low(H)
    hi = high(H)
    lo_open = falses(n)
    hi_open = falses(n)
    return OpenHyperrectangle(lo, hi, lo_open, hi_open)
end


# conversion to Hyperrectangle
function Hyperrectangle(OH::OpenHyperrectangle)
    return Hyperrectangle(low=OH.lo, high=OH.hi)
end


function convert(::Type{<:Hyperrectangle}, OH::OpenHyperrectangle)
    return Hyperrectangle(OH)
end


function convert(::Type{<:OpenHyperrectangle}, H::AbstractHyperrectangle)
    return OpenHyperrectangle(H)
end


function dim(B::OpenHyperrectangle)
    return length(B.lo)
end


# note: may return a value outside the set if open
function low(B::OpenHyperrectangle, i::Int)
    return B.lo[i]
end


# note: may return a value outside the set if open
function high(B::OpenHyperrectangle, i::Int)
    return B.hi[i]
end


function isempty(B::OpenHyperrectangle)
    @inbounds for i in 1:dim(B)
        l, h = B.lo[i], B.hi[i]
        if l > h || (l == h && (B.lo_open[i] || B.hi_open[i]))
            return true
        end
    end
    return false
end


function translate(B::OpenHyperrectangle, v::AbstractVector)
    return OpenHyperrectangle(B.lo + v, B.hi + v, B.lo_open, B.hi_open)
end


function linear_map(A::AbstractMatrix, B::OpenHyperrectangle)
    return linear_map(A, Hyperrectangle(B))
end


@recipe function _plot(B::OpenHyperrectangle)
    return Hyperrectangle(B)
end


@recipe function _plot(V::AbstractVector{<:OpenHyperrectangle})
    return [Hyperrectangle(B) for B in V]
end
