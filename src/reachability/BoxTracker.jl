struct BoxTracker
    sets::Vector{Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}}
end


function BoxTracker()
    return BoxTracker(Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}[])
end


function iterate(T::BoxTracker)
    return iterate(T.sets)
end


function iterate(T::BoxTracker, i::Int)
    return iterate(T.sets, i)
end


function track!(T::BoxTracker, X; verbose::Bool=false)
    # ignored
end


function track!(T::BoxTracker, X::Hyperrectangle; verbose::Bool=false)
    push!(T.sets, X)
    verbose && println("learned $X")
end


function track!(T::BoxTracker, B::IntervalBox; verbose::Bool=false)
    track!(T, convert(Hyperrectangle, B))
end


# partition the first set into boxes that are not covered by any set in the tracker
function partition(H::Union{AbstractHyperrectangle, IntervalBox}, tracker::BoxTracker, dimensions)
    return partition(H, tracker.sets, dimensions)
end
