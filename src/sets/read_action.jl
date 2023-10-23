function int(r::Real)
    return round(Int, r)
end


function read_action(H::Hyperrectangle)
    return int(H.center[end])
end


function read_action(H::IntervalBox)
    return int(mid(H.v[end]))
end


function read_action(X::Vector{<:TaylorModelN})
    return int(X[end].pol.coeffs[1])
end


function read_action(X::Vector{<:TaylorModel1})
    return int(X[end].pol.coeffs[1])
end


function read_action(X::ApproximatedSet)
    return read_action(X.Y)
end
