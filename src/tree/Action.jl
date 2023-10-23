struct Action
    name::String
    value::Int
end


function Base.isless(A1::Action, A2::Action)
    return A1.name < A2.name
end


# The task of `apply(A, X)` is to update the control mode of `X` to `A.value`.
# By convention, the control mode is stored in the last dimension.


function apply(A::Action, X)
    return apply!(A, deepcopy(X))
end


function apply!(A::Action, x::Vector)
    x[end] = A.value
    return x
end


function apply(A::Action, H::Hyperrectangle{N, SV, SV}) where {N,
                                    SV<:IntervalArithmetic.StaticArrays.SVector}
    H2 = Hyperrectangle(_vector(center(H)), _vector(radius_hyperrectangle(H)))
    return apply!(A, H2)
end


function apply!(A::Action, H::Hyperrectangle)
    H.center[end] = A.value
    return H
end


function apply!(A::Action, H::IntervalBox)
    n = dim(H)
    H2 = IntervalBox(vcat(H.v[1:n-1], A.value .. A.value))
    return H2
end


function apply!(A::Action, X::Vector{<:TaylorModelN})
    # update the 0-order coefficient for the input variable
    X[end].pol.coeffs[1] = A.value
    return X
end


function apply!(A::Action, X::Vector{<:TaylorModel1})
    # update the 0-order coefficient for the input variable
    X[end].pol.coeffs[1] = A.value
    return X
end
