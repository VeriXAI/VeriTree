function project(B::IntervalBox, vars::AbstractVector{Int})
    return IntervalBox(getindex.(Ref(B), vars))
end
