# partition the first set into boxes that are not covered by any box in the vector
function partition(H::AbstractHyperrectangle, boxes::Vector{<:Hyperrectangle}, dimensions)
    diff = isnothing(dimensions) ? [H] : [project(H, dimensions)]
    for H2 in boxes
        diff_new = Vector{eltype(diff)}()
        for H1 in diff
            append!(diff_new, array(difference(H1, H2)))
        end
        diff = diff_new
    end
    return diff
end


function partition(B::IntervalBox, boxes::Vector{<:Hyperrectangle}, dimensions)
    return partition(convert(Hyperrectangle, B), boxes, dimensions)
end
