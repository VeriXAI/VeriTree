function merge_dimensions(union::Vector{<:AbstractHyperrectangle}, H, dimensions)
    n = dim(H)
    union_merged = Vector{Hyperrectangle{Float64, Vector{Float64}, Vector{Float64}}}(undef, length(union))
    for (k, X) in enumerate(union)
        j = 1
        l = zeros(n)
        h = zeros(n)
        for i in 1:n
            if i in dimensions
                l[i] = low(X, j)
                h[i] = high(X, j)
                j += 1
            else
                l[i] = low(H, i)
                h[i] = high(H, i)
            end
        end
        union_merged[k] = Hyperrectangle(low=l, high=h)
    end
    return union_merged
end
