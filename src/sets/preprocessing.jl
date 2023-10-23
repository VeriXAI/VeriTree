function preprocessing_given_mapping(mapping)
    return X -> preprocess(X, mapping)
end

function preprocess(x::AbstractVector, preprocessing_mapping)
    return preprocessing_mapping(x)
end

function preprocess(X::ApproximatedSet, preprocessing_mapping)
    # use the precise set (typically a Taylor model except for the first step)
    X_preprocessed = _preprocess(X.X, preprocessing_mapping)
    X_preprocessed_approximated = _overapproximate(X_preprocessed)
    return ApproximatedSet(X_preprocessed, X_preprocessed_approximated,
                           X.X_time, X.time)
end

function _preprocess(X::Vector{<:TaylorModelN}, preprocessing_mapping)
    intervals = map(p -> evaluate(p, p.dom), preprocessing_mapping(X))
    return convert(Hyperrectangle, IntervalBox(intervals))
end

function _preprocess(X::IntervalBox, preprocessing_mapping)
    intervals = preprocessing_mapping(X)  # uses interval arithmetic
    return convert(Hyperrectangle, IntervalBox(intervals))
end

function _preprocess(X::LazySet, preprocessing_mapping)
    B = convert(IntervalBox, box_approximation(X))
    return _preprocess(B, preprocessing_mapping)
end
