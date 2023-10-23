struct UnidimensionalBinaryDecisionTree{N, L, R} <: AbstractDecisionTree
    predicate::OrthogonalHalfSpace{N}
    left::L
    right::R
end


function isleaf(::UnidimensionalBinaryDecisionTree)
    return false
end


function get_leaf(DT::UnidimensionalBinaryDecisionTree, x::AbstractVector)
    branch = x ∈ DT.predicate ? DT.left : DT.right
    return get_leaf(branch, x)
end
