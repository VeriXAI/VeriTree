struct Leaf{P} <:AbstractDecisionTree
    payload::P
end


function isleaf(::Leaf)
    return true
end


function get_leaf(L::Leaf, x)
    return L
end
