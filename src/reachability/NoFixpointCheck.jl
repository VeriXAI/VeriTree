struct NoFixpointCheck
end


function wants_to_learn(::NoFixpointCheck, X, i)
    return false
end


function learn!(::NoFixpointCheck, X; verbose::Bool=false)
    # ignored
end


function check(::NoFixpointCheck, X)
    return false
end


function split(::NoFixpointCheck, X; verbose::Bool=false)
    return X
end
