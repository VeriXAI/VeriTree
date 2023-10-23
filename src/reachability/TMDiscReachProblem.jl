struct TMDiscReachProblem{F, A} <: AbstractReachProblem
    f::F
    alg::A
end


function TMDiscReachProblem(f)
    return TMDiscReachProblem(f, TMJets21b())
end
