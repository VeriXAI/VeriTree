struct TMContReachProblem{F, A} <: AbstractReachProblem
    f::F
    Δt::Float64
    alg::A
end


function TMContReachProblem(f, Δt)
    return TMContReachProblem(f, Δt, TMJets21b())
end
