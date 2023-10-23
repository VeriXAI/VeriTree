struct ContSimProblem{F} <: AbstractSimProblem
    f::F
    Δt::Float64
    δ::Float64
end


function ContSimProblem(f, Δt)
    return ContSimProblem(f, Δt, Δt / 10.0)
end
