function simulate(x0::AbstractVector, DT::AbstractDecisionTree,
                  P::AbstractSimProblem, abort;
                  preprocess_fun=identity,
                  verbose=false)
    x = initialize(P, x0)
    reach_states = [x]
    i = 1
    while true
        if abort(x, i) != no_abort
            verbose && println("reached termination criterion")
            break
        end
        px = preprocess_fun(_state(x))
        action = get_leaf(DT, px).payload
        xs = post(P, x, action)
        append!(reach_states, xs)
        x = xs[end]
        i += 1
    end
    return reach_states
end


function simulate(X0::AbstractVector{<:AbstractVector},
                  DT::AbstractDecisionTree, P::AbstractSimProblem, abort;
                  preprocess_fun=identity, verbose=false)
    states = []
    for (i, x0) in enumerate(X0)
        x = simulate(x0, DT, P, abort; preprocess_fun=preprocess_fun,
                     verbose=verbose)
        push!(states, x)
    end
    return states
end


function simulate(X0::LazySet, DT::AbstractDecisionTree, P::AbstractSimProblem,
                  abort; preprocess_fun=identity, count::Int=10,
                  include_vertices::Bool=false, verbose=false)
    samples = sample(X0, count; include_vertices=include_vertices)
    return simulate(samples, DT, P, abort; preprocess_fun=preprocess_fun,
                     verbose=verbose)
end


function simulate(X0s::UnionSetArray, DT::AbstractDecisionTree, P::AbstractSimProblem,
                  abort; preprocess_fun=identity, count::Int=10,
                  include_vertices::Bool=false, verbose=false)
    n_average = div(count, length(X0s))
    samples = Vector{Vector{Float64}}()
    for X0 in array(X0s)
        append!(samples, sample(X0, n_average; include_vertices=include_vertices))
    end
    return simulate(samples, DT, P, abort; preprocess_fun=preprocess_fun,
                     verbose=verbose)
end

## discrete

function initialize(P::DiscSimProblem, x::AbstractVector)
    return x
end


function _state(x::AbstractVector)
    return x
end


function post(P::DiscSimProblem, x::AbstractVector, action::Action)
    y = apply(action, x)
    z = P.f(y)
    return [z]
end

## continuous


function initialize(P::ContSimProblem, x::AbstractVector)
    return x = (0.0, x)
end

function _state(tx::Tuple)
    return tx[2]
end


function post(P::ContSimProblem, tx::Tuple, action::Action)
    t, x = tx
    y = apply(action, x)
    zs = rk4(P.f, y, P.δ, t + P.Δt; t0=t)
    return zs[2:end]
end


# Euler method with arbitrary but fixed steps
function euler(f!, x0, steps::AbstractVector; t0=0.0)
    x = x0
    results = Vector{Tuple{typeof(t0), typeof(x0)}}(undef, length(steps))
    dx = copy(x)
    p = nothing  # ignore parameters
    for (i, ti) in enumerate(steps)
        if ti > 0
            δ = ti - t0
            f!(dx, x, p, t0)  # in-place update of dx
            x += δ * dx
            t0 = ti
        end
        results[i] = (t0, copy(x))
    end
    return results
end

# Euler method with a step bound
function euler(f!, x0, δ, k::Integer; t0=0.0)
    r = range(t0, step=δ, length=(k + 1))
    return euler(f!, x0, r; t0=t0)
end

# Euler method with a time bound
function euler(f!, x0, δ, T::Float64; t0=0.0)
    r = range(t0, stop=T, step=δ)
    return euler(f!, x0, r; t0=t0)
end


# Runge-Kutta method with arbitrary but fixed steps
function rk4(f!, x0, steps::AbstractVector; t0=0.0)
    x = x0
    results = Vector{Tuple{typeof(t0), typeof(x0)}}(undef, length(steps))
    k1 = copy(x)
    k2 = copy(x)
    k3 = copy(x)
    k4 = copy(x)
    p = nothing  # ignore parameters
    for (i, ti) in enumerate(steps)
        if ti > 0
            δ = ti - t0
            δ_half = δ / 2
            f!(k1, x, p, t0)
            f!(k2, x + k1 * δ_half, p, t0 + δ_half)
            f!(k3, x + k2 * δ_half, p, t0 + δ_half)
            f!(k4, x + k3 * δ, p, t0 + δ)
            x += 1/6 * δ * (k1 + 2 * k2 + 2 * k3 + k4)
            t0 = ti
        end
        results[i] = (t0, copy(x))
    end
    return results
end


# Runge-Kutta method with a step bound
function rk4(f!, x0, δ, k::Integer; t0=0.0)
    r = range(t0, step=δ, length=(k + 1))
    return rk4(f!, x0, r; t0=t0)
end


# Runge-Kutta method with a time bound
function rk4(f!, x0, δ, T::Float64; t0=0.0)
    r = range(t0, stop=T, step=δ)
    return rk4(f!, x0, r; t0=t0)
end
