function initialize(P::TMContReachProblem, X::Union{Vector, IntervalBox})
    return X
end


function initialize(P::TMContReachProblem, X::LazySet)
    return ReachabilityAnalysis._initialize(X, P.alg.orderT, P.alg.orderQ)
end


function post(P::TMContReachProblem, X::ApproximatedSet, action::Action)
    n = dim(X)
    Xt = apply(action, X_time(X))
    ivp = IVP(BlackBoxContinuousSystem(P.f, n), Xt)
    t0 = X.time
    t1 = t0 + P.Δt
    comp_time =  @elapsed sol = ReachabilityAnalysis.post(P.alg, ivp, t0 .. t1)
    sol_end = sol[end]  # reach set (in time) at the next time step
    Y_precise_time = evaluate(set(sol_end), sup(domain(sol_end)))  # set without time
    # convert polynomials back to a TaylorModelN vector
    rem = remainder(sol_end)
    zB = zeroBox(n)
    sB = symBox(n)
    Y_precise = [fp_rpa(TaylorModelN(Y_precise_time[j], rem[j], zB, sB)) for j in 1:n]
    # wrap TaylorModelN vector back into TaylorModel1 vector in time
    Y_precise_time = [TaylorModel1(Taylor1(polynomial(Y_precise[j]), P.alg.orderT), zeroI, zeroI, zeroI) for j in 1:n]
    # box approximation
    Y_approx = _overapproximate(Y_precise)
    return sol, ApproximatedSet(Y_precise, Y_approx, Y_precise_time, t1), comp_time
end
