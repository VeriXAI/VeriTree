function initialize(P::TMDiscReachProblem, X::Union{Vector, IntervalBox})
    return X
end


function initialize(P::TMDiscReachProblem{F, <:IntervalArith}, X::LazySet) where F
    return convert(IntervalBox, box_approximation(X))
end


function initialize(P::TMDiscReachProblem{F, <:TMJets21b}, X::LazySet) where F
    return ReachabilityAnalysis._initialize(X, P.alg.orderT, P.alg.orderQ)
end


function post(P::TMDiscReachProblem, X::ApproximatedSet, action::Action)
    Xt = apply(action, X.X)

    comp_time =  @elapsed sol = _post_discrete(P.alg, P, Xt)
    Y_precise = sol

    for i in eachindex(Y_precise)
        Y_precise[i] = TaylorModel1(Y_precise[i].pol, 0..0, Y_precise[i].x0, Y_precise[i].dom)
    end

    Y_approx = _overapproximate(Y_precise)

    return sol, ApproximatedSet(Y_precise, Y_approx, nothing, X.time + 1), comp_time
end


function _post_discrete(::IntervalArith, P, Xt)
    return P.f(Xt)
end


function _post_discrete(::TMJets21b, P, Xt::Vector)
    return P.f(Xt)
end


function _post_discrete(alg::TMJets21b, P, Xt::AbstractHyperrectangle)
    Yt = convert(TaylorModelReachSet, Xt;
                 orderQ=alg.orderQ, orderT=alg.orderT).X
    return _post_discrete(alg, P, Yt)
end


function _post_discrete(alg::TMJets21b, P, Xt::IntervalBox)
    return _post_discrete(alg, P, convert(Hyperrectangle, Xt))
end
