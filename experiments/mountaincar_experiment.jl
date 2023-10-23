using DTCSR
using DTCSR: ApproximatedSet

# fix random seed
import Random
Random.seed!(1234)

# --- problem definition ---

gravity = 0.0025
force = 0.001
clip_x = [-1.2, 0.6]
clip_v = [-0.07, 0.07]

function mountaincar(s::AbstractVector)
    x = s[1]  # car position
    v = s[2]  # car velocity
    action = s[3]  # action

    v = v + (action - 1) * force - cos(3 * x) * gravity
    x = x + v

    return [x, v, action]
end

function mountaincar(s::IntervalBox)
    x = s[1]  # car position
    v = s[2]  # car velocity
    action = s[3]  # action

    v = v + (action - 1) * force - cos(3 * x) * gravity
    x = x + v

    return IntervalBox([x, v, action])
end

left = Action("left", 0)
none = Action("none", 1)
right = Action("right", 2)

DT_depth_3 = Node(
    OrthogonalHalfSpace(2, 2, 0.0, true, false),
    Node(
        OrthogonalHalfSpace(2, 2, -0.004, true, false),
        Leaf(left),
        Node(
            OrthogonalHalfSpace(1, 2, 0.391, true, false),
            Leaf(right),
            Leaf(left)
        )
    ),
    Node(
        OrthogonalHalfSpace(2, 2, -0.011, true, false),
        Node(
            OrthogonalHalfSpace(1, 2, -0.37, true, false),
            Leaf(right),
            Leaf(none)
        ),
        Leaf(right)
    )
)

DT = DT_depth_3

Ps = DiscSimProblem(mountaincar)
Pr = TMDiscReachProblem(mountaincar, IntervalArith())
Pr = TMDiscReachProblem(mountaincar, TMJets21b())

X0 = Hyperrectangle([-0.475, 0, 0], [0.03, 0, 0])

# predicate for aborting

kmax = 200

function abort_mountaincar(x::AbstractVector, i)
    if i > kmax
        return time_bound
    elseif x[1] >= 0.5
        return other_condition
    else
        return no_abort
    end
end

function abort_mountaincar(X, i)
    if i > kmax
        return time_bound
    else
        return no_abort
    end
end

function abort_mountaincar(X::ApproximatedSet{T1, <:IntervalBox}, i) where T1
    if i > kmax
        return time_bound
    elseif X.Y[1] >= 0.5
        return other_condition
    else
        return no_abort
    end
end

function abort_mountaincar(X::ApproximatedSet{T1, <:Hyperrectangle}, i) where T1
    if i > kmax
        return time_bound
    elseif low(X.Y, 1) >= 0.5
        return other_condition
    else
        return no_abort
    end
end

# predicate for aborting with clipping
function abort_mountaincar_clipped(x::AbstractVector, i)
    res = abort_mountaincar(x, i)
    if res != no_abort
        return res
    end
    if x[1] ∉ Interval(clip_x) || x[2] ∉ Interval(clip_v)
        return other_condition
    else
        return no_abort
    end
end

verbose = false
benchmark = true

# --- simulation ---

sim_mountaincar = simulate(X0, DT, Ps, abort_mountaincar_clipped, count=1000);

# --- reachability ---

fixpoint_check = FixpointCheck(split=true)
if benchmark
    # warm-up run for precompilation
    println("# Warm-up run:")
    reach(X0, DT, Pr, abort_mountaincar, fixpoint_check=deepcopy(fixpoint_check), verbose=false);
end
println("# Verification run:")
@time reach_step_mountaincar, reach_mountaincar =
    reach(X0, DT, Pr, abort_mountaincar, fixpoint_check=deepcopy(fixpoint_check), verbose=verbose);

# --- plotting ---

import Plots
using Plots: plot, plot!, hline!, vline!, colormap, savefig, ylims!, scatter!, font
using Plots.Measures
using LaTeXStrings
!isdir("plots") && mkdir("plots")
all_var_combs = [[0, 1]]
labels = [L"x", L"v"]

function clip_interval(X::Interval, var)
    if var == 1
        return Interval(max(clip_x[1], low(X, 1)), min(clip_x[2], high(X, 1)))
    elseif var == 2
        return Interval(max(clip_v[1], low(X, 1)), min(clip_v[2], high(X, 1)))
    end
    return x
end

scatter_plot = true
if isdefined(@__MODULE__, :reach_step_mountaincar)
    # compute outer bounds
    bounds = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    n = DTCSR.dim(X0)
    for R in reach_step_mountaincar
        l, h = get!(bounds, R.time, (fill(Inf, n), fill(-Inf, n)))
        for i in 1:n
            l[i] = min(l[i], low(R.Y, i))
            h[i] = max(h[i], high(R.Y, i))
        end
    end
    bounds = sort(bounds)

    cols = colormap("RdBu", length(bounds))
    cols = fill(:red, length(bounds))
    for vars in all_var_combs
        plot(leg=:topleft,
             xlab=(vars[1] == 0 ? "time" : labels[vars[1]]),
             ylab=(vars[2] == 0 ? "time" : labels[vars[2]]),
             yguidefontrotation=-90,
             legendfontsize=20,
             tickfont=font(20, "Times"),
             xguidefont=font(22, "Times"),
             yguidefont=font(22, "Times"),
             size=(1200, 800),
             bottom_margin=4mm,left_margin=8mm,top_margin=0mm,right_margin=0mm
            )

        for (i, tXi) in enumerate(bounds)
            t = tXi[1]
            Xi = Hyperrectangle(low=tXi[2][1], high=tXi[2][2])
            if 0 ∈ vars
                @assert vars[1] == 0
                Xi2 = Singleton([t]) × clip_interval(overapproximate(project(Xi, [vars[2]]), LazySets.Interval), vars[2])
            else
                Xi2 = project(Xi, vars)
            end
            plot!(Xi2, vars=vars, c=cols[i], alpha=0.8, ms=2)
        end

        for trajectory in sim_mountaincar
            if 0 ∈ vars
                @assert vars[1] == 0
                xs = [(i-1, x[vars[2]]) for (i, x) in enumerate(trajectory)]
            else
                xs = [(x[vars[1]], x[vars[2]]) for x in trajectory]
            end
            if scatter_plot
                scatter!(xs, c=:black, label="", alpha=1, ms=1)
            else
                plot!(xs, c=:black, alpha=0.4, label="")
            end
        end

        if 1 ∈ vars
            p! = vars[1] == 1 ? vline! : hline!
            p!([0.5], c=:blue, s=:dash, lab="goal")
        end

        savefig("plots/mountaincar_reach")
    end
end
