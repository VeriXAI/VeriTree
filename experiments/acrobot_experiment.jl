using DTCSR
using DTCSR: ApproximatedSet

# fix random seed
import Random
Random.seed!(1234)

LazySets.set_ztol(Float64, 1e-10)  # to avoid plotting artifacts

# --- problem definition ---

g = 9.8  # gravity
l1 = 1.0  # [m]
m1 = 1.0  # [kg] mass of link 1
m2 = 1.0  # [kg] mass of link 2
lc1 = 0.5  # [m] position of the center of mass of link 1
lc2 = 0.5  # [m] position of the center of mass of link 2
I1 = I2 = 1.0  # moments of inertia for both links
tau = 0.2  # control cycle length
clip_θ₁ = [-4π, 4π]
clip_θ₂ = [-9π, 9π]

@taylorize function acrobot!(ds, s, p, t)
    θ₁ = s[1]  # angle 1
    θ₂ = s[2]  # angle 2
    θ₁_dot = s[3]  # angular velocity 1
    θ₂_dot = s[4]  # angular velocity 2
    action = s[5]  # action

    d1 = m1 * lc1^2 +
         m2 * (l1^2 + lc2^2 + 2 * l1 * lc2 * cos(θ₂)) +
         I1 +
         I2
    d2 = m2 * (lc2^2 + l1 * lc2 * cos(θ₂)) + I2
    ϕ2 = m2 * lc2 * g * cos(θ₁ + θ₂ - pi / 2.0)
    ϕ1 = -m2 * l1 * lc2 * θ₂_dot^2 * sin(θ₂) -
           2 * m2 * l1 * lc2 * θ₂_dot * θ₁_dot * sin(θ₂) +
           (m1 * lc1 + m2 * l1) * g * cos(θ₁ - pi / 2.0) +
           ϕ2

    θ₂_dot_dot = (action + d2 / d1 * ϕ1 - m2 * l1 * lc2 * θ₁_dot^2 * sin(θ₂) - ϕ2) /
                 (m2 * lc2^2 + I2 - d2^2 / d1)
    θ₁_dot_dot = -(d2 * θ₂_dot_dot + ϕ1) / d1

    ds[1] = θ₁_dot
    ds[2] = θ₂_dot
    ds[3] = θ₁_dot_dot
    ds[4] = θ₂_dot_dot
    ds[5] = zero(action)
end

mapping = s -> [cos(s[1]), sin(s[1]), cos(s[2]), sin(s[2]), s[3], s[4], s[5]]
preprocess_fun = preprocessing_given_mapping(mapping)
mapping = s -> [acos(s[1]), acos(s[3]), s[5], s[6], s[7]]
preprocess_undo_fun = preprocessing_given_mapping(mapping)

left = Action("left", 0)
none = Action("none", 1)
right = Action("right", 2)

DT_depth_2 = Node(
    OrthogonalHalfSpace(5, 6, -0.097, true, false),
    Node(
        OrthogonalHalfSpace(6, 6, -1.131, true, false),
        Leaf(left),
        Leaf(right)
    ),
    Node(
        OrthogonalHalfSpace(6, 6, 1.582, true, false),
        Leaf(left),
        Leaf(right)
    )
)

DT_depth_3 = Node(
    OrthogonalHalfSpace(5, 6, -0.099, true, false),
    Node(
        OrthogonalHalfSpace(6, 6, -1.351, true, false),
        Node(
            OrthogonalHalfSpace(5, 6, -2.661, true, false),
            Leaf(right),
            Leaf(left)
        ),
        Leaf(right)
    ),
    Node(
        OrthogonalHalfSpace(6, 6, 1.123, true, false),
        Leaf(left),
        Leaf(right)
    )
)

Ps = ContSimProblem(acrobot!, tau)
Pr = TMContReachProblem(acrobot!, tau, TMJets21b(abstol=1e-5, orderT=6, orderQ=1))

X0 = UnionSetArray([
    Hyperrectangle([-0.05, 0, 0, 0, 0], [0.0, 0, 0, 0, 0]),
    Hyperrectangle([-0.07, 0, 0, 0, 0], [0.0005, 0.0005, 0, 0, 0])
])

max_steps = 500

# predicate for aborting
function abort_acrobot(tx::Tuple, i)
    if i > max_steps
        return time_bound
    end
    x = tx[2]
    if -cos(x[1]) - cos(x[2] + x[1]) > 1.0
        return other_condition
    else
        return no_abort
    end
end

function abort_acrobot(X::LazySet, i)
    if i > max_steps
        return time_bound
    else
        return no_abort
    end
end

function abort_acrobot(X::IntervalBox, i)
    if i > max_steps
        return time_bound
    elseif -cos(X[1]) - cos(X[2] + X[1]) > 1.0
        return other_condition
    else
        return no_abort
    end
end

function abort_acrobot(X::ApproximatedSet, i)
    return abort_acrobot(X.X, i)
end

function abort_acrobot(X, i)  # Taylor model
    if i > max_steps
        return time_bound
    end
    p = -cos(X[1]) - cos(X[2] + X[1])  # still a Taylor model
    if evaluate(p, p.dom) > 1.0  # evaluate to check the constraint
        return other_condition
    else
        return no_abort
    end
end

# predicate for aborting with crude clipping
function abort_acrobot_clipped(x::AbstractVector, i)
    res = abort_acrobot(x, i)
    if res != no_abort
        return res
    elseif x[1] ∉ Interval(clip_θ₁) || x[2] ∉ Interval(clip_θ₂)
        return other_condition
    else
        return no_abort
    end
end

verbose = false
benchmark = true

sim_acrobot_2 = nothing
reach_step_acrobot_2 = nothing
reach_acrobot_2 = nothing
sim_acrobot_3 = nothing
reach_step_acrobot_3 = nothing
reach_acrobot_3 = nothing
sim_acrobot_5 = nothing
reach_step_acrobot_5 = nothing
reach_acrobot_5 = nothing
sim_acrobot_10 = nothing
reach_step_acrobot_10 = nothing
reach_acrobot_10 = nothing

fixpoint_check = FixpointCheck(skip_redundancy_check=true, split=true)
for i in [1, 2]  # run experiment for two different decision trees
    if i == 1
        local DT = DT_depth_2
    elseif i == 2
        local DT = DT_depth_3
    end

    if benchmark
        println("running experiment with decision tree #$(i)")
    end

    # --- simulation ---

    sim_acrobot = simulate(X0, DT, Ps, abort_acrobot;
                           preprocess_fun=preprocess_fun, count=100);

    if i == 1
        global sim_acrobot_2 = sim_acrobot
    elseif i == 2
        global sim_acrobot_3 = sim_acrobot
    end

    # --- reachability ---

    if benchmark && i == 1
        # warm-up run for precompilation
        println("# Warm-up run:")
        reach(X0, DT, Pr, abort_acrobot; preprocess_fun=preprocess_fun,
              preprocess_undo_fun=preprocess_undo_fun,
              fixpoint_check=deepcopy(fixpoint_check), verbose=false);
    end
    println("# Verification run:")
    @time reach_step_acrobot, reach_acrobot =
        reach(X0, DT, Pr, abort_acrobot; preprocess_fun=preprocess_fun,
              preprocess_undo_fun=preprocess_undo_fun,
              fixpoint_check=deepcopy(fixpoint_check), verbose=verbose);

    if i == 1
        global reach_step_acrobot_2 = reach_step_acrobot
        global reach_acrobot_2 = reach_acrobot
    elseif i == 2
        global reach_step_acrobot_3 = reach_step_acrobot
        global reach_acrobot_3 = reach_acrobot
    end
end;

# --- plotting ---

import Plots
using Plots: plot, plot!, hline!, vline!, colormap, xlims!, ylims!, savefig, font
using LaTeXStrings
using Plots.Measures
!isdir("plots") && mkdir("plots")
labels = [L"\theta_1", L"\theta_2", L"\dot{\theta}_1", L"\dot{\theta}_2"]
x = (-3π/2:π/256:-π/2)
vars1 = [1, 2]

if isdefined(@__MODULE__, :reach_acrobot_2) && isdefined(@__MODULE__, :reach_acrobot_3) &&
        isdefined(@__MODULE__, :sim_acrobot_2) && isdefined(@__MODULE__, :sim_acrobot_3)
    ntdiv=3
    for i in [1, 2]#, 3, 4]
        if i == 1
            sim_acrobot = sim_acrobot_2
            reach_step_acrobot = reach_step_acrobot_2
            reach_acrobot = reach_acrobot_2
        elseif i == 2
            sim_acrobot = sim_acrobot_3
            reach_step_acrobot = reach_step_acrobot_3
            reach_acrobot = reach_acrobot_3
        end

        plot(leg=:bottomleft,
             xlab=(vars1[1] == 0 ? "t" : labels[vars1[1]]),
             ylab=(vars1[2] == 0 ? "t" : labels[vars1[2]]),
             yguidefontrotation=-90,
             legendfontsize=20,
             tickfont=font(20, "Times"),
             xguidefont=font(22, "Times"),
             yguidefont=font(22, "Times"),
             size=(1200, 800),
             bottom_margin=4mm,left_margin=9mm,top_margin=0mm,right_margin=0mm
            )

        cols = colormap("RdBu", length(reach_acrobot))
        cols = fill(:red, length(reach_acrobot))
        for (i, Fi) in enumerate(reach_acrobot)
            plot!(Fi, vars=vars1, c=cols[i], alpha=0.8)
        end

        xlims!(-3π, π)
        ylims!(-π, 5π)

        if vars1 == [1, 2]
            lab = "goal"
            for i in -1:1:1
                y1 = (x -> acos(-cos(x) - 1) - x + 2π * i).(x)
                y2 = (x -> -acos(-cos(x) - 1) - x + 2π * (i + 1)).(x)
                for j in -1:1:1
                    x2 = [xi + 2j*π for xi in x]
                    shapex = vcat(reverse(x2), x2)
                    shapey = vcat(reverse(y1), y2)
                    plot!(shapex, shapey, fill=true, label=lab, c=:blue, linewidth=0, fillalpha=0.5)
                    lab = ""
                end
                lab = ""
            end
        end

        for trajectory in sim_acrobot
            if 0 ∈ vars1
                @assert vars1[1] == 0
                xs = [(tx[1], tx[2][vars1[2]]) for tx in trajectory]
            else
                xs = [(tx[2][vars1[1]], tx[2][vars1[2]]) for tx in trajectory]
            end
            plot!(xs, color=:black, alpha=0.4, label="")
        end

        if (0 ∉ vars1)
            # the initial set is not visible at the normal scale, so we
            # exaggerate it
            X0p = project(X0, vars1)
            if X0p isa UnionSetArray
                X0q = UnionSetArray([Bloating(X, 0.1, Inf) for X in array(X0p)])
            else
                X0q = Bloating(X0p, 0.1, Inf)
            end
            plot!(X0q, alpha=1, lc=:yellow, c=:yellow, la=0.8, lab="initial")
        end

        savefig("plots/acrobot_reach_$i")
    end
end
