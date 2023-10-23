using DTCSR

# fix random seed
import Random
Random.seed!(1234)

# --- problem definition ---

gravity = 9.8
masscart = 1.0
masspole = 0.1
total_mass = masspole + masscart
pole_length = 0.5  # actually half the pole's length
polemass_length = masspole * pole_length
force_mag = 10.0
tau = 0.02  # control cycle length

@taylorize function cartpole!(ds, s, p, t)
#     x = s[1]  # cart position
    x_dot = s[2]  # cart velocity
    θ = s[3]  # pendulum angle
    θ_dot = s[4]  # pendulum angular velocity
    action = s[5]  # action
    if action == 1
        force = force_mag
    else
        force = -force_mag
    end
    cosθ = cos(θ)
    sinθ = sin(θ)

    temp = (force + polemass_length * θ_dot^2 * sinθ) / total_mass
    θ_acc = (gravity * sinθ - cosθ * temp) / (pole_length * (4.0/3.0 - masspole * cosθ^2 / total_mass))
    x_acc = temp - polemass_length * θ_acc * cosθ / total_mass

    ds[1] = x_dot
    ds[2] = x_acc
    ds[3] = θ_dot
    ds[4] = θ_acc
    ds[5] = zero(action)
end

left = Action("left", 0)
right = Action("right", 1)

DT_depth_2 = Node(
    OrthogonalHalfSpace(3, 4, -0.014, true, false),
    Node(
        OrthogonalHalfSpace(1, 4, -0.221, true, false),
        Leaf(left),
        Leaf(left)
    ),
    Node(
        OrthogonalHalfSpace(4, 4, -0.201, true, false),
        Leaf(left),
        Leaf(right)
    )
)

DT = DT_depth_2

Ps = ContSimProblem(cartpole!, tau)
Pr = TMContReachProblem(cartpole!, tau, TMJets21b(abstol=1e-9, orderT=6, orderQ=1))

X0 = UnionSetArray([Hyperrectangle([0.0, 0, 0, 0, 0], [0.05, 0.05, 0.05, 0.05, 0]),
                    Hyperrectangle([0.0, 0, 0, 0, 0], [0, 0, 0.025, 0.2, 0])])

function abort_cartpole(X, i)  # predicate for aborting
    if i > 25
        return time_bound
    else
        return no_abort
    end
end

function abort_cartpole_longer(X, i)  # predicate for aborting (long run)
    if i > 100
        return time_bound
    else
        return no_abort
    end
end

function learning_predicate_cartpole_all(X, i)
    return true
end

verbose = false
benchmark = true

# --- simulation ---

sim_cartpole = simulate(X0, DT, Ps, abort_cartpole, count=100);
# simulations for separate plot with longer time horizon
sim_cartpole_motivation = simulate(X0, DT, Ps, abort_cartpole_longer, count=100);

# --- reachability ---

fixpoint_check = FixpointCheck(dimensions=3:4, tracker=BoxTracker(),
                               learning_predicate=learning_predicate_cartpole_all,
                               split=true)
if benchmark
    # warm-up run for precompilation
    println("# Warm-up run:")
    reach(X0, DT, Pr, abort_cartpole, fixpoint_check=deepcopy(fixpoint_check), verbose=false);
end
println("# Verification run:")
@time reach_step_cartpole, reach_cartpole =
    reach(X0, DT, Pr, abort_cartpole, fixpoint_check=deepcopy(fixpoint_check), verbose=verbose);

# --- plotting ---

import Plots
using Plots: plot, plot!, scatter!, hline!, vline!, colormap, savefig, font
using Plots.Measures
using LaTeXStrings
!isdir("plots") && mkdir("plots")
labels = [L"p", L"v", L"\theta", L"\omega"]

if isdefined(@__MODULE__, :sim_cartpole_motivation)
    all_var_combs = [[0, 3]]
    for vars in all_var_combs
        plot(leg=:outerright,
             xlab=(vars[1] == 0 ? "time" : labels[vars[1]]),
             ylab=(vars[2] == 0 ? "time" : labels[vars[2]]),
             yguidefontrotation=-90,
             legendfontsize=20,
             tickfont=font(20, "Times"),
             xguidefont=font(22, "Times"),
             yguidefont=font(22, "Times"),
             yticks=[-0.06, 0, 0.06],
             ylims=[-0.06, 0.06],
             size=(1200, 800),
             bottom_margin=4mm,left_margin=8mm,top_margin=1mm,right_margin=0mm
        )

        for trajectory in sim_cartpole_motivation
            if 0 ∈ vars
                @assert vars[1] == 0
                xs = [(tx[1], tx[2][vars[2]]) for tx in trajectory]
            else
                xs = [(tx[2][vars[1]], tx[2][vars[2]]) for tx in trajectory]
            end
            plot!(xs, label="")
        end

        savefig("plots/cartpole_simulations")
    end
end

if isdefined(@__MODULE__, :reach_step_cartpole) && isdefined(@__MODULE__, :sim_cartpole)
    all_var_combs = [[3, 4]]
    ntdiv = 3
    col_act = [:green, :cyan]
    for vars in all_var_combs
        plot(leg=:topright,
             xlab=(vars[1] == 0 ? "time" : labels[vars[1]]),
             ylab=(vars[2] == 0 ? "time" : labels[vars[2]]),
             yguidefontrotation=-90,
             legendfontsize=20,
             tickfont=font(20, "Times"),
             xguidefont=font(22, "Times"),
             yguidefont=font(22, "Times"),
             xticks=[-0.03, 0, 0.03, 0.06],
             size=(1200, 800),
             bottom_margin=4mm,left_margin=8mm,top_margin=0mm,right_margin=0mm
            )
        if (0 ∉ vars)
            Xnormal = box_approximation(ConvexHullArray([X.Y for X in reach_step_cartpole]))
            bloat = Hyperrectangle(zeros(dim(Xnormal)), radius_hyperrectangle(Xnormal) ./ 10)
            Xnormal = minkowski_sum(Xnormal, bloat)
            plot!(ActionSpace(DT, Xnormal; ST=Hyperrectangle), color_palette=col_act, alpha=0.6, vars=vars)
        end

        cols = colormap("RdBu", length(reach_cartpole))
        cols = fill(:red, length(reach_cartpole))
        for (i, Fi) in enumerate(reach_cartpole)
            for Xj in Fi
                plot!(overapproximate(Xj, Hyperrectangle, ntdiv=ntdiv), vars=vars, c=cols[i], alpha=0.2)
            end
        end

        for trajectory in sim_cartpole
            if 0 ∈ vars
                @assert vars[1] == 0
                xs = [(tx[1], tx[2][vars[2]]) for tx in trajectory]
            else
                xs = [(tx[2][vars[1]], tx[2][vars[2]]) for tx in trajectory]
            end
            plot!(xs, c=:black, alpha=0.4, label="")
        end

        if (0 ∉ vars)
            plot!(project(X0, vars), alpha=0, lw=5, lc=:yellow, la=0.8, lab="initial")
        end

        if 3 ∈ vars
            p! = vars[1] == 3 ? vline! : hline!
            p!([0.06], c=:blue, s=:dash, lw=3, lab="avoid")
        end

        savefig("plots/cartpole_reach")
    end
end
