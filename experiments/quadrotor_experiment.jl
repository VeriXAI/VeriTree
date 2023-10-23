using DTCSR

# fix random seed
import Random
Random.seed!(1234)

# for condition on time
function Base.isless(x::Taylor1, y::Number)
    return isless(evaluate(x), y) && !(evaluate(x) ≈ y)
end

# --- problem definition ---

g = 9.81  # gravity
tau = 0.2  # control cycle length

# source: eqs. (23)-(24) in https://coneural.org/florian/papers/05_cart_pole.pdf
@taylorize function quadrotor!(dx, x, p, t)
    if t < 10 * tau
        bx = 0.25
        by = -0.25
        bz = 0.0
    elseif t < 20 * tau
        bx = 0.25
        by = 0.25
        bz = 0.0
    elseif t < 25 * tau
        bx = 0.0
        by = 0.25
        bz = 0.0
    else
        bx = -0.25
        by = 0.25
        bz = 0.0
    end

    # px = x[1]
    # py = x[2]
    # pz = x[3]
    vx = x[4]
    vy = x[5]
    vz = x[6]
    action = x[7]

    # action conversion
    if action == 1
        θ = 0.1
        ϕ = 0.1
        τ = 11.81
    elseif action == 2
        θ = 0.1
        ϕ = 0.1
        τ = 7.81
    elseif action == 3
        θ = 0.1
        ϕ = -0.1
        τ = 11.81
    elseif action == 4
        θ = 0.1
        ϕ = -0.1
        τ = 7.81
    elseif action == 5
        θ = -0.1
        ϕ = 0.1
        τ = 11.81
    elseif action == 6
        θ = -0.1
        ϕ = 0.1
        τ = 7.81
    elseif action == 7
        θ = -0.1
        ϕ = -0.1
        τ = 11.81
    elseif action == 8
        θ = -0.1
        ϕ = -0.1
        τ = 7.81
    end

    dx[1] = vx - bx
    dx[2] = vy - by
    dx[3] = vz - bz
    dx[4] = one(vx) * (g * tan(θ))
    dx[5] = one(vy) * (-g * tan(ϕ))
    dx[6] = one(vz) * (τ - g)
    dx[7] = zero(action)
    return dx
end

action1 = Action("1", 1)
action2 = Action("2", 2)
action3 = Action("3", 3)
action4 = Action("4", 4)
action5 = Action("5", 5)
action6 = Action("6", 6)
action7 = Action("7", 7)
action8 = Action("8", 8)

include("quadrotor_tree.jl")
DT = DT_depth_10

Ps = ContSimProblem(quadrotor!, tau)
Pr = TMContReachProblem(quadrotor!, tau,
                        TMJets21b(abstol=1e-8, orderT=10, orderQ=1, adaptive=false))

# the first two states are uncertain in [-0.05, 0.05]
X0 = Hyperrectangle(zeros(7), [0.05, 0.05, 0, 0, 0, 0, 0])

kmax = 30
Tmax = kmax * tau
function abort_quadrotor(X, i)  # predicate for aborting
    if i > kmax
        return time_bound
    else
        return no_abort
    end
end

verbose = false
benchmark = true

# --- simulation ---

sim_quadrotor = simulate(X0, DT, Ps, abort_quadrotor, count=1,
                         include_vertices=true);

# --- reachability ---

function reference(t::Real)
    x = [0.0, 0.0]
    for (Δt, Δx) in [(10 * tau, [0.25, -0.25]),
                     (10 * tau, [0.25, 0.25]),
                     (5 * tau, [0.0, 0.25]),
                     (5 * tau, [-0.25, 0.25])]
        t_next = min(t, Δt)
        x .+= t_next .* Δx
        t -= t_next
        if t <= 0
            break
        end
    end
    return x
end

goal_area = intersection(overapproximate(Ball2([0.81, 0.5], 0.26), 1e-3), HalfSpace([0, 1.0], 0.5))

# check inclusion in goal area
function specification(reach_quadrotor)
    validated = true
    offset = reference(Tmax)
    for Fi in reach_quadrotor
        if Tmax ∉ tspan(Fi)
            continue
        end
        final_set = translate(set(project(overapproximate(Fi(Tmax), Zonotope), [1, 2])), offset)
        validated &= final_set ⊆ goal_area
    end
    return validated
end

fixpoint_check = NoFixpointCheck()
if benchmark
    # warm-up run for precompilation
    println("# Warm-up run:")
    reach_step_quadrotor, reach_quadrotor = reach(X0, DT, Pr, abort_quadrotor,
            fixpoint_check=deepcopy(fixpoint_check), verbose=verbose);
    validated = specification(reach_quadrotor)
    @assert validated "the specification is violated"
end
println("# Verification run:")
@time begin
    reach_step_quadrotor, reach_quadrotor = reach(X0, DT, Pr, abort_quadrotor,
            fixpoint_check=deepcopy(fixpoint_check), verbose=verbose);
    if specification(reach_quadrotor)
        println("the specification was verified")
    else
        @assert false "the specification is violated"
    end
end

# --- plotting ---

import Plots
using Plots: plot, plot!, scatter!, hline!, vline!, colormap, savefig, font
using Plots.Measures
using LaTeXStrings
!isdir("plots") && mkdir("plots")
labels = [L"p_x", L"p_y", L"p_z", L"v_x", L"v_y", L"v_z"]

@taylorize function quadrotor_ref!(dx, x, p, t)
    if t < 10 * tau
        bx = 0.25
        by = -0.25
        bz = 0.0
    elseif t < 20 * tau
        bx = 0.25
        by = 0.25
        bz = 0.0
    elseif t < 25 * tau
        bx = 0.0
        by = 0.25
        bz = 0.0
    else
        bx = -0.25
        by = 0.25
        bz = 0.0
    end

    dx[1] = bx
    dx[2] = by
    dx[3] = bz
    dx[4] = zero(x[4])
    return dx
end

Ps_ref = ContSimProblem(quadrotor_ref!, tau)
X0_ref = Singleton(zeros(4))
DT_ref = Node(
    OrthogonalHalfSpace(1, 4, 0.0, true, false),
    Leaf(action1),
    Leaf(action1)
)
sim_quadrotor_ref = simulate(X0_ref, DT_ref, Ps_ref, abort_quadrotor,
                             count=1, include_vertices=false);

if isdefined(@__MODULE__, :sim_quadrotor_ref) && isdefined(@__MODULE__, :sim_quadrotor)
    vars = [1, 2]
    plot(leg=:bottomright,
         xlab=labels[1], ylab=labels[2],
         yguidefontrotation=-90,
         legendfontsize=20,
         tickfont=font(20, "Times"),
         xguidefont=font(22, "Times"),
         yguidefont=font(22, "Times"),
         size=(1200, 800),
         bottom_margin=4mm,left_margin=9mm,top_margin=0mm,right_margin=0mm
    )

    plot!(project(X0, vars), alpha=0.8, c=:yellow, lc=:yellow, lab="initial")
    plot!(goal_area, alpha=0.5, c=:cornflowerblue, lc=:cornflowerblue, lab="goal")

    traj_ref = sim_quadrotor_ref[1]
    traj_ref_shortened = traj_ref[1:end-15]
    plot!([(tx[2][1], tx[2][2]) for tx in traj_ref_shortened],
          c=:brown, lw=2, label="reference")
    @assert length(traj_ref) == length(sim_quadrotor[1])
    for traj in sim_quadrotor
        plot!([(traj[i][2][1] .+ traj_ref[i][2][1],
                traj[i][2][2] .+ traj_ref[i][2][2]) for i in eachindex(traj)],
              label="")
    end

    x0 = (0.0, 0.0)
    t1 = 10 * tau
    x1 = x0 .+ t1 .* (0.25, -0.25)
    t2 = 10 * tau
    x2 = x1 .+ t2 .* (0.25, 0.25)
    t3 = 5 * tau
    x3 = x2 .+ t3 .* (0.0, 0.25)
    t4 = 5 * tau
    x4 = x3 .+ t4 .* (-0.25, 0.25)

    for o in [-0.32, 0.32]
        plot!([x0 .+ (0, o), x1 .+ (0, o), x2 .+ (-o, 0.0),
               x3 .+ (-o, 0), x4 .+ (-o, 0.0)],
               alpha=0, la=1, lc=:blue, lw=2, ls=:dash,
               label=(o == -0.32 ? "boundary" : ""))
    end

    savefig("plots/quadrotor_simulation_absolute")
end

if isdefined(@__MODULE__, :sim_quadrotor_ref) && isdefined(@__MODULE__, :sim_quadrotor) &&
         isdefined(@__MODULE__, :reach_step_quadrotor)
    vars = [1, 2]
    plot(leg=:bottomright,
         xlab=labels[1], ylab=labels[2],
         yguidefontrotation=-90,
         legendfontsize=20,
         tickfont=font(20, "Times"),
         xguidefont=font(22, "Times"),
         yguidefont=font(22, "Times"),
         size=(1200, 800),
         bottom_margin=4mm,left_margin=9mm,top_margin=0mm,right_margin=0mm
    )

    ntdiv = 5
    local cols = fill(:red, length(reach_quadrotor))
    for (i, Fi) in enumerate(reach_quadrotor)
        for Xj in Fi
            for B in overapproximate(Xj, Hyperrectangle, ntdiv=ntdiv)
                t = mid(B.Δt)
                offset = reference(t)
                BB = translate(project(set(B), vars), offset)
                plot!(BB, vars=vars, c=cols[i], alpha=0.2)
            end
        end
    end

    plot!(project(X0, vars), alpha=0.8, c=:yellow, lc=:yellow, lab="initial")
    plot!(goal_area, alpha=0.5, c=:cornflowerblue, lc=:cornflowerblue, lab="goal")

    traj_ref = sim_quadrotor_ref[1]
    traj_ref_shortened = traj_ref[1:end-15]
    plot!([(tx[2][1], tx[2][2]) for tx in traj_ref_shortened],
          c=:brown, lw=2, label="reference")
    @assert length(traj_ref) == length(sim_quadrotor[1])
    for traj in sim_quadrotor
        plot!([(traj[i][2][1] .+ traj_ref[i][2][1],
                traj[i][2][2] .+ traj_ref[i][2][2]) for i in eachindex(traj)],
              label="")
    end

    x0 = (0.0, 0.0)
    t1 = 10 * tau
    x1 = x0 .+ t1 .* (0.25, -0.25)
    t2 = 10 * tau
    x2 = x1 .+ t2 .* (0.25, 0.25)
    t3 = 5 * tau
    x3 = x2 .+ t3 .* (0.0, 0.25)
    t4 = 5 * tau
    x4 = x3 .+ t4 .* (-0.25, 0.25)

    for o in [-0.32, 0.32]
        plot!([x0 .+ (0, o), x1 .+ (0, o), x2 .+ (-o, 0.0),
               x3 .+ (-o, 0), x4 .+ (-o, 0.0)],
               alpha=0, la=1, lc=:blue, lw=2, ls=:dash,
               label=(o == -0.32 ? "boundary" : ""))
    end

    savefig("plots/quadrotor_reach_absolute")
end
