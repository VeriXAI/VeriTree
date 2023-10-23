using ReachabilityAnalysis

X0 = BallInf([0.4, 0.4], 0.2)
A1 = [1.739401  0.241926;
      0.049354  0.971989]
T = 1.0
P = @ivp(x' = A1 * x, x(0) ∈ X0)
sol = solve(P, T=T);

X1 = linear_map(exp(A1*T), X0)
H = HalfSpace([1.0, 0], 3.0)
HC = HalfSpace([-1.0, 0], -3.0)

Z1 = intersection(X1, H)
Z2 = intersection(X1, HC)

BZ1 = box_approximation(Z1)
BZ2 = box_approximation(Z2)

using Plots, Plots.Measures, LaTeXStrings

plot(X0, c=:yellow, alpha=0.7, lab=L"X_0")
plot!(sol[1], vars=[1, 2], c=:cornflowerblue, alpha=0.8, lab=L"Y")
plot!(sol, vars=[1, 2], c=:cornflowerblue, alpha=0.1, lab="")
plot!(X0, c=:yellow, alpha=0.7, lab="")
plot!(X1, c=:purple, alpha=0.7, lab=L"Z")
plot!(X1, alpha=0, linealpha=1, lc=:black, lw=3, ls=:dash, lab="")
plot!(leg=:topleft,
      xlab=L"x", ylab=L"y",
      yguidefontrotation=-90,
      legendfontsize=13,
      tickfont=font(13, "Times"),
      xguidefont=font(15, "Times"),
      yguidefont=font(15, "Times"),
      xlims=(0.2, 4.05), ylims=(0.2, 1.75),
      bottom_margin=0mm,left_margin=4mm)
savefig("plots/algorithm_example_1")

plot(BZ1, c=:cyan, lab="box("*L"\,Z_1\,"*")")
plot!(BZ2, c=:orange, lab="box("*L"\,Z_2\,"*")")
plot!(Z1, c=:green, alpha=0.8, lab=L"Z_1")
plot!(Z2, c=:brown, alpha=0.8, lab=L"Z_2")
plot!(X1, alpha=0, linealpha=1, lc=:black, lw=3, ls=:dash)
plot!([3.0, 3.0], [0.2, 1.75], ls=:dash, lc=:red, lw=3, lab="")
plot!(leg=:topleft,
      xlab=L"x", ylab=L"y",
      yguidefontrotation=-90,
      legendfontsize=13,
      tickfont=font(13, "Times"),
      xguidefont=font(15, "Times"),
      yguidefont=font(15, "Times"),
      xlims=(0.2, 4.05), ylims=(0.2, 1.75),
      bottom_margin=0mm,left_margin=4mm)
savefig("plots/algorithm_example_2")
