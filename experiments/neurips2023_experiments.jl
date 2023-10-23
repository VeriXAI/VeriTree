println("Plots are created in the folder `plots`.")
mkpath("plots")
println("Verification is always run twice to trigger precompilation.\n")

# create Figure 3
module Figure
    println("### Creating Figure 3")
    include("algorithm_example.jl")
    println("-----\n")
end

# run cart/pole example; create Figure 5(a) and Figure 8
module CartPole
    println("### Cart/Pole Experiment")
    include("cartpole_experiment.jl")
    println("-----\n")
end

# run acrobot example; create Figure 6
module Acrobot
    println("### Acrobot Experiment")
    include("acrobot_experiment.jl")
    println("-----\n")
end

# run mountain/car example; create Figure 5(b)
module MountainCar
    println("### Mountain/Car Experiment")
    include("mountaincar_experiment.jl")
    println("-----\n")
end

# run quadrotor example; create Figure 1(a) and Figure 4
module Quadrotor
    println("### Quadrotor Experiment")
    include("quadrotor_experiment.jl")
    println("-----\n")
end

println("all done!")
