module DTCSR

using Reexport, RecipesBase
@reexport using ReachabilityAnalysis
using ReachabilityAnalysis: zeroBox, symBox, zeroI, TaylorModel1, TaylorModelN, fp_rpa
import ReachabilityAnalysis.LazySets:
    dim, complement, intersection, isdisjoint, ⊆, low, high, isempty, translate,
    minkowski_sum, linear_map, project, ∈, HalfSpace, Hyperrectangle
using ReachabilityAnalysis.LazySets: @commutative
import Base: convert, iterate

export OpenHyperrectangle, OrthogonalHalfSpace,
       preprocessing_given_mapping,
       Action, Leaf, UnidimensionalBinaryDecisionTree, Node,
       ActionSpace,
       ContSimProblem, DiscSimProblem, simulate,
       TMContReachProblem, TMDiscReachProblem, reach,
       fixpoint_check,
       IntervalArith,
       no_abort, time_bound, other_condition,
       BoxTracker, FixpointCheck, NoFixpointCheck

include("util/IntervalArith.jl")
include("util/Abort.jl")

include("sets/OrthogonalHalfSpace.jl")
include("sets/OpenHyperrectangle.jl")
include("sets/ApproximatedSet.jl")
include("sets/overapproximate.jl")
include("sets/isdisjoint.jl")
include("sets/issubset.jl")
include("sets/intersection.jl")
include("sets/split.jl")
include("sets/minkowski_sum.jl")
include("sets/project.jl")
include("sets/read_action.jl")
include("sets/preprocessing.jl")
include("sets/merge_dimensions.jl")
include("sets/partition.jl")

include("tree/Action.jl")
include("tree/AbstractDecisionTree.jl")
include("tree/Leaf.jl")
include("tree/UnidimensionalBinaryDecisionTree.jl")
include("tree/AbstractUnidimensionalBinaryDecisionTree.jl")
include("tree/ActionSpace.jl")

include("simulation/AbstractSimProblem.jl")
include("simulation/ContSimProblem.jl")
include("simulation/DiscSimProblem.jl")
include("simulation/simulate.jl")

include("reachability/AbstractReachProblem.jl")
include("reachability/TMContReachProblem.jl")
include("reachability/TMDiscReachProblem.jl")
include("reachability/BoxTracker.jl")
include("reachability/NoFixpointCheck.jl")
include("reachability/FixpointCheck.jl")
include("reachability/continuous.jl")
include("reachability/discrete.jl")
include("reachability/reach.jl")

end  # module
