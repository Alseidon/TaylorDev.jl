module TaylorDev
export TDev, my_const, order, to_tdev, compute, deriv, epsilon
export MTDev, epsilons, get_coeffs  
export taylor_method

include("core.jl")
include("core_multi.jl")
include("functions_single.jl")
include("functions_multi.jl")
include("printing.jl")
include("integration.jl")

using .TDevCore
using .MTDevCore

my_const = 4# randn()

end
