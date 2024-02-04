module TaylorDev
export my_const
export TDev,  order,  to_tdev,  epsilon              # TDev
export MTDev, orders, to_mtdev, epsilons, get_coeffs # MTDev 
export compute, deriv                                # Commons
export identity_jet, taylor_method

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
