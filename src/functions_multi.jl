using .MTDevCore


orders(mtdev::MTDev) = size(mtdev.dev) .- 1
Base.ndims(mtdev::MTDev) = ndims(mtdev.dev)
Base.copy(mtdev::MTDev) = MTDev(mtdev.dev)
Base.eltype(mtdev::MTDev) = eltype(mtdev.dev)
Base.zero(mtdev::MTDev) = MTDev(eltype(mtdev), orders(mtdev))
Base.zero(::Type{MTDev}) = MTDev([0.])
function Base.one(mtdev::MTDev)
    res = zero(mtdev)
    res.dev[1] = one(eltype(mtdev))
    return res
end
Base.one(::Type{MTDev}) = MTDev([1.])


function epsilon(t::Type, orders::Tuple)
    res = MTDev(t, orders)
    ndimensions = ndims(res)
    for i in 1:ndimensions
        res.dev[(i == d ? 2 : 1 for d in 1:ndimensions)...] = 1
    end
    return res
end

epsilon(orders::Tuple) = epsilon(Float64, orders)
epsilon(ndimensions::Int, orders::Int=1) = epsilon(
    Tuple(orders for _ in 1:ndimensions))
epsilon(t::Type, ndimensions::Int, orders::Int=1) = epsilon(t,
    Tuple(orders for _ in 1:ndimensions))


function to_tdev(nb::Number, orders::Tuple)
    res = MTDev(typeof(nb), orders)
    res.dev[1] = nb
    return res
end

#=

function compute(mtdev::MTDev, eps::Vector)
    res = zero(eltype(mtdev))
    for i in get_coeffs(mtdev.dev)
        res += mtdev.dev[i] * eps^(i-1)
    end
    return res
end

function deriv(tdev::MTDev)
    @assert order(tdev) >= 1
    res = MTDev(eltype(tdev), order(tdev)-1)
    for i in eachindex(res.dev)
        res.dev[i] = tdev.dev[i+1] * i
    end
    return res
end

# BASE OPERATIONS
Base.:(==)(a::MTDev, b::MTDev) = (a.dev == b.dev)
function Base.:+(a::MTDev, b::MTDev)
    return MTDev(a.dev + b.dev)
end

function Base.:-(a::MTDev, b::MTDev)
    return MTDev(a.dev - b.dev)
end

function Base.:-(a::MTDev)
    return MTDev(-a.dev)
end

function Base.:*(a::MTDev, b::MTDev)
    res = MTDev(eltype(a), order(a))
    for i in eachindex(a.dev)
        for j in 1:(order(a)-i+2)
            res.dev[i+j-1] += a.dev[i] * b.dev[j]
        end
    end
    return res
end

function Base.:/(a::MTDev, b::MTDev)
    @assert b.dev[1] != 0
    res = copy(a)#TDev(eltype(a), order(a))
    for i in eachindex(res.dev) # k=i-1; c_k = a_k
        for j in 2:i
            res.dev[i] -= b.dev[j] * res.dev[i-j+1]
        end
        res.dev[i] /= b.dev[1]
    end
    return res
end

# INTERACTION WITH NUMBERS
Base.:(==)(nb::Number, b::MTDev) = false
Base.:(==)(b::MTDev, nb::Number) = false

Base.isless(a::MTDev, nb::Number) = a.dev[1] < nb
Base.isless(nb::Number, a::MTDev) = nb < a.dev[1]


function Base.:+(nb::Number, b::MTDev)
    return to_tdev(nb, order(b)) + b
end
Base.:+(b::MTDev, nb::Number) = nb + b

function Base.:-(nb::Number, b::MTDev)
    return to_tdev(nb, order(b)) - b
end
function Base.:-(b::MTDev, nb::Number)
    return b - to_tdev(nb, order(b))
end

function Base.:*(nb::Number, b::MTDev)
    return MTDev(nb * b.dev)
end

Base.:*(b::MTDev, nb::Number) = nb * b

Base.:/(a::MTDev, nb::Number) = MTDev(a.dev ./ nb)

function Base.:^(a::MTDev, nb::Integer)
    if nb == 1
        return a
    elseif nb == 0
        return MTDev(eltype(a), order(a))
    end
end

# BASE FUNCTIONS
function Base.sin(a::MTDev)
end


# IN-PLACE OPERATIONS
function add!(target::MTDev, b::MTDev)
    #@assert order(target) == order(b)
    target.dev .= target.dev + b.dev
    return
end
=#