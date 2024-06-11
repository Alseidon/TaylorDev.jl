using .MTDevCore

using StaticArrays: Size

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


function epsilons(t::Type, orders::Tuple)
    res = MTDev(t, orders)
    ndimensions = ndims(res)
    for i in 1:ndimensions
        res.dev[(i == d ? 2 : 1 for d in 1:ndimensions)...] = 1
    end
    return res
end

epsilons(orders::Tuple) = epsilons(Float64, orders)
epsilons(ndimensions::Int, orders::Int=1) = epsilons(
    Tuple(orders for _ in 1:ndimensions))
epsilons(t::Type, ndimensions::Int, orders::Int=1) = epsilons(t,
    Tuple(orders for _ in 1:ndimensions))

function _lowest_nonzero_order_i(mt::MTDev, i::Int)
    for ord in 0:orders(mt)[i]
        if any(mt.dev[[k == i ? ord+1 : (:) for k  in 1:ndims(mt)]...] .!= 0)
            return ord
        end
    end
    return order(mt) + 1
end

lowest_nonzero_orders(mt::MTDev) = map(
    i->_lowest_nonzero_order_i(mt, i),
    1:ndims(mt))

constant_term(a::MTDev) = a.dev[1]


function to_mtdev(nb::Number, orders::Tuple)
    res = MTDev(typeof(nb), orders)
    res.dev[1] = nb
    return res
end

function compute(mtdev::MTDev, eps::Vector)
    res = zero(eltype(mtdev))
    for i in get_coeffs(mtdev)
        res += mtdev.dev[i...] * prod(eps .^ (i .- 1))
    end
    return res
end

function deriv(mtdev::MTDev, nvar::Int)
    ord = orders(mtdev) |> collect
    ord[nvar] -= 1
    @assert ord[nvar] >= 0
    res = MTDev(eltype(mtdev), Tuple(ord))
    for i in get_coeffs(res)
        j = collect(i)
        j[nvar] += 1
        res.dev[i...] = mtdev.dev[j...] * i[nvar]
    end
    return res
end

# BASE OPERATIONS
Base.:(==)(a::MTDev, b::MTDev) = (a.dev == b.dev)

function Base.:+(a::MTDev{Ns}, b::MTDev{Ms}) where {Ns, Ms}
    idx = map((i, j)->1:min(i,j), Tuple(Size(Ns)), Tuple(Size(Ms)))
    return MTDev(a.dev[idx...] + b.dev[idx...])
end

function Base.:-(a::MTDev{Ns}, b::MTDev{Ms}) where {Ns, Ms}
    idx = map((i, j)->1:min(i,j), Tuple(Size(Ns)), Tuple(Size(Ms)))
    return MTDev(a.dev[idx...] - b.dev[idx...])
end

function Base.:-(a::MTDev)
    return MTDev(-a.dev)
end

function Base.:*(a::MTDev, b::MTDev)
    # TODO different orders
    ord = orders(a)
    @assert ord == orders(b)
    res = MTDev(eltype(a), ord)
    for i in get_coeffs(res)
        for j in collect(Iterators.product(map(c->1:c, i)...))
            res.dev[i...] += a.dev[j...] * b.dev[(i .- j .+ 1)...]
        end
    end
    return res
end


function Base.:/(a::MTDev, b::MTDev)
    # TODO
    @assert b.dev[1] != 0
    @error "Not implemented"
end


# INTERACTION WITH NUMBERS
Base.:(==)(nb::Number, b::MTDev) = false
Base.:(==)(b::MTDev, nb::Number) = false

Base.isless(a::MTDev, nb::Number) = constant_term(a) < nb
Base.isless(nb::Number, a::MTDev) = nb < constant_term(a)


function Base.:+(nb::Number, b::MTDev)
    return to_mtdev(nb, orders(b)) + b
end
Base.:+(b::MTDev, nb::Number) = nb + b

function Base.:-(nb::Number, b::MTDev)
    return to_mtdev(nb, orders(b)) - b
end
function Base.:-(b::MTDev, nb::Number)
    return b - to_mtdev(nb, orders(b))
end

function Base.:*(nb::Number, b::MTDev)
    return MTDev(nb * b.dev)
end

Base.:*(b::MTDev, nb::Number) = nb * b

Base.:/(a::MTDev, nb::Number) = MTDev(a.dev ./ nb)

Base.:/(nb::Number, a::TDev) = to_mtdev(nb, orders(a)) / a

#=
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

# INTEG UTILS
function identity_jet(v::Vector, ord::Int=1)
    ndimensions = length(v)
    res = zeros(MTDev, ndimensions)
    ords = Tuple(ord for i in 1:ndimensions)
    for i in 1:ndimensions
        res[i] = v[i] + MTDev(eltype(v), ords)
        j = [dim==i ? 2 : 1 for dim in 1:ndimensions]
        res[i].dev[j...] = 1
    end
    return res
end