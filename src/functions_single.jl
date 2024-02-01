using .TDevCore

order(tdev::TDev) = length(tdev.dev) - 1
Base.copy(tdev::TDev) = TDev(tdev.dev)
Base.eltype(tdev::TDev) = eltype(tdev.dev)
Base.zero(tdev::TDev) = TDev(eltype(tdev), order(tdev))
Base.zero(::Type{TDev}) = TDev([0.])
function Base.one(tdev::TDev)
    res = zero(tdev)
    res.dev[1] = one(eltype(res))
    return res
end
Base.one(::Type{TDev}) = TDev([1.])

epsilon(order::Int=1) = TDev([1*(i==2) for i in 1:(order+1)])
epsilon(t::Type, order::Int=1) = TDev([(i==2 ? one(t) : zero(t)) for i in 1:(order+1)])

function to_tdev(nb::Number, ord::Int)
    t = TDev([(i==1 ? nb : zero(nb)) for i in 1:(ord+1)])#TDev(typeof(nb), ord)
    #t.dev[1] = nb
    return t
end

function compute(tdev::TDev, eps::Number)
    res = zero(eltype(tdev))
    for i in eachindex(tdev.dev)
        res += tdev.dev[i] * eps^(i-1)
    end
    return res
end

function deriv(tdev::TDev)
    @assert order(tdev) >= 1
    res = TDev(eltype(tdev), order(tdev)-1)
    for i in eachindex(res.dev)
        res.dev[i] = tdev.dev[i+1] * i
    end
    return res
end

# BASE OPERATIONS
Base.:(==)(a::TDev, b::TDev) = (a.dev == b.dev)
function Base.:+(a::TDev, b::TDev)
    return TDev(a.dev + b.dev)
end

function Base.:-(a::TDev, b::TDev)
    return TDev(a.dev - b.dev)
end

function Base.:-(a::TDev)
    return TDev(-a.dev)
end

function Base.:*(a::TDev, b::TDev)
    res = TDev(eltype(a), order(a))
    for i in eachindex(a.dev)
        for j in 1:(order(a)-i+2)
            res.dev[i+j-1] += a.dev[i] * b.dev[j]
        end
    end
    return res
end

function Base.:/(a::TDev, b::TDev)
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
Base.:(==)(nb::Number, b::TDev) = false
Base.:(==)(b::TDev, nb::Number) = false

Base.isless(a::TDev, nb::Number) = a.dev[1] < nb
Base.isless(nb::Number, a::TDev) = nb < a.dev[1]


function Base.:+(nb::Number, b::TDev)
    return to_tdev(nb, order(b)) + b
end
Base.:+(b::TDev, nb::Number) = nb + b

function Base.:-(nb::Number, b::TDev)
    return to_tdev(nb, order(b)) - b
end
function Base.:-(b::TDev, nb::Number)
    return b - to_tdev(nb, order(b))
end

function Base.:*(nb::Number, b::TDev)
    return TDev(nb * b.dev)
end

Base.:*(b::TDev, nb::Number) = nb * b

Base.:/(a::TDev, nb::Number) = TDev(a.dev ./ nb)

function Base.:^(a::TDev, nb::Integer)
    if nb == 1
        return a
    elseif nb == 0
        return TDev(eltype(a), order(a))
    end
end

# BASE FUNCTIONS
function Base.sin(a::TDev)
end


# IN-PLACE OPERATIONS
function add!(target::TDev, b::TDev)
    #@assert order(target) == order(b)
    target.dev .= target.dev + b.dev
    return
end