using .TDevCore

order(tdev::TDev{N}) where {N} = N - 1 #length(tdev.dev) - 1
Base.copy(tdev::TDev) = TDev(tdev.dev)
Base.eltype(tdev::TDev) = eltype(tdev.dev)
Base.zero(tdev::TDev) = TDev(eltype(tdev), order(tdev))
Base.zero(::Type{TDev{N, T}}) where {N, T} = TDev(T, N-1)
function Base.one(tdev::TDev)
    res = zero(tdev)
    res.dev[1] = one(eltype(res))
    return res
end
Base.one(::Type{TDev{N, T}}) where {N, T} = TDev([
    i == 1 ? one(T) : zero(0) for i in 1:N
])

epsilon(order::Int=1) = TDev([1*(i==2) for i in 1:(order+1)])
epsilon(t::Type, order::Int=1) = TDev([(i==2 ? one(t) : zero(t)) for i in 1:(order+1)])

function lowest_nonzero_order(t::TDev)
    for i in eachindex(t.dev)
        if t.dev[i] != 0
            return i-1
        end
    end
    return order(t) + 1
end

constant_term(a::TDev) = a.dev[1]

function to_tdev(nb::Number, ord::Int)
    t = TDev([(i==1 ? nb : zero(nb)) for i in 1:(ord+1)])
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

function Base.:+(a::TDev{N, T}, b::TDev{M, U}) where {N, T, M, U}
    return TDev( a.dev[1:min(N, M)] + b.dev[1:min(N, M)] )
end

function Base.:-(a::TDev{N, T}, b::TDev{M, U}) where {N, T, M, U}
    return TDev( a.dev[1:min(N, M)] - b.dev[1:min(N, M)] )
end

function Base.:-(a::TDev)
    return TDev(-a.dev)
end

function Base.:*(a::TDev, b::TDev)
    mina = lowest_nonzero_order(a)
    minb = lowest_nonzero_order(b)
    ord1 = order(a) + minb
    ord2 = order(b) + mina
    if ord2 < ord1
        return b*a
    end
    res = TDev(promote_type(eltype(a), eltype(b)), ord1)
    for o in (mina+minb):ord1#(mina+1):length(a.dev)
        for i in max(mina, o-order(b)):min(order(a), o-minb)#(minb+1):(order(b)-i+2)
            res.dev[o+1] += a.dev[i+1] * b.dev[o-i+1]
        end
    end
    return res
end

function Base.:/(a::TDev, b::TDev)
    @assert constant_term(b) != 0
    res = copy(a)#TDev(eltype(a), order(a))
    for i in eachindex(res.dev) # k=i-1; c_k = a_k
        for j in 2:i
            res.dev[i] -= b.dev[j] * res.dev[i-j+1]
        end
        res.dev[i] /= constant_term(b)
    end
    return res
end

# INTERACTION WITH NUMBERS
Base.:(==)(nb::Number, b::TDev) = false
Base.:(==)(b::TDev, nb::Number) = false

Base.isless(a::TDev, nb::Number) = constant_term(a) < nb
Base.isless(nb::Number, a::TDev) = nb < constant_term(a)


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

Base.:/(nb::Number, a::TDev) = to_tdev(nb, order(a)) / a

#= function Base.:^(a::TDev, nb::Integer)
    if nb == 1
        return a
    elseif nb == 0
        return TDev(eltype(a), order(a))
    end
end =#

# BASE FUNCTIONS
function Base.abs(a::TDev)
    @assert constant_term(a) != 0.
    return TDev(map(i->(i==1 ? abs(a.dev[1]) : a.dev[i]), 1:(order(a)+1)))
end

function Base.promote_rule(nb::Number, a::TDev)
    return (to_tdev(nb, order(a)), a)
end

function Base.sin(a::TDev)
    res = to_tdev(sin(constant_term(a)), order(a))

    remain = copy(a)
    remain.dev[1] = zero(eltype(remain))
    remain_i = copy(remain)
    vals = [sin(constant_term(a)), cos(constant_term(a)),
            -sin(constant_term(a)), -cos(constant_term(a)),]

    for i in 1:(order(a))
        res += vals[i%4 + 1] * remain_i / factorial(i)
        remain_i *= remain
    end
    return res
end

function Base.cos(a::TDev)
    res = to_tdev(cos(constant_term(a)), order(a))

    remain = copy(a)
    remain.dev[1] = zero(eltype(remain))
    remain_i = copy(remain)
    vals = [cos(constant_term(a)), -sin(constant_term(a)),
            -cos(constant_term(a)), sin(constant_term(a)),]

    for i in 1:(order(a))
        res += vals[i%4 + 1] * remain_i / factorial(i)
        remain_i *= remain
    end
    return res
end

function Base.sqrt(a::TDev)
    res = to_tdev(sqrt(constant_term(a)), order(a))

    remain = copy(a)
    remain.dev[1] = zero(eltype(remain))
    remain_i = copy(remain)
    val = .5/sqrt(constant_term(a))

    for i in 1:(order(a))
        res += val * remain_i / factorial(i)
        val *= -(2*(i-1)+1)/2 / constant_term(a)
        remain_i *= remain
    end
    return res
end

function sin_exp(val, k)
    x = BigFloat(val % 2π)
    res = zero(x)
    for i in 0:(k ÷ 2)
        res += (-1)^i * x^(2i+1) / factorial(big(2i+1))
    end
    return res
end

# IN-PLACE OPERATIONS
function add!(target::TDev, b::TDev)
    #@assert order(target) == order(b)
    target.dev .= target.dev + b.dev
    return
end