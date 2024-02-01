#include("core.jl")

using .TDevCore, .MTDevCore

const subscript_digits = [c for c in "₀₁₂₃₄₅₆₇₈₉"]

const superscript_digits = [c for c in "⁰¹²³⁴⁵⁶⁷⁸⁹"]

function subscriptify(n::Int)
    dig = reverse(digits(n))
    return join([subscript_digits[i+1] for i in dig])
end

function superscriptify(n::Int)
    if n == 1
        return ""
    end
    dig = reverse(digits(n))
    return join([superscript_digits[i+1] for i in dig])
end

function inherits_from(t1::Type, t2::Type)
    if t1 == Any
        return t2 == Any
    else
        return t1 == t2 || inherits_from(supertype(t1), t2)
    end
end

function inherits_from(t1::Type, t2s::Array)
    if t1 == Any
        return Any ∈ t2s
    else
        return t1 ∈ t2s || inherits_from(supertype(t1), t2s)
    end
end

function stringify(nb::Number)
    type = typeof(nb)
    if inherits_from(type, [Real, Integer])
        return string(nb)
    else
        return "(" * string(nb) * ")"
    end
end

function pretty_print(tdev::TDev)
    type = eltype(tdev)
    #=if order(tdev) == 0
        return ""
    end=#
    str = stringify(tdev.dev[1])
    for i in 1:order(tdev)
        if tdev.dev[i+1] != 0
            if hasmethod(isless, (eltype(tdev), Int)) && tdev.dev[i+1] < 0
                str *= " - " *
                       stringify(-tdev.dev[i+1]) * "ε" * superscriptify(i)
            else
                str *= " + " * stringify(tdev.dev[i+1]) * "ε" * superscriptify(i)
            end
        end
    end
    str *= " + 𝒪(ε" * superscriptify(order(tdev) + 1) * ")"
    return str
end

function pretty_print(mtdev::MTDev)
    type = eltype(mtdev)
    ord = orders(mtdev)
    ndimensions = ndims(mtdev)
    str = stringify(mtdev.dev[1])
    for tp in get_coeffs(mtdev)[2:end]
        val = mtdev.dev[tp...]
        if val != 0
            if hasmethod(isless, (type, Int)) && val < 0
                str *= " - " *
                    stringify(-val)
            else
                str *= " + " * stringify(val)
            end
            for i in 1:ndimensions
                if tp[i] > 1
                    str *= "ε" * subscriptify(i) * superscriptify(tp[i] - 1)
                end
            end
        end
    end
    str *= " + 𝒪("
    for i in 1:ndimensions
        if i != 1
            str *= ","
        end
        str *= "ε" * subscriptify(i) * superscriptify(ord[i] + 1)
    end
    str *= ")"
    return str
end

function Base.show(io::IO, a::T) where T <: Union{TDev, MTDev}
    return print(io, pretty_print(a))
end