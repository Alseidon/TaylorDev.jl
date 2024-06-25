module TDevCore
export TDev#, order, to_tdev, compute, deriv, epsilon

using StaticArrays

struct TDev{N,T} <: Number
    dev::MVector{N,T}
    TDev(v::A) where {A<:Union{Vector,SVector,MVector}} = (
        new{length(v),eltype(v)}(MVector{length(v),eltype(v)}(v)))
end

TDev(ord::Int) = TDev(zeros(ord + 1))
TDev(t::DataType, ord::Int) = TDev(zeros(t, ord + 1))

TDev{N, T}(x::Number) where {N, T} = TDev(
    [(i==1 ? convert(T, x) : zero(T)) for i in 1:N]
)

end