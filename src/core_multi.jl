module MTDevCore
export MTDev, get_coeffs

using StaticArrays

struct MTDev{S, T, N, L} <: Number
    dev::MArray{S, T, N, L}
    MTDev(v::A) where {A<:Union{Array,SArray,MArray}} = (
        new{Tuple{size(v)...}, eltype(v), ndims(v), prod(size(v))}(
            MArray{Tuple{size(v)...}, eltype(v)}(v)
        ))
    MTDev{S, T}() where {S, T} = (
        sizer = Size(S); new{S, T, length(sizer), prod(sizer)}(
            MArray{S, T}(zeros(T, Tuple(sizer)))
        )
    ) # use this better
end

MTDev(orders::Tuple) = MTDev(zeros(orders .+ 1))
MTDev(t::DataType, orders::Tuple) = MTDev(zeros(t, orders .+ 1))


get_coeffs(mtdev::MTDev{N}) where {N} = collect(Iterators.product(axes(mtdev.dev)...)) :: MArray{N}
#= 
lower_order(mtdev::MTDev{N}, ::Type{M}) where {N, M} = MTDev(
    mtdev.dev[:]
)
=#

end