module MTDevCore
export MTDev, get_coeffs, lower_order, highest_common_order

using StaticArrays

struct MTDev{S, T, N, L} <: Number
    dev::MArray{S, T, N, L}
    MTDev(v::A) where {A<:AbstractArray} = (
        new{Tuple{size(v)...}, eltype(v), ndims(v), prod(size(v))}(
            MArray{Tuple{size(v)...}, eltype(v)}(v)
        ))
    MTDev{S, T}() where {S, T} = (
        sizer = Size(S); new{S, T, length(sizer), prod(sizer)}(
            MArray{S, T}(zeros(T, Tuple(sizer)))
        )
    )
    MTDev{S, T}(v::A) where {S, T, A<:AbstractArray} = (
        sizer = Size(S); new{S, T, length(sizer), prod(sizer)}(
            MArray{S, T}(v)
        )
    )
     # use this better
end

MTDev(orders::Tuple) = MTDev(zeros(orders .+ 1))
MTDev(t::DataType, orders::Tuple) = MTDev(zeros(t, orders .+ 1))


get_coeffs(mtdev::MTDev{N}) where {N} = collect(Iterators.product(axes(mtdev.dev)...)) :: MArray{N}

lower_order(mtdev::MTDev{N}, S) where {N} = MTDev(
    view(mtdev.dev, map(i->1:i, S)...)
)

highest_common_order(::Type{N}, ::Type{M}) where {N, M} = Tuple(
    map(min, Tuple(Size(N)), Tuple(Size(M)))
)


end