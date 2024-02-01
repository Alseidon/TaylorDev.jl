module MTDevCore
export MTDev #, order, to_tdev, compute, deriv, epsilon

using StaticArrays

struct MTDev{S, T, N, L} <: Number
    dev::MArray{S, T, N, L}
    MTDev(v::A) where {A<:Union{Array,SArray,MArray}} = (
        new{Tuple{size(v)...}, eltype(v), ndims(v), prod(size(v))}(
            MArray{Tuple{size(v)...}, eltype(v)}(v)
        ))
end

MTDev(orders::Tuple) = MTDev(zeros(orders .+ 1))
MTDev(t::DataType, orders::Tuple) = MTDev(zeros(t, orders .+ 1))

end