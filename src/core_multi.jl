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

#=
struct MTDev{S1, S2, T, L} <: Number
    dev::MMatrix{S1, S2, T, L}
    
    MTDev(m::A) where {A<:Union{Matrix,SMatrix,MMatrix}} = (
            new{size(m)..., eltype(m), prod(size(m))}(
                MMatrix{size(m)..., eltype(m)}(m)
            ))
end

MTDev(ndims::Int, order::Int) = MTDev(zeros(ndims, order+1))
MTDev(t::DataType, ndims::Int, order::Int) = MTDev(zeros(t, (ndims, order+1)))
=#

end