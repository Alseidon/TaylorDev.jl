module MTDevCore
export MTDev #, order, to_tdev, compute, deriv, epsilon

using StaticArrays

struct MTDev{S1, S2, T, L} <: Number
    dev::MMatrix{S1, S2, T, L}
    
    MTDev(m::A) where {A<:Union{Matrix,SMatrix,MMatrix}} = (
            new{size(m)..., eltype(m), prod(size(m))}(
                MMatrix{size(m)..., eltype(m)}(m)
            ))
end

MTDev(ndims::Int, order::Int) = MTDev(zeros(ndims, order+1))
MTDev(t::DataType, ndims::Int, order::Int) = MTDev(zeros(t, (ndims, order+1)))

end