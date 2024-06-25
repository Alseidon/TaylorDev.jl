using .TDevCore, .MTDevCore

function derivative(f, val)
    return f(val + epsilon()).dev[2]
end


function gradient(
		f::Function, pt::AbstractVector,
		::Type{T}, params=nothing) where T
    res = zeros(T, length(pt))
	if params === nothing
		for i in eachindex(res)
			res[i] = f(pt + [j==i ? epsilon(eltype(pt)) : zero(eltype(pt)) for j in eachindex(pt)]).dev[2]
		end
	else
		for i in eachindex(res)
			res[i] = f(pt + [j==i ? epsilon(eltype(pt)) : zero(eltype(pt)) for j in eachindex(pt)], params).dev[2]
		end
	end
    return res
end


function hamiltonian_vec_field(
		ham::Function, p::AbstractVector, q::AbstractVector,
		::Type{T}, params=nothing) where T
	n = length(p)
	vf = gradient(ham, vcat(p, q), T, params)
	vf[1:n], vf[(n+1):end] = -vf[(n+1):end], vf[1:n]
	return vf
end