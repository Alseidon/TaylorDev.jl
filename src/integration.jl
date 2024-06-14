using .TDevCore

# SINGLE VARIABLE
function taylor_method_step!(f::Function, initial_pt::Number,
    initial_t::Number, final_t::Number,
    err::Number, ord::Int)
    x = to_tdev(initial_pt, ord)
    aux = zero(x)
    for i in 1:ord
        aux = f(x, initial_t)
        x.dev[i+1] = aux.dev[i] / i
    end
    max_step = final_t - initial_t
    if max_step == 0
        return compute(x, 0), initial_t, true
    end
    step = (err / x.dev[ord+1])^(1 / ord) * sign(max_step)
    if abs(step) > abs(max_step)
        return compute(x, max_step), initial_t + max_step, true
    else
        return compute(x, step), initial_t + step, false
    end
end

function taylor_method(f::Function, initial_pt::Number,
    initial_t::Number, final_t::Number;
    err::Number=1e-16, ord::Int=15, max_steps::Int=100)
    pt = copy(initial_pt)
    t = initial_t
    for i in 1:max_steps
        pt, t, is_over = taylor_method_step!(f, pt, t, final_t, err, ord)
        if is_over
            return pt
        end
    end
    @error "couldn't finish" pt t
    @assert false
end


# MULTIPLE VARIABLES
function taylor_method_step!(f::Function, initial_pt::AbstractVector,
    initial_t::Number, final_t::Number,
    err::Number, ord::Int)
    x = to_tdev.(initial_pt, ord)
    aux = zero.(x)
    for i in 1:ord
        aux = f(x, initial_t)
        for j in eachindex(x)
            x[j].dev[i+1] = aux[j].dev[i] / i
        end
    end
    max_step = final_t - initial_t
    if max_step == 0
        return compute.(x, 0), initial_t, true
    end
    # TODO compute order
    step = (err / abs(sum(map(c->c.dev[ord+1], x))))^(1 / ord) * sign(max_step)
    if abs(step) > abs(max_step)
        return compute.(x, max_step), initial_t + max_step, true
    else
        return compute.(x, step), initial_t + step, false
    end
end

function taylor_method(f::Function, initial_pt::AbstractVector,
    initial_t::Number, final_t::Number;
    err::Number=1e-16, ord::Int=15, max_steps::Int=100)
    pt = copy(initial_pt)
    t = initial_t
    for i in 1:max_steps
        pt, t, is_over = taylor_method_step!(f, pt, t, final_t, err, ord)
        if is_over
            return pt
        end
    end
    @error "couldn't finish" pt t
    @assert false
end

function taylor_method(f::Function, initial_pt::AbstractVector,
    times::AbstractVector;
    err::Number=1e-16, ord::Int=15, max_steps::Int=100)
    res = zeros(eltype(initial_pt), (length(times), length(initial_pt)))
    previous_t = zero(eltype(times))
    previous_pt = copy(initial_pt)
    for i in eachindex(times)
        res[i,:] = taylor_method(f, previous_pt, previous_t, times[i];
                err=err, ord=ord, max_steps=max_steps)
        previous_pt[:] = res[i,:]
        previous_t = times[i]
    end
    return res
end