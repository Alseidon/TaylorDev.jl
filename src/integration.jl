using .TDevCore

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
            @info "finished" i
            return pt
        end
    end
    @error "couldn't finish" pt t
    @assert false
end