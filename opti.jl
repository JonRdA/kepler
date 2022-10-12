"""
Methods for starter search & evaluation with optimization techniques.
"""

# Ansatz 1 - Starter error with real solution
"""
Compute the cost for lambda optimization by cuadratic error sum of starter.

Only parameters is `lambdas` length 9 array of `Float64`, due to package optim.
The following variables must be set in function scope: 
- `ee` iterator of `e` values.
- `MM` iterator of `M` values.
- `sol_g:Array{Float64, 2}` real solution grid for `ee` and `MM` values.
- `starter_g:Array{Float64, 2}` starter grid for `ee` and `MM` values.
- `error_g:Array{Float64, 2}` error grid for `ee` and `MM` values.

Catches 3 solution cubics and penalizes them with a error val of `1e6`.
"""
function cost(lambdas)
    apply_vsolver!(starter_g, vstarter, ee, MM, lambdas)
    apply_grids!(error_g, errorf_01, starter_g, sol_g)
    replace!(error_g, NaN=>1e6)
    c = reduce_grid(error_g, "sum")
    return c
end

# Ansatz 1 - Starter evaluation with ϵ values
"Calculate ϵ1, ϵ2, ϵ3 from starter value"
function epsilon(E0, e, M)
    c = e * cos(E0)
    if c == 1. && M == 0. return 0., 0., 0. end
    s = e * sin(E0)
    aux = 1 - c
    d1 = (E0 - s - M) / aux
    e1 = abs(d1)
    e2 = abs((s * d1) / aux)
    e3 = sqrt(abs(c) * d1^2 / aux)
    return e1, e2, e3
end

"Calculate the sum of ϵ1, ϵ2, ϵ3 from entire grid."
function epsilon_sum(g, ee, MM)
    e1sum, e2sum, e3sum = 0., 0., 0.        # epsilon values
    ne, nM = length(ee), length(MM)
    for j in 1:nM
        for i in 1:ne
            e1, e2, e3 = epsilon(g[i, j], ee[i], MM[j])
            e1sum += e1
            e2sum += e2
            e3sum += e3
        end
    end
    return e1sum, e2sum, e3sum
end

"Calculate maximum of ϵ1, ϵ2, ϵ3 from entire grid."
function epsilon_max(g, ee, MM)
    e1max, e2max, e3max = 0., 0., 0.        # epsilon values
    ne, nM = length(ee), length(MM)
    for j in 1:nM
        for i in 1:ne
            e1, e2, e3 = epsilon(g[i, j], ee[i], MM[j])
            if e1 > e1max e1max = e1 end
            if e2 > e2max e2max = e2 end
            if e3 > e3max e3max = e3 end
        end
    end
    return e1max, e2max, e3max
end

"Calculate ϵ1, ϵ2, ϵ3 values for entire grid and save on `ep_grids`."
function epsilon_grids!(g, epsilon_gs, ee, MM)
    e1g, e2g, e3g = epsilon_gs
    ne, nM = length(ee), length(MM)
    for i in 1:ne
        for j in 1:nM
            e1, e2, e3 = epsilon(g[i, j], ee[i], MM[j])
            e1g[i, j] = e1
            e2g[i, j] = e2
            e3g[i, j] = e3
        end
    end
end

"""
Compute the cost for lambda optimization by epsilon minimization of starter.

Only parameters is `lambdas` length 9 array of `Float64`, due to package optim.
The following variables must be set. 
- `ee` iterator of `e` values.
- `MM` iterator of `M` values.
- `starter_g:Array{Float64, 2}` starter grid for `ee` and `MM` values.
- `epsilon_reduce::Function` either `epsilon_sum` or `epsilon_max`.

Catches 3 solution cubics and penalizes them with a error val of `1e6`.
"""
function cost2(lambdas)
    apply_vsolver!(starter_g, vstarter, ee, MM, lambdas)
    replace!(starter_g, NaN=>1e3)
    eps = epsilon_reduce(starter_g, ee, MM)
    #return eps[1]/8 + eps[2] + eps[3]
    return eps[1] ^ 2 + eps[2] ^ 2 + eps[3] ^ 2
end

# Ansatz 1 - Starter evaluation with γ
"Create needed objects for `cost3` function evaluation."
function prep_cost3(ee, MM)
    ne, nM = length(ee), length(MM)
    c_g = Array{Float64, 2}(undef, (ne, nM))
    E_g = Array{Float64, 2}(undef, (ne, nM))
    x_g = Array{Float64, 2}(undef, (ne, nM))
    cost_g = Array{Float64, 2}(undef, (ne, nM))
    apply_solver!(c_g, cij, ee, MM)
    apply_solver!(E_g, realsol, ee, MM)
    apply_solver!(x_g, xij, ee, MM)
    return c_g, E_g, x_g, cost_g
end

"Compute η and ν values for a given solution `E`."
function deltaetanu(e, M, E)
    if e == 1. &&  E == 0. return 0., 0., 0. end
    df3 = e * cos(E)
    df2 = e * sin(E)
    df = 1 - df3
    f = E - df2 - M
    delta = - f / df
    eta = df2 / df
    nu = df3 / df
    return delta, eta, nu
end

"Compute γ function value from `eta` and `nu`."
function gamma(eta, nu)
    return abs(eta / 24 + eta * nu / 4 - eta ^ 3 / 8)
end

"Compute as BigFloat the value of `c` for given `e` and `M`."
function cij(e, M, E)
    delta, eta, nu = deltaetanu(e, M, E)
    gm = gamma(eta, nu)
    return sqrt(gm)
end

"Compute as BigFloat the value of `x`(E - M) for given `e` and `M`."
function xij(e, M, E)
    return e * sin(E)
end

"Compute delta hat for ansatz coefs, checking Δ and reduced P(E)/f'(E)."
function deltahat(e, M, E, x, lambdas, coef_f, discriminant=true, reduced=false)
    r, q, p = coef_f(e, M, lambdas)
    aux = r^2 + q^3
    if discriminant
        if aux < 0.
            return NaN
        elseif aux == 0.
            return 0.
        end
    end
    if e == 1. && M == 0. return 0. end
    y = x - p
    aux1 = y ^ 3 + 3 * q * y - 2 * r
    reduced ?  aux2 = 1 - e * cos(E) : aux2 = 3 * y ^ 2 + 3 * q
    return aux1 / aux2
end

"""
Objective function for lambda by error estimation of solving procedure.

Estimates the error of a variable starter (lambdas) and 2 Newton iterations,
the second one being a simplified one.

Due to optim package function call scope must have variables:
- `reduced` bool to used reduced δ, P(E)/f'(E).
- `discriminant` bool to check 3 solution cubics and apply δ = 1e6.
- `ee` iterator of `e` values.
- `MM` iterator of `M` values.
- `c_g:Array{Float64, 2}` grid of c values (c = sqrt(γ)).
- `E_g:Array{Float64, 2}` grid of real solution E.
- `x_g:Array{Float64, 2}` grid of E - M values.
- `error_g:Array{Float64, 2}` error grid for `ee` and `MM` values.

So far catches 3 solution cubics and penalizes them with a `delta` val of `1e6`.
"""
function cost3(lambdas)
    cst = 0
    ne, nM = length(ee), length(MM)
    for j in 1:nM
        for i in 1:ne
            delta = deltahat(ee[i], MM[j], E_g[i, j], x_g[i, j], lambdas,
                            ansatz_coefs, discriminant, reduced)
            isnan(delta) ? cs = 1e6 : cs = c_g[i, j] * delta ^ 2
            cst += cs
        end
    end
    return cst
end

"Inplace version of `cost3` with parameters for obtaining grid and efficiency."
function cost3!(cost_g, ee, MM, c_g, E_g, x_g, lambdas,
        discriminant=true, reduced=false)
    ne, nM = length(ee), length(MM)
    for j in 1:nM
        for i in 1:ne
            delta = deltahat(ee[i], MM[j], E_g[i, j], x_g[i, j], lambdas,
                            ansatz_coefs, discriminant, reduced)
            isnan(delta) ? cs = 1e6 : cs = c_g[i, j] * delta ^ 2
            cost_g[i, j] = cs
        end
    end
end

# Ansatz 2 - Starter evaluation with γ
"""
Objective function for lambda (Ansatz2) by error estimation of solver.

Estimates the error of a variable starter (lambdas) and 2 Newton iterations,
the second one being a simplified one.

Due to optim package function call scope must have variables:
- `reduced` bool to used reduced δ, P(E)/f'(E).
- `discriminant` bool to check 3 solution cubics and apply δ = 1e6.
- `ee` iterator of `e` values.
- `MM` iterator of `M` values.
- `c_g:Array{Float64, 2}` grid of c values (c = sqrt(γ)).
- `E_g:Array{Float64, 2}` grid of real solution E.
- `x_g:Array{Float64, 2}` grid of E - M values.
- `error_g:Array{Float64, 2}` error grid for `ee` and `MM` values.

Catches 3 solution cubics and penalizes them with a `delta` val of `1e6`.
"""
function cost4(lambdas)
    cst = 0
    ne, nM = length(ee), length(MM)
    for j in 1:nM
        for i in 1:ne
            delta = deltahat(ee[i], MM[j], E_g[i, j], x_g[i, j], lambdas,
                             ansatz2_coefs, discriminant, reduced)
            isnan(delta) ? cs = 1e6 : cs = c_g[i, j] * delta ^ 2
            cst += cs
        end
    end
    return cst
end

"Inplace version of `cost4` with parameters for obtaining grid and efficiency."
function cost4!(cost_g, ee, MM, c_g, E_g, x_g, lambdas,
        discriminant=true, reduced=false)
    ne, nM = length(ee), length(MM)
    for j in 1:nM
        for i in 1:ne
            delta = deltahat(ee[i], MM[j], E_g[i, j], x_g[i, j], lambdas,
                             ansatz2_coefs, discriminant, reduced)
            isnan(delta) ? cs = 1e6 : cs = c_g[i, j] * delta ^ 2
            cost_g[i, j] = cs
        end
    end
end

# Ansatz 2 - Starter evaluation ee, EE grid
"Create needed objects for `cost5` function evaluation."
function prep_cost5(ee, EE)
    ne, nE = length(ee), length(EE)
    e_g, E_g = fill_grid(ee, EE)
    M_g = Array{Float64, 2}(undef, (ne, nE))
    c_g = Array{Float64, 2}(undef, (ne, nE))
    x_g = Array{Float64, 2}(undef, (ne, nE))
    cost_g = Array{Float64, 2}(undef, (ne, nE))
    apply_grid!(M_g, M_solver, e_g, e_g, E_g)
    apply_grid!(c_g, cij, e_g, M_g, E_g)
    apply_grid!(x_g, xij, e_g, M_g, E_g)

    return e_g, M_g, E_g, c_g, x_g, cost_g
end

"""
Objective function for lambda (Ansatz2) by error estimation of solver.

Estimates the error of a variable starter (lambdas) and 2 Newton iterations,
the second one being a simplified one.

Due to optim package function call scope must have variables:
- `reduced` bool to used reduced δ, P(E)/f'(E).
- `discriminant` bool to check 3 solution cubics and apply δ = 1e6.
- `e_g:Array{Float64, 2}` grid of `e` values.
- `M_g:Array{Float64, 2}` grid of `M` values.
- `E_g:Array{Float64, 2}` grid of real solution E.
- `x_g:Array{Float64, 2}` grid of x values.
- `c_g:Array{Float64, 2}` grid of x values.

Catches 3 solution cubics and penalizes them with a `delta` val of `1e6`.
"""
function cost5(lambdas)
    cst = 0
    for i in 1:length(e_g)
        delta = deltahat(e_g[i], M_g[i], E_g[i], x_g[i], lambdas,
                         ansatz2_coefs, discriminant, reduced)
        #isnan(delta) ? cs = 1e6 : cs = c_g[i, j] * delta ^ 2
        cs = c_g[i] * delta ^ 2
        cst += cs
    end
    return cst
end

function cost5!(cost_g, lambdas, e_g, M_g, E_g, c_g, x_g, discriminant, reduced)
    for i in 1:length(e_g)
        delta = deltahat(e_g[i], M_g[i], E_g[i], x_g[i], lambdas,
                         ansatz2_coefs, discriminant, reduced)
        #isnan(delta) ? cs = 1e6 : cs = c_g[i, j] * delta ^ 2
        cost_g[i] = c_g[i] * delta ^ 2
    end
end
