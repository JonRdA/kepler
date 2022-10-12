"""
Kepler's equation solver methods.

Evaluate in grids for different values of `e` and `M`.
"""

# Starters
"RAE's best starter, fails for e = 1, M = 0, it has no conditions."
function starter(e, M)
    e1 = 1 - e
    q = 2 * e1 / e
    r = 3 * M / e
    s = (sqrt(r^2 + q^3) + r)^(1 / 3)
    E00 = 2 * r / (s^2 + q + (q / s)^2)
    E0 = (M^2 + (pi - M) * E00) / pi
end

# Variable starters
"Calculate Ansatz 1 `r`, `q` and `p` coefficients."
function ansatz_coefs(e, M, lambdas)
    l0, l1, l2, l3, l4, l5, l6, l7, l8 = lambdas
    a, b = cubic_ab(e, M)
    r = (l0 + l1 * b) / (1 + l2 * b) * a
    q = ((l3 + l4 * b) * (1 - b) + (l5 + l6 * b) * a ^ 2) / (1 + l2 * b)
    p = (l7 + l8 * b) / (1 + l2 * b) * a
    return r, q, p
end

"Ansatz 1 variable starter depending on `lambdas` vector of 9 elements."
function vstarter(e, M, lambdas, trace=false)
    r, q, p = ansatz_coefs(e, M, lambdas)
    E0 = cubic_depressed(r, q, p) + M
    return E0
end

"Calculate Ansatz 2 `r`, `q` and `p` coefficients."
function ansatz2_coefs_expanded(e, M, lambdas)
    l1, l2, l3, l4, l5, l6, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = lambdas

    a, b = cubic_ab(e, M)
    a2 = a^2
    a3 = a^3
    b2 = b^2
    b3 = b^3

    l8 = fma(-3, l10, fma(-2, l9, -2))
    l7 = -l9 - l8 - l10
    l0 = 3 - l1 - l2 - l3

    r = fma(fma(l3, b3, fma(l2, b2, fma(l1, b, l0))), a, (fma(l6, b2, fma(l5, b, l4)) * a3))
    q = fma(fma(l13, b2, fma(l12, b, l11)), a2, fma(l10, b3, fma(l9, b2, fma(l8, b, l7))))
    p = fma(fma(l17, b3, fma(l16, b2, fma(l15, b, l14))), a, fma(l20, b2, fma(l19, b, l18)) * a3)
    return r, q, p
end

"Calculate Ansatz 2 `r`, `q` and `p` coefficients."
function ansatz2_coefs_reduced(e, M, lambdas)
    l1, l2, l3, l4, l7, l8, l9, l10, l11, l12, l13, l14 = lambdas
    l0 = 3 - l1 - l2
    l5 = l7 + 2
    l6 = -2 * (1 + l7)
    a, b = cubic_ab(e, M)
    r = (l0 + l1 * b + l2 * b^2) * a + (l3 + l4 * b) * a^3
    q = (l5 + l6 * b + l7 * b^2) + (l8 + l9 * b) * a^2
    p = (l10 + l11 * b + l12 * b^2) * a + (l13 + l14 * b) * a^3
    return r, q, p
end

"Ansatz 2 variable starter with `lambdas` vector of 15 elements."
function vstarter2(e, M, lambdas)
    r, q, p = ansatz2_coefs(e, M, lambdas)
    E0 = cubic_depressed(r, q, p) + M
    return E0
end

# Solvers
"Kepler equation RAE procedure (EKEPL), `NaN` at e = 0, M = 0 due to starter."
function ekepl(e, M)
    E = starter(e, M)
    for i in 1:2
        df2 = e * sin(E)
        df3 = e * cos(E)
        f = (E - df2) - M
        df = 1 - df3
        
        # Halley's delta
        d = f * df / (f * df2 / 2 - df^2)

        # Rectified f and f'
        f = f + d * (df + d / 2 * (df2 + d * df3 / 3))
        df = df + d * (df2 + d * df3 / 2)

        # Halley's delta + Newton's
        E = E + d - f / df
    end
    return E
end

"Solve following Ander Murua iteration procedure. Uses RAE starter `starter`."
function kepler2(e, M)
    E = starter(e, M)
    for i in 1:2
        df2 = e * sin(E)
        df3 = e * cos(E)
        f = E - df2 - M
        df = 1 - df3

        d = -f / df
        a = df2 / df
        b = df3 / df

        E = E + d *(1 + d / 2 * (d * (a^2 - b/3) - a))
    end
    return E
end

"Solve Kepler's equation for e and E."
function M_solver(e, M, E)
    return E - e * sin(E)
end

"Solve Kepler's equation iterating until desired tolerance `tol` reached."
function realsol(e, M, tol=1e-50)::BigFloat
    rrM = mod2pi(M)     # range reduce M
    if rrM < -pi rrM += 2 * pi end
    if rrM > pi rrM -= 2 * pi end
    if e == 0. || rrM == 0.  return M end

    E = BigFloat(vstarter(e, rrM,
        [5.452, -2.754, -0.074, 2.283, -0.508, -0.166, -0.138, -0.598, -0.082]))
    for i in 1:100
        E0 = E
        E = iter(10, e, rrM, E, d_halley)
        if abs(E0 - E) < tol
            return E + M - rrM
        end
    end
    error(string("No exact solution found at [e, M]: ", [e M]))
end

# Variable solvers
"Solve using 2 Newton simplified iterations, variable starter 1, Ansatz 1."
function vkepler(e, M, lambdas)
    E0 = vstarter(e, M, lambdas)
    return newton_simple2(e, M, E0)
end

"Solve using 2 Newton simplified iterations, variable starter 2, Ansatz 2."
function vkepler2(e, M, lambdas)
    E0 = vstarter2(e, M, lambdas)
    return newton_simple2(e, M, E0)
end

# Wrappers
"Construct grid of pairs of `(e, M)` values and apply function `f` to them."
function apply_solver!(g, f::Function, e_g, M_g)
    @assert size(g) == size(e_g) == size(M_g) "Wrong grid size."
    for i in 1:length(e_g)
        g[i] = f(e_g[i], M_g[i])
    end
end

"Construct grid of `(e, M)` pairs and apply variable starter to them."
function apply_vsolver!(g, f::Function, e_g, M_g, lambdas)
    @assert size(g) == size(e_g) == size(M_g) "Wrong grid size."
    for i in 1:length(e_g)
        g[i] = f(e_g[i], M_g[i], lambdas)
    end
end

"Apply function `f` to grids `g1` and `g2`, write results in new grid (matrix)."
function apply_grids!(g, f::Function, g1, g2)
    @assert size(g1) == size(g2) "Size of grids does not match."
    for i in 1:length(g1)
        g[i] = f(g1[i], g2[i])
    end
end

"Apply function `f` to grid."
function apply_grid!(g, f::Function, e_g, M_g, E_g)
    @assert size(g) == size(e_g) == size(M_g) == size(E_g) "Wrong grid size."
    for i in 1:length(g)
        g[i] = f(e_g[i], M_g[i], E_g[i])
    end
end


function kepler(e, M, lambdas)
    # Lambda parameter values

    # Variable change
    a = e * sin(M)
    b = e * cos(M)

    a2 = a^2
    a3 = a^3
    b2 = b^2
    b3 = b^3

    l1, l2, l3, l4, l5, l6, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = lambdas
    
    # Calculate r, q, p
    l8 = fma(-3, l10, fma(-2, l9, -2))
    l7 = -l9 - l8 - l10
    l0 = 3 - l1 - l2 - l3
    r = fma(fma(l3, b3, fma(l2, b2, fma(l1, b, l0))), a, (fma(l6, b2, fma(l5, b, l4)) * a3))
    q = fma(fma(l13, b2, fma(l12, b, l11)), a2, fma(l10, b3, fma(l9, b2, fma(l8, b, l7))))
    p = fma(fma(l17, b3, fma(l16, b2, fma(l15, b, l14))), a, fma(l20, b2, fma(l19, b, l18)) * a3)

    # Solve cubic
    aux = r^2 + q^3
    s = cbrt((sqrt(aux) + r))
    y = 2 * r / (s^2 + q + (q / s)^2)
    x = p + y
    x2 = x ^ 2

    # x - sin(x) Taylor polynomial, level 18
    k1 = 0.16666666666666666
    k2 = -0.008333333333333333
    k3 = 0.0001984126984126984
    k4 = -2.7557319223985893e-6
    k5 = 2.505210838544172e-8
    k6 = -1.6059043836821613e-10
    k7 = 7.647163731819816e-13
    k8 = -2.8114572543455206e-15
    k9 = 8.22063524662433e-18
    v8 = fma(x2, k9, k8)
    v7 = fma(x2, v8, k7)
    v6 = fma(x2, v7, k6)
    v5 = fma(x2, v6, k5)
    v4 = fma(x2, v5, k4)
    v3 = fma(x2, v4, k3)
    v2 = fma(x2, v3, k2)
    v1 = fma(x2, v2, k1)
    x_sinx = v1 * x * x2

    # 1 - cos(x) Taylor polynomial, level 18
    w1 = 0.5
    w2 = -0.041666666666666664
    w3 = 0.001388888888888889
    w4 = -2.48015873015873e-5
    w5 = 2.755731922398589e-7
    w6 = -2.08767569878681e-9
    w7 = 1.1470745597729725e-11
    w8 = -4.779477332387385e-14
    w9 = 1.5619206968586225e-16
    v8 = fma(x2, w9, w8)
    v7 = fma(x2, v8, w7)
    v6 = fma(x2, v7, w6)
    v5 = fma(x2, v6, w5)
    v4 = fma(x2, v5, w4)
    v3 = fma(x2, v4, w3)
    v2 = fma(x2, v3, w2)
    v1 = fma(x2, v2, w1)
    bat_cosx = v1 * x2

    # 1 - b or 1 - e * cos(M)
    bat_b = 1 - e + e * (2 * (sin(M / 2))^2)

    # Iterative method
    sinx = sin(x)
    cosx = cos(x)
    f = (bat_b) * x - a + a * bat_cosx + b * x_sinx
    df = (bat_b) + b * bat_cosx + a * x - a * x_sinx
    df2 = b * sinx + a * cosx
    df3 = b * cosx - a * sinx

    delta = - f / df
    eta = df2 / df
    nu = df3 / df

    aux1 = eta * delta
    aux2 = 1 + aux1 / 2 - nu * delta ^ 2 / 6
    d = delta * aux2 / (1 + aux1)
    E = M + x + d
    return E
end
