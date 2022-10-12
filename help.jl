"""
Helper functions.
"""

using JLD2
include("kepl.jl")

LL = [5.452627468351345,
    -2.7573257341445694,
    -0.07602563473606634,
    2.2832251849342495,
    -0.5011496788005338,
    -0.16672728969322886,
    -0.13844909489132098,
    -0.5998550692761098,
    -0.08943355865879812]

# Iteration delta values
" Calculate a Newton method's iteration delta value."
function d_newton(E, e, M)
    df2 = e * sin(E)
    df3 = e * cos(E)
    f = (E - df2) - M
    df = 1 - df3
    return -f / df
end

" Calculate a Newton method's iteration delta value."
function d_halley(E, e, M)
    df2 = e * sin(E)
    df3 = e * cos(E)
    f = (E - df2) - M
    df = 1 - df3
    return -2 * f * df / (2 * df^2 - f * df2)
end

# Iteration
"Perform `n` iterations of specified `delta_f` method."
function iter(n, e, M, E0, delta_f, trace=false)
    for i in 1:n
        E0 = E0 + delta_f(E0, e, M)
        if trace println(Float64(E0)) end
    end
    return E0
end

"Perform 2 Newton iterations being the second one a simplified version."
function newton_simple2(e, M, E0)
    delta, eta, nu = deltaetanu(e, M, E0)
    aux1 = eta * delta
    aux2 = 1 + aux1 / 2 - nu * delta ^ 2 / 6
    d = delta * aux2 / (1 + aux1)
    E = E0 + d
    return E
end

# Grid operations
"Reduce grid by the operation `op` (max, sum) ignoring `NaN` values."
function reduce_grid(g, op)
    mxi, sma = 0, 0   # maximum, summation
    for i in g
        if !isnan(i)
            sma += i
            if i > mxi mxi = i end
        end
    end

    if op == "sum"
        return sma
    elseif op == "max"
        return mxi
    else
        error(op, " not supported")
    end
end

# Error
errorf(x, y) = abs(x - y)
errorf_00(x, y) = (x - y)
errorf_01(x, y) = (x - y)^2
errorf_02(x, y) = (x - y)^2 / y
errorf_03(x, y) = abs(x - y)
errorf_04(x, y) = abs(x - y) / y
errorf_05(x, y) = (x - y)^4

# Helpers
"Join 2 ranges in array."
function join_ranges(start, n1, stop1, n2, stop2)
    stp1 = stop1 - (stop1 - start) / n1
    r1 = range(start, stp1, n1)
    r2 = range(stop1, stop2, n2)
    arr = vcat(r1, r2)
    return arr
end

"Create array of `e` values with more close to value 1."
function variable_ee(ne)
    ee = Array{Float64}(undef, ne)
    for i in 1:ne
        ee[i] = sqrt(sin((i - 1) * pi / (2 * ne - 2)))
    end
    return ee
end

"Create array of `M` values with more close to value 1."
function variable_MM(nM)
    ee = variable_ee(nM)
    ee1 = 1 .- ee
    return reverse(ee1) .* pi
end

"Coefficients for cubic approximation."
function cubic_ab(e, M)
    a = e * sin(M)
    b = e * cos(M)
    return a, b
end

"Root of depressed cubic x^3 + 3qx -2r = 0 with E0 = M + x."
function cubic_depressed(r, q, p)
#    if trace
#        c = 3 * q
#        d = - 2 * r
#        Δ = -( 4 * c ^ 3 + 27 * d ^ 2)
#        println("r = ", r, "\nq = ", q, "\np = ", p)
#        println("f(x) = a*x^3 + c*x + d")
#        println("a = 1\nc = ", c, "\nd = ", d, "\nΔ = ", Δ, "\n")
#        #f(E0) = (E0 - p)^3 + 3 * q * (E0 - p) - 2 * r
#        #df(E0) = 3 * (E0 - p)^2 + 3 * q
#    end

    # Filter by equivalent to discriminant if it has 3 roots (Δ = -108 * aux).
    aux = r^2 + q^3
    if aux < 0.  return NaN
    else
        if aux == 0. && q == 0.     # f(x) = x^3
            y = 0
        else                        # f(x) = x^3 - kx | f(x) = x^3 +cx + d
            s = cbrt((sqrt(aux) + r))
            y = 2 * r / (s^2 + q + (q / s)^2)
        end
        x = p + y
        return x
    end
end

"Print Anstaz 2 λ values conditions."
function lambda_conditions(lambdas)
    l1, l2, l3, l4, l7, l8, l9, l10, l11, l12, l13, l14 = lambdas
    l0 = 3 - l1 - l2
    l5 = l7 + 2
    l6 = -2 * (1 + l7)
    println("λ0 + λ1 + λ2 = 3 \t -->\t", round(l0 + l1 + l2, digits=3))
    println("λ5 + λ6 + λ7 = 0 \t -->\t", round(l5 + l6 + l7, digits=3))
    println("λ6 + 2 * λ7 = -2 \t -->\t", round(l6 + 2 * l7, digits=3))
    println("λ10 + λ11 + λ12 = -1 \t -->\t", round(l10 + l11 + l12, digits=3))
end

"Create corner grid for ee = [emin, 1.] & MM = [0., Mmax] focusing at [1., 0]."
function create_corner(n, emin, Mmax)
    ne, nM = n + 1, n
    ee = variable_ee(ne)
    MM = variable_MM(nM)

    ee = (ee .* (1 - emin)) .+ emin
    MM = MM .* (Mmax / pi)
    return ee, MM
end

"Fill grid values to go from parameter vector to grid (ee to e_g)."
function fill_grid(p1, p2)
    np1, np2 = length(p1), length(p2)
    p1_g = Array{Float64, 2}(undef, (np1, np2))
    p2_g = Array{Float64, 2}(undef, (np1, np2))

    for i in 1:np2
        p1_g[:, i] = p1
    end
    for i in 1:np1
        p2_g[i, :] = p2
    end
    return p1_g, p2_g
end

# Saving
"Calculate a mixed resolution real solution grid."
function calculate_realsol()
    ee = join_ranges(0, 1000, .7, 2001, 1)
    MM = join_ranges(0, 2000, .6*pi, 1000, pi)
    sol_g = Array{Float64, 2}(undef, (length(ee), length(MM)))
    apply_solver!(sol_g, realsol, ee, MM)
    return sol_g, ee, MM
end

"Save calculated real solution grid."
function save_realsol(fname, sol_g, ee, MM)
    fpath = string("sols/", fname, ".jld2")
    jldsave(fpath; sol_g, ee, MM)
end
