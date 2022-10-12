"""
Result analysis and evalutation methods.
"""

using JLD2
using Plots

RSOL = "sols/rsol_eval.jld2"

# Error grid evaluations
"Search for max value in grid with NaNs, print info and return val with index."
function errorg_info(g, ee, MM)
    mxvl, emxi, Mmxi = -Inf, -Inf, -Inf
    r, c = size(g)
    for i in 1:r
        for j in 1:c
            v = g[i, j]
            if !isnan(v) && mxvl < v
                mxvl = v
                emxi, Mmxi = i, j
            end
        end
    end
    
    nnan = sum(isnan.(g))
    if nnan > 0 println("NaN values = ", nnan) end
    println("Max error = ", round(mxvl, sigdigits=3),
            #" at e = ", round(ee[emxi], digits=3),
            #", M = ", round(MM[Mmxi], digits=3), "\n" )
            " at e = ", ee[emxi],
            ", M = ", MM[Mmxi], "\n" )
    return mxvl, (emxi, Mmxi)
end

# Plotting
"Plot a cubic function in the specified range"
function plot_cubic(f::Function, x0=-5, x1=5; v=0)
    x = range(x0, x1, 151)
    y = f.(x)
    p = plot(x, y, legend=false)
    hline!(p, [0])
    vline!(p, [v])
    display(p)
end

"Plot grid values in 3d and maximum's place cuts."
function plot_grid(g, ee, MM, zlabel="", titl="")
    emxi, Mmxi = errorg_info(g, ee, MM)[2]
    p1 = surface(MM, ee, g,
        xlabel="M", ylabel="e", zlabel="ϵ",
        title=titl)
    p2 = contourf(MM, ee, g,
        xlabel="M", ylabel="e", zlabel="ϵ",
        title=titl)
    p3 = plot(MM, g[emxi, :],
        xlabel="M", ylabel="ϵ",
        label=string("e = ",round(ee[emxi], digits=3)), legend=:topright,
        title = titl)
    p4 = plot(ee, g[:, Mmxi],
        xlabel="e", ylabel="ϵ",
        label=string("M = ",round(MM[Mmxi], digits=3)), legend=:topleft,
        title = titl)
    return p1, p3, p4
    #return p1, p2, p3, p4
end

# Helpers
"Solve cubic equation by Newton iteration."
function csol(e, M, f::Function, df::Function)
    x = ekepl(e, M)
    for i in 1:100
        x0 = x
        x = x0 - f(x0) / df(x0)
        if abs(x - x0) < 1e-14
            return x0
        end
    end
    println("No cubic solution found")
    return NaN
end

# Lambda evaluation
"Evaluate λs by plotting and printing error values on loaded evaluation grid."
function eval_lambda(lambdas; plot=false, error_f=errorf, vstarter_f=vstarter)
    sol_g, ee, MM = values(load(RSOL))
    starter_g = Array{Float64, 2}(undef, size(sol_g))
    error_g = Array{Float64, 2}(undef, size(sol_g))

    apply_vsolver!(starter_g, vstarter_f, ee, MM, lambdas)
    apply_grids!(error_g, error_f, starter_g, sol_g)
    epsilons = epsilon_max(starter_g, ee, MM)

    println("λs = ", round.(lambdas, sigdigits=4))
    #println("ϵs = ", round.(epsilons, digits=4))
    plts = plot_grid(error_g, ee, MM, "ϵ")
    if plot display.(plts) end
    return plts
end

"Evaluate iteration method by plotting and printing error values on eval grid."
function eval_solver(sol_f::Function; plt=false, error_f=errorf)
    sol_g, ee, MM = values(load(RSOL))
    e_g, M_g = fill_grid(ee, MM)
    iter_g = Array{Float64, 2}(undef, size(sol_g))
    error_g = Array{Float64, 2}(undef, size(sol_g))

    apply_solver!(iter_g, sol_f, e_g, M_g)
    apply_grids!(error_g, error_f, iter_g, sol_g)

    println("Method = '", string(sol_f), "'")
    plts = plot_grid(error_g, ee, MM, "ϵ")
    if plot display.(plts) end
    return plts
end

"Evaluate variable iteration method by plotting and printing error values."
function eval_vsolver(sol_f::Function, lambdas;
        sol=NaN, plt=false, error_f=errorf)
    if length(sol) != 3
        sol_g, ee, MM = values(load(RSOL))
    else
        sol_g, ee, MM = sol
    end
    e_g, M_g = fill_grid(ee, MM)
    replace!(e_g, 0. => NaN)
    replace!(M_g, 0. => NaN)
    iter_g = Array{Float64, 2}(undef, size(sol_g))
    error_g = Array{Float64, 2}(undef, size(sol_g))

    apply_vsolver!(iter_g, sol_f, e_g, M_g, lambdas)
    apply_grids!(error_g, error_f, iter_g, sol_g)

    println("Method = '", string(sol_f), "'")
    println("λs = ", round.(lambdas, sigdigits=4))
    plts = plot_grid(error_g, ee, MM, "ϵ")
    if plt display.(plts) end
    return iter_g, error_g
end

"Evaluate variable iteration method by plotting and printing error values."
function eval_vsolver_bf(sol_f::Function, lambdas;
        sol=NaN, plt=false, error_f=errorf)
    if length(sol) != 3
        sol_g, ee, MM = values(load(RSOL))
    else
        sol_g, ee, MM = sol
    end
    e_g, M_g = fill_grid(ee, MM)
    e_g = BigFloat.(e_g)
    M_g = BigFloat.(M_g)
    lambdas = BigFloat.(lambdas)
    replace!(e_g, 0. => NaN)
    iter_g = Array{Float64, 2}(undef, size(sol_g))
    error_g = Array{Float64, 2}(undef, size(sol_g))

    apply_vsolver!(iter_g, sol_f, e_g, M_g, lambdas)
    apply_grids!(error_g, error_f, iter_g, sol_g)

    println("Method = '", string(sol_f), "'")
    println("λs = ", round.(lambdas, sigdigits=4))
    plts = plot_grid(error_g, ee, MM, "ϵ")
    if plt display.(plts) end
    return iter_g, error_g
end

"Plot error in corner for e = 1."
function plot_e1(lambdas)
    sol_corner = values(load("sols/rsol_corner.jld2"))
    E_g, ee, MM = sol_corner
    iter_g = similar(E_g)
    error_g = similar(E_g)
    apply_vsolver!(iter_g, vkepler2, ee, MM, l2)
    apply_grids!(error_g, errorf, iter_g, E_g)

    x = MM
    y = error_g[end, :]

    xt = 10. .^ Array(-20:-2)
    yt = 10. .^ Array(-18:10)
    replace!(y, 0=>NaN)
    replace!(x, 0=>1e-10)
    p = scatter(x, y,
        xlabel="M", ylabel="ϵ", label="e = 1",
        xscale=:log10, xticks=xt, yscale=:log10, yticks=yt)
    display(p)
    return p
end

"Plot error in corner for M = 0."
function plot_M0(lambdas)
    sol_corner = values(load("sols/rsol_corner.jld2"))
    E_g, ee, MM = sol_corner
    iter_g = similar(E_g)
    error_g = similar(E_g)
    apply_vsolver!(iter_g, vkepler2, ee, MM, l2)
    apply_grids!(error_g, errorf, iter_g, E_g)

    x = 1 .- ee
    y = error_g[:, 2]

    xt = 10. .^ Array(-20:-2)
    yt = 10. .^ Array(-26:10)
    replace!(y, 0=>NaN)
    replace!(x, 0=>1e-10)
    p = scatter(x, y,
        xlabel="1 - e --> from e = 1 to e = .99", label="M = 0", ylabel="ϵ",
        xscale=:log10, xticks=xt, yscale=:log10, yticks=yt)
    display(p)
    return p
end
