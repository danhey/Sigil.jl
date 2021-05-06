using DifferentialEquations
using Optim

const X_spline = X_spl()
const Y_spline = Y_spl()
const K_spline = K_spl()
# const rho_spline = rho_spl()

function init_guess(M, X, Y)
    #=
    Gives a reasonable guess for Pc, Tc, R, and L for a star of mass M
    =#
    yhalf = [Rsun * 10 ^ 0.062544715341505616, Lsun * 10 ^ -0.36630780815983022, 10 ^ 19.700356657315954, 10 ^ 7.2863105891013289]
    y2 = [Rsun * 10 ^ 0.36072469919806871, Lsun * 10 ^ 1.3021335874028637, 10 ^ 17.186879800538229, 10 ^ 7.3557279873577759]
    y5 = [Rsun * 10 ^ 0.59322648100945885, Lsun * 10 ^ 2.8743140385067893, 10 ^ 16.752833342638592, 10 ^ 7.4547058342916337]
    y4 = [Rsun * 10 ^ 0.52962615387538248, Lsun * 10 ^ 2.50275963393719, 10 ^ 16.855407680647915, 10 ^ 7.429616124439673]
    y10 = [Rsun * 10 ^ 0.68191099104456998, Lsun * 10 ^ 3.8447263619458569, 10 ^ 16.517451455410988, 10 ^ 7.5012272142698366]
    if M == Msun
        return [Rsun, Lsun, 2.4e17, 1.6e7]
    elseif M == 0.5 * Msun
        @info("init_guess: from MESA-Web")
        return yhalf
    elseif M == 2 * Msun
        @info("init_guess: from MESA-Web")
        return y2
    elseif (M >= 4 * Msun) & (M <= 5 * Msun)
        @info("init_guess: from MESA-Web")
        dm = M - 4 * Msun
        dydm = (y5 - y4) / Msun
        return y4 + dydm * dm
    elseif M == 10 * Msun
        @info("init_guess: from MESA-Web")
        return y10
    end
    mu = mu_from_composition(X, Y)

    # mass-radius relation
    if M / Msun > 1.3 # then CNO burning dominates
        z1 = 0.8
        z2 = 0.64
    else # then pp chain dominates
        z1 = 0.5
        z2 = 0.1
    end # if
    R = Rsun * (M / Msun) ^ z1 * (mu / 0.61) ^ z2
    # mass-luminosity relation
    L = Lsun * (M / Msun) ^ 3 * (mu / 0.61) ^ 4
    # Pc, Tc guesses, essentially dimensional analysis with fudges
    Pc = 2 * G * M ^ 2 / (pi * R ^ 4)
    Tc = 8 * 0.02 * mu * m_H * G * M / (k_B * R)
    # Tc = 1.57e7
    return [R, L, Pc, Tc]
end # init_guess

# function load1(m, Pc, Tc, X, Y)
#     #=
#     Subroutine for shootf.
#     Loads the (r, l, P, T) center (m = 0) values for starting the outward
#     integration.
#     =#
#     mu = mu_from_composition(X, Y)
#     # m = 0.01 * Msun # needs tuning
#     rhoc = rho_ideal(Pc, Tc, mu)
#     rc = (3 * m / (4 * pi * rhoc)) ^ (1 / 3)
#     lc = dldm(m, 0, 0, 0, Tc, rhoc, X, Y) * m
#     # @info("load1: center values = ", [rc, lc, Pc, Tc])
#     return [rc, lc, Pc, Tc]
# end # load1

function load2(M, Rs, Ls, X, Y; max_iterations=1000)
    #=
    Subroutine for shootf.
    Loads the (r, l, P, T) surface (m = M) values for starting the inward
    integration.
    =#
    mu = mu_from_composition(X, Y)
    converged = false
    rho_guess = 1e-8 # need a less arbitrary guess?
    Ts = try
        T_surface(Rs, Ls)
    catch
        @warn("load2: Rs = ", Rs)
        @warn("load2: Ls = ", Ls)
    end
    Ps = 0
    count = 1
    while ~converged & (count < 100)
        # kappa = opacity(logkappa_spl, rho_guess, Ts)
        kappa = modelS["kappa"][end] # Kappa at the surface..
        # kappa = K_spline(M)
        # kappa = K_spl(M)
        Ps = P_surface(M, Rs, kappa)
        rho = rho_ideal(Ps, Ts, mu)
        # rho = modelS["rho"][end]
        relerr = (rho - rho_guess) / rho
        if abs(relerr) < 1e-6 # need less arbitrary convergence condition
            converged = true
        else
            rho_guess = rho_guess + relerr * rho / 6
            count += 1
        end # if
    end # while
    if count >= max_iterations
        @warn("load2: Exceeded ", max_iterations, " iterations!")
    end # if
    @info("load2: surface values = ", [Rs, Ls, Ps, Ts])
    return [Rs, Ls, Ps, Ts]
end # load2

function profiles(m, M, Rs, Ls, Pc, Tc, X, Y, mf)
    #=
    Integrates outward from m to mf and inward from M to mf to construct
    r, l, P, and T profiles from the given initial conditions.
    =#

    # X, Y = X_spli ne(m), Y_spline(m)
    # yc = load1(m, Pc, Tc, X_spline(0), Y_spline(0))
    yc = [1e-5, 1e-5, Pc, Tc] # <-- here I fix the interior BC
    prob = ODEProblem(deriv!,yc,(0, mf), (X, Y))
    sol_out = solve(prob, Tsit5());
    sol_out_v = hcat(sol_out.u...)'

    ys = load2(M, Rs, Ls, X_spline(M), Y_spline(M))
    # If I fix the outer BC then the solution is unstable
    # the radius always wants to exceed ~1.5
    # ys = [Rs, Ls, 1, 1]
    # ys[1] = Rsun
    # @info(ys)
    # ys = [Rs, Ls, 1e-5, 1e-5]
    prob = ODEProblem(deriv!,ys,(M, mf), [X, Y])

    sol_in = solve(prob, Tsit5());
    sol_in_v = hcat(sol_in.u...)'
    return sol_out.t, sol_in.t, sol_out_v, sol_in_v
end

function profiles(star::Star, Rs, Ls, Pc, Tc; n=1000)
    return profiles(star.m, star.M, Rs, Ls, Pc, Tc, star.X, star.Y, star.mf)
end

function deriv!(du,u,p,m)
    r, l, P, T = u
    X, Y = X_spline(m), Y_spline(m)
    # @info(X, Y)
    mu = mu_from_composition(X, Y)
    rho = rho_ideal(P, T, mu)

    # kappa = try
    #     K_spline(m)
    # catch
    #     @info("mmmm")
    # kappa =opacity(logkappa_spl, rho, T)
    # rho = rho_spline(m) # this doesn't work.. units???
    kappa = K_spline(m)

    du[1] = drdm(m, r, l, P, T, rho, mu)
    du[2] = dldm(m, r, l, P, T, rho, X, Y)
    du[3] = dPdm(m, r, l, P, T)
    du[4] = dTdm(m, r, l, P, T, mu, kappa)
end

function score(m, M, Rs, Ls, Pc, Tc, X, Y, mf)
    #=
    Function to zero.
    =#
    ss = try
        sol_out_t, sol_in_t, sol_out_v, sol_in_v = profiles(m, M, Rs, Ls, Pc, Tc, X, Y, mf)
        (sol_out_v[end, :] - sol_in_v[end, :]) ./ sol_in_v[end, :]
    catch
        [Inf]
    end

    return (sum(ss.^2))
end # score

function shootf(m, M, X, Y; fixed_point=0.333)
    #=
    Shooting to a fixed point
    =#
    X, Y = X_spline(m), Y_spline(m)
    y0 = init_guess(M, X, Y)
    mf = fixed_point * M
    f(y) = score(m, M, y[1], y[2], y[3], y[4], X, Y, mf)

    res = optimize(f, y0,
                g_tol = 1e-15,
                iterations = 10000,
                allow_f_increases=true,
                store_trace = false,
                show_trace = false
    )
    return res
end # shootf

function shootf(star::Star)
    return shootf(star.m, star.M, star.X, star.Y)
end

function shootf!(star::Star)
    res = shootf(star)
    bc = Optim.minimizer(res)
    sol_out_t, sol_in_t, sol_out_v, sol_in_v = profiles(star, bc...)
    star.grid = vcat(sol_out_t, reverse(sol_in_t))
    star.r = vcat(sol_out_v[:,1], reverse(sol_in_v[:,1]))
    star.l = vcat(sol_out_v[:,2], reverse(sol_in_v[:,2]))
    star.P = vcat(sol_out_v[:,3], reverse(sol_in_v[:,3]))
    star.T = vcat(sol_out_v[:,4], reverse(sol_in_v[:,4]));
    return res
end
