
const logkappa_spl = opal_test_spl()

function drdm(m, r, l, P, T, rho, mu)
    #=
    Derivative of r with respect to m, from mass conservation.
    =#
    return 1 / (4 * pi * r .^ 2 .* rho)
end # drdm

function dPdm(m, r, l, P, T)
    #=
    Derivative of P with respect to m, from hydrostatic equilibrium.
    =#
    return -G * m / (4 * pi * r .^ 4)
end # dPdm

function dldm(m, r, l, P, T, rho, X, Y)
    #=
    Derivative of l with respect to m, from energy conservation.
    =#
    return epsilon_pp(rho, T, X) + epsilon_cno(rho, T, X, Y)
end # dldm

function dTdm(m, r, l, P, T, mu, kappa)
    #=
    Derivative of T with respect to m, from... thermodynamics.
    =#
    grad_rad = ∇_rad(m, l, P, T, kappa)
    grad_ad = ∇_ad()
    stable =  grad_rad < grad_ad
    if stable
        # @debug("dTdm: stable to convection")
        grad = grad_rad
    else
        # @debug("dTdm: unstable to convection")
        grad = grad_ad
    end # if
    factor = - G .* m .* T ./ (4 * pi .* r .^ 4 .* P)
     return factor .* grad
end # dTdm