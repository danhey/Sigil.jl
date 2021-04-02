#=
Miscellaneous physics
=#

# provides c0, G, m_H, k_B, sigma_sb, a_rad, Msun, Rsun, Lsun
# include("constants.jl")

function mu_from_composition(X, Y)
    #=
    Calculates the dimensionless mean molecular weight for a fully ionized gas.
    
    Parameters
    ----------
    X : hydrogen mass fraction
    Y : helium mass fraction
    
    Returns
    -------
    mu : mean molecular weight
    =#
    return 2 / (1 + 3 * X + 0.5 * Y)
end # mu_from_composition

function beta(Pgas, T)
    #=
    Calculates the ratio between gas pressure and total pressure.
    Parameters
    ----------
    Pgas : gas pressure in dyne cm^-2
    T : temperature in K
    Returns
    -------
    beta : Pgas / P
    =#
    Prad = a_rad * T .^ 4 / 3
    return Pgas / (Pgas + Prad)
end # beta

function rho_ideal(Pgas, T, mu)
    #=
    Calculates the density of an ideal fully ionized gas.
    Parameters
    ----------
    Pgas : gas pressure in dyne cm^-2
    T : temperature in K
    mu : dimensionless mean molecular weight
    Returns
    -------
    rho : density in g cm^-3
    =#
    Prad = a_rad * T .^ 4 / 3
    P = Pgas .+ Prad
    rho = P .* m_H .* mu ./ (k_B * T)
    return rho
end # rho_ideal

function acceleration(m, r)
    #=
    Calculates the local acceleration due to gravity.
    =#
    return G * m ./ (r .^ 2)
    
end # acceleration

function scale_height(m, r, P, T, mu)
    #=
    Calculates the pressure scale height.
    =#
    g = acceleration(m, r)
    rho = rho_ideal(P, T, mu)
    return P ./ (rho .* g)
end # scale_height
    
function T_surface(R, L_surface)
    #=
    Calculates the effective temperature, given the surface luminosity and the 
    stellar radius.
    =#
    return (L_surface ./ (4 * pi * R .^ 2 * sigma_sb)) .^ (1 / 4)
end # T_surface
    
function P_surface(M, R, κ)
    #=
    Calculates the surface pressure, given the total mass, radius, and surface
    opacity.
    =#
    g = acceleration(M, R)
    return 2 * g ./ (3 * κ)
end # P_surface

function ∇_ad()
    return 0.4
end # nabla_ad

function ∇_rad(m, l, P, T, κ)
    #=
    Gradient (d logP / d logT) in fully radiative case.
    =#
    return 3 * κ .* l .* P ./ (16 * pi * a_rad * c0 * G * m .* T .^ 4)
end # nabla_rad

function epsilon_pp(ρ, T, X)
    #=
    Proton-proton chain luminosity density
    =#
    T9 = T / 1e9
    if T / 1e7 >= abs(2 - 0.1)
        psi = 2.0
    else
        psi = 1.5
    end
    f11 = 1 # Weak screening approx
    g11 = 1. + 3.82 .* T9 + 1.151 .* T9.^2 + 0.144 .* T9.^3 - 0.0114.*T9.^4
    ϵ = (2.57e4) .* psi .* f11 .* g11 .* ρ .* (X ^ 2) .* (T9 .^ (-2/3)) .* exp(-3.381 ./ (T9 .^ (1/3)))
    return ϵ
end

function epsilon_cno(ρ, T, X, Y)
    #=
    CNO luminosity density
    =#
    T9 = T / 1e9
    X_CNO = 0.7 * (1 - X - Y)
    g141 = 1 - 2 .* T9 + 3.41 .* T9 .^ 2 - 2.43 .* T9 .^ 3
    ϵ = 8.24e25 .* g141 .* X_CNO .*  X .* ρ .* (T9 .^ (-2 / 3)) .* exp(-15.231 .* T9 .^ (-1 / 3) - (T9 / 0.8) .^ 2)
    return ϵ
end