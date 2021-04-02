using Plots

function plot_star(star::Star)
    p1 = plot(star.r / Rsun, star.grid / Msun, ylabel="\$(M/M_{sun})\$", xlabel="\$(R/R_{sun})\$", color="blue")
    p2 = plot(star.r / Rsun, star.l / Lsun, xlabel="\$(R/R_{sun})\$", ylabel="\$(L/L_{sun})\$", color="blue")
    p3 = plot(star.r / Rsun, star.P, xlabel="\$(R/R_{sun})\$", ylabel="\$(dyne/cm^2)\$", color="blue")
    p4 = plot(star.r / Rsun, star.T, xlabel="\$(R/R_{sun})\$", ylabel="\$T\$", color="blue")
    plot(p1, p2, p3, p4,layout=(2,2), legend=false)
end