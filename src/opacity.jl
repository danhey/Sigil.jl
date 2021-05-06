using Dierckx
using DelimitedFiles
using DataFrames
using CSV

modelS = DataFrame(CSV.File("data/modelS.dat"))
sort!(Sigil.modelS, ["m"])

function opal_test_spl()
    # function opal_test_spl()
    #     opal_all = myarray=(open(readdlm,"data/opal_all.txt"))
    #     i = 78
    #     table = opal_all[182+ (i-1)*69 + (i-1)*4:182+i*(69)+(i-1)*4,:]
    #     replace!(table, 9.999 => 0);
    #     # replace!(table, "" => 0);
    #     logT = table[:,1]
    table = myarray=(open(readdlm,"data/solar.txt"))
    replace!(table, "" => 0);
    logT = table[:,1]
    logR = Array(range(-8, 1., step=0.5))
    logkappa = table[:,2:end]

    spline = Spline2D(logT, logR, logkappa)
    return spline
end

function opacity(spl, rho, T)
    R = rho ./ (T / 1e6) .^ 3
    # @info(log10(T), log10(R))
    return 10 .^ spl(log10.(T), log10.(R))
end

function X_spl()
    spline = Spline1D(modelS["m"], modelS["X"])
    return spline
end

function Y_spl()
    spline = Spline1D(modelS["m"], modelS["Y"])
    return spline
end

function K_spl()
    spline = Spline1D(modelS["m"], modelS["kappa"])
    return spline
end

function rho_spl()
    spline = Spline1D(modelS["m"], modelS["rho"])
    return spline
end
