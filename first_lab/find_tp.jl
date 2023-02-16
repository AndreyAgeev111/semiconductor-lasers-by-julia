import Pkg; Pkg.add("Plots"); Pkg.add("DifferentialEquations")

using DifferentialEquations
using Plots

# Constants 

J_c = [13.05, 13.3, 13.55, 13.8, 14.05, 14.3, 14.55, 14.8, 15.05, 15.3]
J_thr = 13
J = J_c./ J_thr
J_thr_step = 1.5
t_p = 8 * 10 ^ -12

F_c = [76.14, 257.76, 356.48, 433.26, 498.34, 555.85, 607.94, 655.9, 700.58, 742.58]
F = F_c.* 10 ^ 6

# Function

function frequency(j, t)
    return 1 / (4 * 10 ^ -9 * pi) * sqrt(4 * (10 ^ -9) * (j - 1) / t - (j ^ 2))
end

# Plots 

F_exp = frequency.(J, t_p)

plot(J, [F F_exp], label = ["F выданного варианта" "F после вычисленного tp"])
xlabel!("J, мА")
ylabel!("F, Гц")
savefig("Релаксационная частота")

# Calculation

relax_freq = frequency(J_thr_step, t_p)
println(relax_freq)