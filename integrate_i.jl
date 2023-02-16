import Pkg; Pkg.add("Plots"); Pkg.add("DifferentialEquations")

using DifferentialEquations
using Plots


# Constants 

t_p = 8 * 10 ^ -12
t_c = 10 ^ -9
t_begin = 0
t_end = 1500
t_span = (t_begin, t_end)
t_peak = 1.4862020543941983
I_init = 3
N_init = 1

# Function

function ode_fn(du, u, p, t)
    I, N = u
    du[1] = (N - 1) * I
    du[2] = t_p / t_c * (1.5 - N * (1 + I))
end

prob = ODEProblem(ode_fn, [I_init, N_init], t_span)
num_sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

I_num_sol = [u[1] for u in num_sol.u]

# Plots 

plot!(num_sol.t, I_num_sol, label = "I(t)")
xlabel!("t, c")
ylabel!("I")
savefig("Нормированная интенсивность")

# Calculation

relax_freq = 1 / (num_sol.t[findfirst(isequal(t_peak), I_num_sol)]) / t_p
print(relax_freq)