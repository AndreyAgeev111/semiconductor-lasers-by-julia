import Pkg; Pkg.add("Plots"); Pkg.add("DifferentialEquations")

using DifferentialEquations
using Plots


# Constants 

t_p = 15 * 10 ^ -12
t_d = 0.7 * 10 ^ -9
t_w = 10 ^ -9
t_begin = 0
t_end = 1500
t_span = (t_begin, t_end)
t_peak_first = 1.7067164798813914
t_peak_second = 0.8370686916195039
g = 1.13
B = 500
I_init = 3
P_init = 0
N_init = 1
J_thr = ((g + 1) * (2 * B * (g - 1) + 2 * g)) / (B * (2 * g) * (g - 1))
J = J_thr * 1.5

# Function

function ode_fn(du, u, p, t)
    I, P, N = u
    du[1] = (g * (2 * P - 1) - 1) * I
    du[2] = t_p / t_d * (B * N * (1 - P) - P - I * (2 * P - 1))
    du[3] = t_p / t_w * (J - N - 2 * (B * N * (1 - P)))
end

prob = ODEProblem(ode_fn, [I_init, P_init, N_init], t_span)
num_sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

I_num_sol = [u[1] for u in num_sol.u]

# Plots 

plot!(num_sol.t, I_num_sol, label = "I(t)")
xlabel!("t")
ylabel!("I")
savefig("Нормированная интенсивность 2")

# Calculation
relax_freq = 1 / (num_sol.t[findfirst(isequal(t_peak_second), I_num_sol)] - num_sol.t[findfirst(isequal(t_peak_first), I_num_sol)]) / t_p
print(relax_freq)