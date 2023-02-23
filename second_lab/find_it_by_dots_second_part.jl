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
g = 1.13
B = 500
I_init = 3
P_init = 0
N_init = 1
J_thr = ((g + 1) * (2 * B * (g - 1) + 2 * g)) / (B * (2 * g) * (g - 1))
J = J_thr * 1.5 * 0.75
delta_A = J * 0.1
f_a = 1.2 * 10 ^ 9
f_a_norm = f_a * t_p
length_of_range = 25
f_array = range(f_a_norm * 0.05, f_a_norm * 1.25, length_of_range)
f_max = zeros(length_of_range)
f_min = zeros(length_of_range)

# Function

function solve_freq(f)

    function ode_fn(du, u, p, t)
        I, P, N = u
        du[1] = (g * (2 * P - 1) - 1) * I
        du[2] = t_p / t_d * (B * N * (1 - P) - P - I * (2 * P - 1))
        du[3] = t_p / t_w * (J + delta_A * sin(2 * pi * f * t) - N - 2 * (B * N * (1 - P)))
    end
    
    prob = ODEProblem(ode_fn, [I_init, P_init, N_init], t_span)
    num_sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    
    I_num_sol = [u[1] for u in num_sol.u]
    I_num_sol = I_num_sol[trunc(Int, size(I_num_sol)[1] * 5 / 6):size(I_num_sol)[1]]
    I_max = maximum(I_num_sol)
    I_min = minimum(I_num_sol)

    return I_max, I_min

end    

# Calculation

for (i, f) in enumerate(f_array)
    (I_max, I_min) = solve_freq(f)
    f_max[i] = I_max
    f_min[i] = I_min
end

# Plots 

scatter!(f_array./ t_p / 10 ^ 9, f_max, color = "blue", label = "I_max(f)")
plot!(f_array./ t_p / 10 ^ 9, f_max, color = "blue", label = "I_max(f)")
scatter!(f_array./ t_p / 10 ^ 9, f_min, color = "red", label = "I_min(f)")
plot!(f_array./ t_p / 10 ^ 9, f_min, color = "red", label = "I_min(f)")
xlabel!("f, ГГц")
ylabel!("I")
savefig("Максимальная и минимальная интенсивность в модели квантовых точек с изменнным J")