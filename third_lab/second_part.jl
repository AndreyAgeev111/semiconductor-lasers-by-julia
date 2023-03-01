import Pkg; Pkg.add("Plots"); Pkg.add("DifferentialEquations")

using DifferentialEquations
using Plots


# Constants 

t_p = 8 * 10 ^ -12
t_c = 10 ^ -9
t_q = 0.8 * t_p
t_begin = 0
t_end = 1000
t_span = (t_begin, t_end)
gamma = 2
s = 2
I_init = 1
N_one_init = 1
N_two_init = 0
J_thr = 1 + gamma
J = J_thr * 3.9

# Function

function model_modified(J_mod, s_mod, t_q_mod)
    
    function ode_fn(du, u, p, t)
        I, N_one, N_two = u
        du[1] = (N_one + N_two - 1) * I
        du[2] = t_p / t_c * (J_mod - N_one * (1 + I))
        du[3] = t_p / t_q_mod * (-gamma - N_two * (1 + s_mod * I))
    end

    function findlocalmaxima(signal::Vector)
        indexs = Int[]
        if length(signal) > 1
            if signal[1] > signal[2]
                push!(indexs, 1)
            end
            for i = 2 : length(signal) - 1
                if signal[i - 1] < signal[i] > signal[i + 1]
                    push!(indexs, i)
                end
            end
            if signal[end] > signal[end - 1]
                push!(indexs, length(signal))
            end
        end
        indexs
    end

    prob = ODEProblem(ode_fn, [I_init, N_one_init, N_two_init], t_span)
    num_sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

    I_num_sol = [u[1] for u in num_sol.u]
    N_one_num_sol = [u[2] for u in num_sol.u]
    N_two_num_sol = [u[3] for u in num_sol.u]

    local_maximums = findlocalmaxima(I_num_sol) 
    t_first_maximum = num_sol.t[local_maximums[2]] * t_p * 10 ^ 9
    t_second_maximum = num_sol.t[local_maximums[3]] * t_p * 10 ^ 9

# Plots 

    plot!(num_sol.t.* (t_p * 10 ^ 9), I_num_sol, label = "I(t)")
    plot!(num_sol.t.* (t_p * 10 ^ 9), N_one_num_sol, label = "N_1(t)")
    plot!(num_sol.t.* (t_p * 10 ^ 9), N_two_num_sol, label = "N_2(t)")
    xlabel!("t, нс")
    ylabel!("I, N_1, N_2")
    savefig("Нормированная интенсивность и концентрация от исходных данных exp")

# Calculation

    I_cut = I_num_sol[local_maximums[2]:size(I_num_sol)[1]]
    t_cut = num_sol.t[local_maximums[2]:size(I_num_sol)[1]]

    for (i, I) in enumerate(I_cut)

        if (I_cut[1] / 2) / I > 1
            println((t_cut[i] - t_cut[1]) * 2 * t_p * 10 ^ 12)
            break
        end  

    end       

    relax_freq = 1 / (t_second_maximum - t_first_maximum)
    println(relax_freq)

end

# Proccessing 

J_modified = [0.75, 1.25, 1.5, 2].* J
s_modified = s / 2
t_q_modified = t_q * 2

model_modified(J_modified[4], s, t_q)