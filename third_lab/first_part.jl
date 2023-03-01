import Pkg; Pkg.add("Plots"); Pkg.add("InverseFunctions")

using Plots
using InverseFunctions

# Constants 

length_of_range = 100
I_range = range(0, 50, length_of_range) 

# Function

function J_s(I_s, s, gamma)
    1 + gamma * (1 + I_s) / (1 + s * I_s) + I_s
end    

# Plots

plot!(J_s.(I_range, 2, 2), I_range, label = "Is, s = 2, gamma = 2")
plot!(J_s.(I_range, 1, 2), I_range, label = "Is, s = 1, gamma = 2")
plot!(J_s.(I_range, 2, 0), I_range, label = "Is, s = 2, gamma = 0")
xlabel!("J")
ylabel!("I")
xlims!(0, 50)
ylims!(0, 50)
savefig("Зависимость интенсивности от накачки")