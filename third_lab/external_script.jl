import Pkg; Pkg.add("Plots")

using Plots

t_imp = [19.24, 19.1, 19.18, 19.34, 19.8]
f_imp = [3.69, 1.68, 2.22, 2.8, 6.95]

gamma = 2
J_thr = 1 + gamma
J = J_thr * 3.9
J_imp = [0.75, 1, 1.25, 1.5, 2].* J

plot(J_imp, t_imp, label = "t_imp(J)")
xlabel!("J")
ylabel!("t, пс")
savefig("Зависимость времени импульса от параметра накачки")

plot(J_imp, f_imp, label = "f_imp(J)")
xlabel!("J")
ylabel!("f, ГГц")
savefig("Зависимость частоты генерации от параметра накачки")

