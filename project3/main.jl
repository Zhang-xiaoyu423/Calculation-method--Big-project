############ The iterative method ####################################################
# author: Zhang Xiaoyu
# 2024/11/22

using LinearAlgebra
using Plots

include("binary.jl")
include("simple.jl")
include("newton.jl")
include("secant.jl")
include("func.jl")

A = [-20; 7; -7; 1; 3; -5; 1]
n = 6

Left = -1
Right = 5
ϵ1 = 1E-8
ϵ2 = 0.1

p = plot()

x_constraint = Binary_constraint(Left, Right, ϵ1, ϵ2, A, n)
x_simple = Simple(x_constraint, ϵ1, A, n)
x_newton = Newton(x_constraint, ϵ1, A, n)
x_start = 4
x_secant = Secant(x_constraint, x_start, ϵ1, A, n)

display(p)
# savefig("Plots.png")
