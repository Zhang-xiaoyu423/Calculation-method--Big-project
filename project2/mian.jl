############ Solving Least Squares Fitting Problem ###############################################
# author: Zhang Xiaoyu
# 2024/11/3

using LinearAlgebra
using Plots

include("Least-square Method.jl")

x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
y = [5.1234, 5.3057, 5.5687, 5.9375, 6.4370, 7.0978, 7.9493, 9.0253, 10.3627]

Least_square_Method(x, y, 4)
