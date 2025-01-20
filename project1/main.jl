############ Solving large-scale sparse equations ######################################################
# author: Zhang Xiaoyu
# 2024/10/21

using LinearAlgebra
using Plots

include("Conjugate Gradient Method.jl")
include("Gradient Descent Method.jl")
# include("Quasi Newton Method.jl")

############ Initial amount setting ####################################################################
n = 400                             # n is the matrix order
A = zeros(n, n)                     # Create a zero matrix A 
for i in 1:n - 1                    
    A[i, i] = -2                    # the primary diagonal element being -2
    A[i + 1, i] = 1                 # Set all elements around the main diagonal element to 1
    A[i, i + 1] = 1                 
end
A[n, n] = -2                        # the primary diagonal element being -2
b = zeros(n)                        # Create a zero vector b
b[1] = -1                           # the first and last elements as 1 and the rest as 0
b[end] = -1
ϵ = 1e-6                            # epsilon is the precision

Conjugate_Gradient_Method(A, n, b, ϵ)
Gradient_Descent_Method(A, n, b, ϵ)
# Quasi_Newton_Method(A, n, b, ϵ)

savefig("n400.png")
