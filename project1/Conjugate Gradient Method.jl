############ Conjugate Gradient Method #################################################################
# author: Zhang Xiaoyu
# 2024/10/21

using LinearAlgebra
using Plots

############ Initial amount setting ####################################################################
n = 100                             # n is the matrix order
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

############ The Conjugate Gradient Method #############################################################
x_cg = zeros(n, 1)                      # Initialize the solution vector
r_cg = b - A * x_cg                     # Calculate the initial residual
d_cg = r_cg                             # Set the initial search direction
rsold_cg = dot(r_cg, r_cg)              # Calculate the squared norm of the initial residual
k_cg = 0                                # Calculate the Iteration
err_cg = []                             # err is the error vector matrix

while k_cg < n                          # 2 norm > the precision, or k == n
    k_cg += 1
    Ad_cg = A * d_cg
    α_cg = rsold_cg / dot(d_cg, Ad_cg)  # Simplified formula
    x_cg += α_cg * d_cg
    r_cg = b - A * x_cg

    if norm(r_cg) <= ϵ
        break
    end

    rsnew_cg = dot(r_cg, r_cg)
    β_cg = rsnew_cg / rsold_cg          # Calculate the step size，Simplified formula
    d_cg = r_cg + β_cg * d_cg           # Update the search direction
    #push!(err_cg, norm(r_cg))
    push!(err_cg, log10(norm(r_cg)))
    rsold_cg = rsnew_cg
end

#=
############ The Quasi Newton Method ###################################################################
x_qn = zeros(n)
r_qn = b - A * x_qn
d_qn = r_qn
rsold_qn = dot(r_qn, r_qn)
k_qn = 0
err_qn = []

while norm(r_qn) >= ϵ && k_qn < n       # 2 norm > the precision, or k == n
    k_qn += 1
    Ad_qn = A * d_qn
    α_qn = rsold_qn / dot(d_qn, Ad_qn)  # Simplified formula
    x_qn += α_qn * d_qn
    r_qn -= α_qn * Ad_qn
    rsnew_qn = dot(r_qn, r_qn)
    d_qn = r_qn + (rsnew_qn / dot(r_qn, Ad_qn)) * d_qn
    #push!(err_qn, norm(r_qn))
    push!(err_qn, log10(norm(r_qn)))
    rsold_qn = rsnew_qn
end
=#

############ The Gradient Descent Method ################################################################
x_gd = zeros(n, 1)                      # Initialize the solution vector
α_gd = 0.1                              # α is the fixed value
k_gd = 0
err_gd = []

while k_gd < n
    k_gd += 1
    gradient = b - A * x_gd             # Calculate the gradient
    x_gd -= α_gd * gradient             # Update the solution vector
    #push!(err_gd, norm(gradient))
    push!(err_gd, log10(norm(gradient)))
    # 如果误差足够小，停止迭代
    if norm(gradient) < ϵ
        break        
    end
end

############ Draw the plot ##############################################################################
plot(err_cg, label="Conjugate Gradient", xlabel="Iteration", ylabel="Log10 of Error")#, title="Error vs. Iteration for CG and Quasi-Newton Methods")
plot!(err_gd, label="Gradient Descent")
#savefig("n100.png")
#=
plot!(err_qn, label="Quasi-Newton")
xaxis!("Iteration")
yaxis!("Log10 of Error")

savefig("plot.png")
=#