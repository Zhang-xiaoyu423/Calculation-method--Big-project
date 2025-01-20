############ Conjugate Gradient Method #################################################################
# author: Zhang Xiaoyu
# 2024/10/21

############ The Conjugate Gradient Method #############################################################
# 共轭梯度法
# A ：系数矩阵，n行n列
# n : 矩阵A的阶数
# b : 右端项
# ϵ : 误差限度
function Conjugate_Gradient_Method(A, n, b, ϵ)
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
        push!(err_cg, log10(norm(r_cg)))
        rsold_cg = rsnew_cg
    end
    # Draw the plot
    plot(err_cg, label="Conjugate Gradient", xlabel="Iteration", ylabel="Log10 of Error") 
end
