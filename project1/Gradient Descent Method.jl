############ The Gradient Descent Method ################################################################
# author: Zhang Xiaoyu
# 2024/10/21

############ The Gradient Descent Method ################################################################
# 梯度下降法
# A ：系数矩阵，n行n列
# n : 矩阵A的阶数
# b : 右端项
# ϵ : 误差限度
function Gradient_Descent_Method(A, n, b, ϵ)
    x_gd = zeros(n, 1)                      # Initialize the solution vector
    α_gd = 0.1                              # α is the fixed value
    k_gd = 0
    err_gd = []

    while k_gd < n
        k_gd += 1
        gradient = b - A * x_gd             # Calculate the gradient
        x_gd -= α_gd * gradient             # Update the solution vector
        push!(err_gd, log10(norm(gradient)))
        if norm(gradient) < ϵ               # If the error is small enough, stop the iteration
            break        
        end
    end
    # Draw the plot
    plot!(err_gd, label="Gradient Descent")
end

