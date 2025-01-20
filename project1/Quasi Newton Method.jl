############ Quasi Newton Method #######################################################################
# author: Zhang Xiaoyu
# 2024/10/21

############ The Quasi Newton Method ###################################################################
############ The Gradient Descent Method ################################################################
# 拟牛顿法
# A ：系数矩阵，n行n列
# n : 矩阵A的阶数
# b : 右端项
# ϵ : 误差限度
function Quasi_Newton_Method(A, n, b, ϵ)
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
        push!(err_qn, log10(norm(r_qn)))
        rsold_qn = rsnew_qn
    end
    # Draw the plot
    plot!(err_qn, label="Quasi Newton")
end
