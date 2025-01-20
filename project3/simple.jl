############ The Simple method ####################################################
# author: Zhang Xiaoyu
# 2024/11/22

include("func.jl")

# 简单迭代法 Simple 
# x_constraint  ：近似值
# ϵ1            : 函数边界
# A             ：方程的系数矩阵，n+1行1列
# n             : 最高阶次
function Simple(x_constraint, ϵ1, A, n)  
    A_simple = -A[1 : n]                                                       
    x_simple = 0
    k = 1
    while true
        x_simple = (F(A_simple, n-1, x_constraint))^(1/n) 
        if abs(x_simple - x_constraint) < ϵ1
            return x_constraint
            break
        end
        x_constraint = x_simple   
        scatter!([k], [x_constraint], label=false, color=:red, xlabel="k", ylabel="x")
        k += 1    
    end   
end

