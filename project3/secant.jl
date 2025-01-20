############ The Secant method ####################################################
# author: Zhang Xiaoyu
# 2024/11/22

include("func.jl")

# 弦割法 Secant 
# x_constraint  ：近似值
# ϵ1            : 函数边界
# A             ：方程的系数矩阵，n+1行1列
# n             : 最高阶次
function Secant(x_constraint, x_start, ϵ1, A, n)                                                       
    x_secant = 0
    k = 1
    while true
        x_secant = x_constraint - (F(A, n, x_constraint)*(x_constraint - x_start))/(F(A, n, x_constraint) - F(A, n, x_start)) 
        if abs(x_secant - x_constraint) < ϵ1
            return x_constraint
            break
        end
        x_constraint = x_secant  
        scatter!([k], [x_constraint], label=false, color=:blue)    
        k += 1
    end
    
end

