############ The Binary method ####################################################
# author: Zhang Xiaoyu
# 2024/11/22

include("func.jl")

# 二分法隔离 Binary_constraint
# Left  ：左边界
# Right : 右边界
# ϵ1    : 函数边界
# ϵ2    : 区间限度
# A     ：方程的系数矩阵，n+1行1列
# n     : 最高阶次
function Binary_constraint(Left, Right, ϵ1, ϵ2, A, n)                                                         
    a = Left
    b = Right
    while true
        x_k = (a + b)/2
        if abs(F(A, n, x_k)) < ϵ1 || abs(a - b) < ϵ2
            x = x_k
            return x
            break
        elseif F(A, n, a)*F(A, n, x_k) < 0
            b = x_k
        else a = x_k
        end    
    end
    
end

