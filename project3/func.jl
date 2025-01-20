############ Function and derivative construction ####################################################
# author: Zhang Xiaoyu
# 2024/11/22

# 目标方程Fx
# A ：方程的系数矩阵，n+1行1列
# n : 最高阶次
# x : 自变量
# F_x : 因变量
function F(A, n, x)                                                         
    F_x = 0
    for i in 1:n+1
        F_x += A[i] * (x ^ (i-1))
    end
    return F_x
end

# A ：方程的系数矩阵，n+1行1列
# n : 最高阶次
# x : 自变量
# dF_x : 因变量
function dF(A, n, x)                                                         
    dF_x = 0
    A_aux = A[2 : n+1]
    for i in 1:n
        dF_x += A_aux[i] * i * (x ^ (i-1))
    end
    return dF_x
end