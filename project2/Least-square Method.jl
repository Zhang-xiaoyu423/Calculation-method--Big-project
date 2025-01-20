############ Least square Method ##########################################################
# author: Zhang Xiaoyu
# 2024/11/3

# 拟合函数Fx
# A ：拟合多项式的系数矩阵，n+1行1列
# n : 拟合次数
# x : 自变量
# F_x : 因变量
function F(A, n, x)                                                         
    F_x = 0
    for i in 1:n+1
        F_x += A[i] * (x ^ (i-1))
    end
    return F_x
end

# 最小二乘法Least_square_Method
# x : 横坐标
# y : 纵坐标
# n : 拟合次数
function Least_square_Method(x, y, n)
    scatter(x, y, color=:red, label="Point", xlabel="x", ylabel="y")    # 绘制拟合点
    #-----------------------------最小二乘拟合-------------------------------#
    m = length(y)
    X = zeros(m, n+1)
    Y = zeros(m, 1)
    Y = y
    A = zeros(n+1, 1)
    Fx = zeros(1, m)
    err_Ls = 0.0
    for i in 1:m
        for j in 1:n+1
            X[i, j] = x[i] ^ (j-1)
        end
    end
    A = (X'*X)\(X'*Y)
    #-----------------------------拟合函数Fx--------------------------------#
    for i in 1:m
        Fx[i] = sum(A[j] * x[i]^(j-1) for j in 1:n+1)
        err_Ls += (Fx[i] - y[i]) ^ 2
    end
    #-----------------------------绘制函数Fx--------------------------------#
    F_x = [F(A, n, xi) for xi in x]
    plot!(x, F_x, label="Fx", color=:green, lw=2)
    # savefig("Least_square_Method.png")
end
    