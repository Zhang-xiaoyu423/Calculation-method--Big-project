############ Finite Element Analysis ####################################################
# author: Zhang Xiaoyu
# 2024/12/02

function FE(nelx, nely, x, penal, KE)
    # Input  ：水平单元数nelx, 竖直单元数nely, 设计变量x, 惩罚因子penal, 单元刚度矩阵KE
    # Output ：全局节点位移U
    # 总体刚度矩阵的稀疏矩阵
    K = spzeros(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1))
    # 力矩阵的稀疏矩阵
    F = spzeros(2*(nely+1)*(nelx+1), 1)
    # U清零，用来保存更新的全局节点位移
    U = zeros(2*(nely+1)*(nelx+1), 1)
    for elx = 1 : nelx
        for ely = 1 : nely
            # 计算单元左上角、右上角节点编号
            n1 = (nely+1)*(elx-1) + ely
            n2 = (nely+1)* elx + ely
            # 同上主程序，计算单元4个节点8个自由度
            edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2]
            # 将单元刚度矩阵KE组装成 总体刚度矩阵K
            K[edof,edof] = K[edof,edof] + x[ely,elx]^penal*KE
        end
    end
    # 施加载荷，本算例应用了一个在右下角的垂直单元力
    F[2*(nelx+1)*(nely+1), 1] = -1
    # 施加约束，消除线性方程中的固定自由度来实现支承结构，本算例左边第一列固定
    fixeddofs = collect(1 : 2*(nely+1))
    # 剩下的不加约束的节点自由度，setdiff()从..中除去..
    alldofs   = collect(1 : 2*(nely+1)*(nelx+1))
    freedofs  = setdiff(alldofs, fixeddofs)
    # 求解线性方程组，得到各节点自由度的位移值储存在U中
    U[freedofs, :] = K[freedofs, freedofs] \ Matrix(F[freedofs, :])
    # 受约束节点固定自由度位移值为0
    U[fixeddofs, :] .= 0
    return U
end

function lk(E, μ)
    k=[ 1/2-μ/6,    1/8+μ/8,    -1/4-μ/12,  -1/8+3*μ/8,
        -1/4+μ/12,  -1/8-μ/8,   μ/6,        1/8-3*μ/8   ]
    KE = E/(1-μ^2).*[   k[1] k[2] k[3] k[4] k[5] k[6] k[7] k[8]
                        k[2] k[1] k[8] k[7] k[6] k[5] k[4] k[3]
                        k[3] k[8] k[1] k[6] k[7] k[4] k[5] k[2]
                        k[4] k[7] k[6] k[1] k[8] k[3] k[2] k[5]
                        k[5] k[6] k[7] k[8] k[1] k[2] k[3] k[4]
                        k[6] k[5] k[4] k[3] k[2] k[1] k[8] k[7]
                        k[7] k[4] k[5] k[2] k[3] k[8] k[1] k[6]
                        k[8] k[3] k[2] k[5] k[4] k[7] k[6] k[1] ]
    return KE
end
