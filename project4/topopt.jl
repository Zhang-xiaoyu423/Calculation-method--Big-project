############ The topology optimization method ###########################################
# author: Zhang Xiaoyu
# 2024/12/02

include("problem.jl")
include("fea.jl")
include("filter.jl")
include("oc.jl")

function TopOpt(nelx, nely, volfrac, penal, rmin, E, μ)
    x, isHole = Problem(nelx, nely, volfrac)
    # loop储存迭代次数
    loop = 0
    # change储存每次迭代之后目标函数的改变值，用以判断是否收敛
    change = 1
    # 输出loop = 0图像（初始结构图象）
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())
    ax.title = "Topology Optimization"
    ax.xlabel = "The x axial grid"
    ax.ylabel = "The y axial grid"
    ax.xaxisposition = :top
    ax.yaxisposition = :left
    ax.subtitle = "Iteration = $loop"
    ax.subtitlesize = 15
    xlims!(ax, 0, 90)
    ylims!(ax, 30, 0)
    CairoMakie.heatmap!(ax, -x'; colormap=:grays)  # 使用 -x' 作为热力图数据，你可以根据需要调整
    display(fig)
    save("iteration_$loop.png", fig)
    # 当目标函数改变量<=0.01时说明收敛，结束迭代
    while change > 0.01
        loop = loop + 1
        # 将前一次的设计变量赋值给xold，x用来储存这一次的结果，之后还要比较它们以判断是否收敛
        xold = x
        # 计算单元刚度矩阵
        KE = lk(E, μ)
        # 每次迭代都进行一次有限元分析，计算结点位移，并储存在全局位移数组U中
        U = FE(nelx, nely, x, penal, KE)
        # c是用来储存目标函数的变量
        c = 0
        dc = zeros(nely, nelx)
        # 遍历设计域矩阵元素，从左上角第一个元素开始，一列一列
        for ely = 1:nely
            for elx = 1:nelx
                # 节点位移储存在U中，如果想获得节点位移必须先知道节点编号，进行索引
                # 节点编号可以根据当前单元在设计域中的位置算出
                # n1是左上角节点编号，n2是右上角节点编号
                n1 = (nely+1)*(elx-1) + ely
                n2 = (nely+1)*elx + ely
                # 局部位移数组Ue储存4个节点共8个自由度位移，每个节点分别有x、y两个自由度
                # 因为是矩形单元，所以根据n1、n2两个节点的编号可以推演出单元所有节点的自由度编号
                # 顺序是：[左上x；左上y；右上x；右上y；右下x；右下y；左下x；左下y]
                # 只适用于矩形单元划分网格;
                Ue = U[[2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2], 1]
                # SIMP模型，将设计变量x从离散型变成指函数型，指数就是惩罚因子：x(ely,elx)^penal
                # 计算总体结构的柔度值，这里目标函数是柔度最小，参见论文中公式（1）
                c = c + x[ely, elx]^penal * Ue' * KE * Ue
                # 计算总体结构的敏度值，实际上dc就是c对Xe的梯度，参见论文中公式（4）
                dc[ely, elx] = -penal*x[ely,elx]^(penal-1)*Ue'*KE*Ue
            end
        end
        # 无关网格敏度过滤
        dc = check(nelx, nely, rmin, x, dc);
        # 采用优化准则法（OC）求解当前模型，得出满足体积约束的结果，更新设计变量
        x = OC(nelx, nely, x, volfrac, dc, isHole)
        # 更新目标函数改变值
        change = maximum( max.( abs.(x - xold) ) );
        # 打印迭代信息: It.迭代次数，Obj.目标函数，Vol.材料体积比，ch.迭代改变量
        println(" It.: ", loop, ", Obj.= ", c, ", Vol.= ", sum(sum(x))/(nelx*nely), ", ch.= ", change)
        # 输出迭代图像
        ax.subtitle = "Iteration = $loop"
        CairoMakie.heatmap!(ax, -x'; colormap=:grays)  # 使用 -x' 作为热力图数据，你可以根据需要调整
        display(fig)
        save("iteration_$loop.png", fig)       
    end
end


