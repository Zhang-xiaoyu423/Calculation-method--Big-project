############ Optimality Criteria method ####################################################
# author: Zhang Xiaoyu
# 2024/12/02

function OC(nelx, nely, x, volfrac, dc, isHole)
    # Input：  水平单元数nelx, 竖直单元数nely, 设计变量x, 材料体积比volfrac, 目标函数灵敏度dc
    # Output： 更新后的设计变量xnew
    # 定义一个取值区间，二分法，得到满足体积约束的拉格朗日算子
    l1 = 0
    l2 = 100000
    # 正向最大位移
    move = 0.2
    xnew = zeros(nely, nelx)
    while (l2-l1 > 1e-4)
        # 二分法，取区间中点
        lmid = 0.5*(l2 + l1)
        # sqrt(-dc./lmid)对应公式中Be^eta（eta=1/2），eta阻尼系数是为了确保计算的收敛性
        safe1 = min.(x .+ move, x .* sqrt.(-dc ./ lmid))
        # rows, cols = size(safe1)
        # println(" 行rows.= ", rows, " 列cols.= ", cols)
        safe2 = min.(1, safe1)
        safe3 = max.(x .- move, safe2)
        xnew = max.(0.001, safe3)
        xnew[ findall( !=(0), isHole ) ] .= 0.001
        # sum(sum(xnew))是更新后的材料体积, volfrac*nelx*nely是优化目标，用它们的差值判断是否收敛
        if sum( sum( xnew ) ) - volfrac*nelx*nely > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end
    return xnew
end