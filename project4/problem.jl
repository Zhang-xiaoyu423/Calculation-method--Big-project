############ The design domain (problem) ##########################################
# author: Zhang Xiaoyu
# 2024/12/02

# x : 设计变量，给设计域内的每个单元一个初始相对密度，值为volfrac
function Problem(nelx, nely, volfrac)
    x = zeros(nely, nelx)
    isHole = zeros(nely, nelx)
    x .= volfrac;
    for ely = 1 : nely
        for elx = 1 : nelx
            if sqrt((ely-nely/2)^2 + (elx-nelx/3)^2) < nely/3
                isHole[ely, elx] = true
                x[ely, elx] = 0.001
            else
                isHole[ely,elx] = false
            end
        end
    end
    return x, isHole
end