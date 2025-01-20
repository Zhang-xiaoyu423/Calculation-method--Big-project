############ Sensitivity filtering method ####################################################
# author: Zhang Xiaoyu
# 2024/12/02

# Unrelated mesh sensitivity filtering
function check(nelx, nely, rmin, x, dc)
    # Input  ：水平单元数nelx, 竖直单元数nely, 敏度过滤半径rmin, 设计变量x, 总体结构敏度dc
    # Output ：过滤后的目标函数敏度dcn
    # dcn清零，用来保存更新的目标函数灵敏度
    dcn = zeros(nely, nelx);
    for i = 1 : nelx
        for j = 1 : nely
            sum = 0.0;
            # % 在过滤半径定义的范围内遍历
            for k = max(i - Int(floor(rmin)), 1) : min(i + Int(floor(rmin)), nelx)
                for l = max(j - Int(floor(rmin)), 1) : min(j + Int(floor(rmin)), nely)
                    fac = rmin - sqrt((i-k)^2 + (j-l)^2);
                    sum = sum + max(0, fac);
                    dcn[j, i] = dcn[j, i] + max(0, fac)*x[l, k]*dc[l, k];
                end
            end
            dcn[j, i] = dcn[j, i] / (x[j, i]*sum);
        end
    end
    return dcn
end