############ The SMIP method ####################################################
# author: Zhang Xiaoyu
# 2024/12/02

using LinearAlgebra
using SparseArrays
using CairoMakie
CairoMakie.activate!()

include("topopt.jl")

# nelx    : 水平方向上的离散单元数
# nely    : 竖直方向上的离散单元数
# volfrac : 容积率，材料体积与设计域体积之比，对应的工程问题就是"将结构减重到百分之多少"
# penal   : 惩罚因子，SIMP方法是在0-1离散模型中引入连续变量x、系数p及中间密度单元，从而将离
#           散型优化问题转换成连续型优化问题，并且令0≤x≤1，p为惩罚因子，通过设定p>1对中间密
#           度单元进行有限度的惩罚，尽量减少中间密度单元数目，使单元密度尽可能趋于0或1
#           合理选择惩罚因子的取值，可以消除多孔材料，从而得到理想的拓扑优化结果
# rmin    : 敏度过滤半径，防止出现棋盘格现象
# E       : 材料杨氏弹性模量
# μ       : 材料泊松比

nelx = 90
nely = 30
volfrac = 0.5
penal = 3.0
rmin = 1.5
E = 1
μ = 0.3

TopOpt(nelx, nely, volfrac, penal, rmin, E, μ)
