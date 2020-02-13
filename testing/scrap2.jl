using JuMP, Ipopt
include("../wind/struct.jl")
include("../cost/struct.jl")#
include("../eqp/struct.jl")#
include("../optimization/struct.jl")#
include("../layout/struct.jl")#

m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m, x)
@variable(m, y)
@variable(m, lamda)
eps=1e-9

#farthest
#x1=25.05#C
#y1=25.05
xy1=xy()
xy1.x=30#B
xy1.y=25


#closest
#x0=30#B
#y0=25
xy0=xy()
xy0.x=10#A
xy0.y=20
PCC=xy()
PCC.x=0
PCC.y=0
xys=[PCC,xy0,xy1]

R=sqrt((xy0.x-PCC.x)^2+(xy0.y-PCC.y)^2)-sqrt((250/6)/pi)
@NLconstraint(m, sum((x-xys[i].x)/(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys))==2*x*lamda)



@NLconstraint(m, (y-PCC.y)/(sqrt((PCC.x-x)^2+(PCC.y-y)^2)+eps)+(y-xy0.y)/(sqrt((xy0.x-x)^2+(xy0.y-y)^2)+eps)+(y-xy1.y)/(sqrt((xy1.x-x)^2+(xy1.y-y)^2)+eps)==2*y*lamda)
@constraint(m, x >= 0)
@constraint(m, y >= 0)

@NLconstraint(m, x^2+y^2 == R^2)
optimize!(m)
JuMP.value.((x, y, lamda))
