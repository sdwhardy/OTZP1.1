using JuMP, Ipopt
include("../wind/struct.jl")
include("../cost/struct.jl")#
include("../eqp/struct.jl")#
include("../optimization/struct.jl")#
include("../layout/struct.jl")#

model = Model(with_optimizer(Ipopt.Optimizer))
@variable(model, x)
@variable(model, y)
@variable(model, lamda)

xy1=xy()
xy1.x=30
xy1.y=25
xy2=xy()
xy2.x=25
xy2.y=35
xyPCC=xy()
xyPCC.x=0
xyPCC.x=0

xys=Array{xy,1}()
#for owpp in ocean.owpps
for owpp in [xy1,xy2,xyPCC]
    x_y=xy()
    x_y.x=owpp.x
    x_y.y=owpp.y
    push!(xys,deepcopy(x_y))
end
f(x,y) = sum((x-xy.x[i])/(sqrt((xy.x-x[i])^2+(xy.y-y[i])^2)+epsi) for (i,xy) in enumerate(xys))
register(model, :f, length(xys), f; autodiff = true)
@NLconstraint(model, f(x,y)==2*x*lamda)

g(x,y) = sum((y[i]-xy.y)/(sqrt((xy.x-x[i])^2+(xy.y-y[i])^2)+epsi) for (i,xy) in enumerate(xys))
register(model, :g, length(xys), f; autodiff = true)
@NLconstraint(model, g(x,y)==2*y*lamda)

@constraint(model, x >= 0)
@constraint(model, y >= 0)
@NLconstraint(model, (x-xy1.x)^2+(y-xy1.y)^2 <= 7.365^2)
@NLconstraint(model, (x-xy2.x)^2+(y-xy2.y)^2 <= 7.365^2)

optimize!(model)
getvalue.((x, y, lamda))



m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m, x)
@variable(m, y)
@variable(m, lamda)
eps=1e-9

x1=30
y1=25
x2=25
y2=35
xPCC=0
yPCC=0
#(x-xPCC)/(sqrt((xPCC-x)^2+(yPCC-y)^2)+eps)+(x-x1)/(sqrt((x1-x)^2+(y1-y)^2)+eps)+(x-x2)/(sqrt((x2-x)^2+(y2-y)^2)+eps)
register(m, :f, 2, f; autodiff = true)
@NLconstraint(m, sum()==2*x*lamda)
@NLconstraint(m, (y-yPCC)/(sqrt((xPCC-x)^2+(yPCC-y)^2)+eps)+(y-y1)/(sqrt((x1-x)^2+(y1-y)^2)+eps)+(y-y2)/(sqrt((x2-x)^2+(y2-y)^2)+eps)==2*y*lamda)
@constraint(m, x >= 0)
@constraint(m, y >= 0)

@NLconstraint(m, (x-x1)^2+(y-y1)^2 <= 7.365^2)
@NLconstraint(m, (x-x2)^2+(y-y2)^2 <= 7.365^2)
optimize!(m)
JuMP.value.((x, y, lamda))

################################# 2 owpp together
using JuMP, Ipopt
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
xys=[xy0,xy1,PCC]
@expression((x-xPCC)/(sqrt((xPCC-x)^2+(yPCC-y)^2)+eps)
    +(x-x0)/(sqrt((x0-x)^2+(y0-y)^2)+eps)
    +(x-x1)/(sqrt((x1-x)^2+(y1-y)^2)+eps))

JuMP.register(m,:f, 2, f, autodiff=true)
R=sqrt((x0-xPCC)^2+(y0-yPCC)^2)-sqrt((250/6)/pi)
@NLconstraint(m, test==2*x*lamda)
@NLconstraint(m, (y-yPCC)/(sqrt((xPCC-x)^2+(yPCC-y)^2)+eps)+(y-y0)/(sqrt((x0-x)^2+(y0-y)^2)+eps)+(y-y1)/(sqrt((x1-x)^2+(y1-y)^2)+eps)==2*y*lamda)
@constraint(m, x >= 0)
@constraint(m, y >= 0)

@NLconstraint(m, x^2+y^2 == R^2)
optimize!(m)
JuMP.value.((x, y, lamda))

################################# 3 owpp together
using JuMP, Ipopt
m = Model(solver=IpoptSolver())
@variable(m, x)
@variable(m, y)
@variable(m, lamda)
eps=1e-9

#farthest
x2=25#C
y2=35
x1=30#B
y1=25


#closest
#x0=30#B
#y0=25
x0=10#A
y0=20
xPCC=0
yPCC=0
R=sqrt((x0-xPCC)^2+(y0-yPCC)^2)-sqrt((250/6)/pi)
@NLconstraint(m, (x-xPCC)/(sqrt((xPCC-x)^2+(yPCC-y)^2)+eps)+(x-x0)/(sqrt((x0-x)^2+(y0-y)^2)+eps)+(x-x1)/(sqrt((x1-x)^2+(y1-y)^2)+eps)+(x-x2)/(sqrt((x2-x)^2+(y2-y)^2)+eps)==2*x*lamda)
@NLconstraint(m, (y-yPCC)/(sqrt((xPCC-x)^2+(yPCC-y)^2)+eps)+(y-y0)/(sqrt((x0-x)^2+(y0-y)^2)+eps)+(y-y1)/(sqrt((x1-x)^2+(y1-y)^2)+eps)+(y-y2)/(sqrt((x2-x)^2+(y2-y)^2)+eps)==2*y*lamda)
@constraint(m, x >= 0)
@constraint(m, y >= 0)

@NLconstraint(m, x^2+y^2 == R^2)
solve(m)
getvalue.((x, y, lamda))
############################################################

using JuMP, Ipopt
m = Model(with_optimizer(Ipopt.Optimizer))
@variable(m, x)
@variable(m, y)
@variable(m, lamda)
eps=1e-9

#farthest
x2=25#C
y2=35
x1=30#B
y1=25


#closest
#x0=30#B
#y0=25
x0=10#A
y0=20
xPCC=0
yPCC=0
R=sqrt((x0-xPCC)^2+(y0-yPCC)^2)-sqrt((250/6)/pi)
@NLobjective(m,Min,sqrt((x0-x)^2+(y0-y)^2)+sqrt((x1-x)^2+(y1-y)^2)+sqrt((x2-x)^2+(y2-y)^2)+sqrt((xPCC-x)^2+(yPCC-y)^2))
@constraint(m, x >= 0)
@constraint(m, y >= 0)

@NLconstraint(m, x^2+y^2 == R^2)
optimize!(m)
getvalue.((x, y, lamda))
