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
@variable(m, t12)
@variable(m, t13)
@variable(m, t21)
@variable(m, t22)
@variable(m, t23)
@variable(m, t31)
@variable(m, t32)
@variable(m, t33)
@variable(m, lamda1)
@variable(m, p2)
@variable(m, p3)
@variable(m, p4)

t11=1
w0=1
w1=1
w2=1
x0=1
x1=-1
x2=1
y0=-1
y1=-1
y2=1
z0=-1
z1=1
z2=1

W0=14
W1=12
W2=1
X0=11
X1=2
X2=1
Y0=7
Y1=4
Y2=1
Z0=12
Z1=11
Z2=1
ocean.nogos[1]

@constraint(m, t11*w0+t12*w1+t13*w2 == p1*W0)
@constraint(m, t21*w0+t22*w1+t23*w2 == p1*W1)
@constraint(m, t31*w0+t32*w1+t33*w2 == p1*W2)

@constraint(m, t11*x0+t12*x1+t13*x2 == p2*X0)
@constraint(m, t21*x0+t22*x1+t23*x2 == p2*X1)
@constraint(m, t31*x0+t32*x1+t33*x2 == p2*X2)

@constraint(m, t11*y0+t12*y1+t13*y2 == p3*Y0)
@constraint(m, t21*y0+t22*y1+t23*y2 == p3*Y1)
@constraint(m, t31*y0+t32*y1+t33*y2 == p3*Y2)

@constraint(m, t11*z0+t12*z1+t13*z2 == p4*Z0)
@constraint(m, t21*z0+t22*z1+t23*z2 == p4*Z1)
@constraint(m, t31*z0+t32*z1+t33*z2 == p4*Z2)


optimize!(m)
vals=JuMP.value.((t12, t13, t21, t22, t23, t31, t32, t33,p1,p2,p3,p4))
vals[6]
T=[1 vals[1] vals[2];vals[3] vals[4] vals[5];vals[6] vals[7] vals[8]]
P=[vals[9] 0 0 0;0 vals[10] 0 0;0 0 vals[11] 0;0 0 0 vals[12]]
W=[W0 W1 W2;X0 X1 X2;Y0 Y1 Y2;Z0 Z1 Z2]
w=[w0 x0 y0 z0;w1 x1 y1 z1;w2 x2 y2 z2]
r0=P*W
r1=T*w
Tinv=inv(T)
c=[]

T0=[6445 1027 -1442;-29123 4447 -2468;-11035 -721 2054]/6445
T0*[w0;w1;w2]
############################################################
using JuMP, ProxSDP, LinearAlgebra,Polyhedra, SetProg, CDDLib
#polygon constraints
factory = with_optimizer(ProxSDP.Optimizer)
m = Model(factory)
simplex=opt_makeHalfSpace(ocean.nogos[1])
cheby_center, cheby_radius = chebyshevcenter(simplex, factory)
interior_point = SetProg.InteriorPoint(cheby_center)
model = Model(factory)
@variable(model, S, Ellipsoid(point=interior_point))
@constraint(model, S ⊆ simplex)
@objective(model, Max, nth_root(volume(S)))
optimize!(model)

using Plots
plot(polyhedron(simplex), xticks = 4:2:30,yticks = 0:2:30, ratio=1)
plot!(value(S))
############################################################
@variable(m, x)
@variable(m, y)
#western
l=ocean.nogos[1].wbnd[1]
@constraint(m, -x+l.m_findx*y <= -l.b_findx)
#southern
l=ocean.nogos[1].sbnd[1]
@constraint(m, l.m_findy*x-y <= -l.b_findy)
#northern
l=ocean.nogos[1].nbnd[1]
@constraint(m, -l.m_findy*x+y <= l.b_findy)
#eastern
l=ocean.nogos[1].ebnd[1]
@constraint(m, x-l.m_findx*y <= l.b_findx)
#simplex = HalfSpace([0, -1], 0) ∩ HalfSpace([-1, -1], -1) ∩ HalfSpace([-1, 1], 3) ∩ HalfSpace([1, 0], 3) ∩ HalfSpace([1, 2], 9)



@variable(m, S,Ellipsoid(symmetric=true))
@variable(m, b)
@variable(m, c)
@variable(m, d)
@variable(m, e)
@variable(m, f)
B=[a b;c d]
d=[e;f]
@objective(m, Max, B)
ph=node()
pt=node()
pt.xy.x=3
pt.xy.y=3
ph.xy.x=1
ph.xy.y=4
l=lof_getStr8line(ph,pt)
