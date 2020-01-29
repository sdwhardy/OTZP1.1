using Plots, JuMP, LinearAlgebra, Polyhedra, SetProg, Mosek, MosekTools, TypedPolynomials, MathOptInterface
#solvers
mosek_solver = with_optimizer(Mosek.Optimizer, QUIET=false)

#construct various test polyhedron
#unit circle centred at (0,0)
unit_sqr=HalfSpace(Float32[0.0, -1.0], 0.5) ∩ HalfSpace(Float32[-1.0, 0.0], 0.5) ∩ HalfSpace(Float32[0.0, 1.0], 0.5) ∩ HalfSpace(Float32[1.0, 0.0], 0.5)
#unit circle offset to (0.5,0.5)
off_set_unit_sqr=HalfSpace(Float32[0.0, -1.0], -0.0) ∩ HalfSpace(Float32[-1.0, 0.0], -0.0) ∩ HalfSpace(Float32[0.0, 1.0], 1.0) ∩ HalfSpace(Float32[1.0, 0.0], 1.0)
#1x2 rectangle centred at (0,0)
rec_1x2=HalfSpace(Float32[0.0, -1.0], 0.5) ∩ HalfSpace(Float32[-1.0, 0.0], 1.0) ∩ HalfSpace(Float32[0.0, 1.0], 0.5) ∩ HalfSpace(Float32[1.0, 0.0], 1.0)
#5 sided offset, rotated polygon
complex_shape=HalfSpace(Float32[0.5, -1.0], 1.0) ∩ HalfSpace(Float32[-0.5, -1.0], -0.0) ∩ HalfSpace(Float32[-1.0, 0.5], -0.0) ∩ HalfSpace(Float32[-0.33333334, 1.0], 0.8333333) ∩ HalfSpace(Float32[1.0, 1.0], 3.5)

#Select which shape to optimize
#poly_gone = unit_sqr
#poly_gone = off_set_unit_sqr
#poly_gone = rec_1x2
poly_gone = complex_shape

#Find largest area ellipse within polygone
Elips=findLargestEllipse(poly_gone,mosek_solver)
#plot polygone
plot(polyhedron(poly_gone), xticks = -1:0.25:3,yticks = -1:0.25:1.5, ratio=1)
plot!(value(Elips))

#Extras
#find A,B,C,D,E,F constants
A,B,C,D,E,F=findABCDEF(Elips)
#find x0 and y0
x0=getX0(A,B,C,D,E)
y0=getY0(A,B,C,D,E)
#find x radius and y radius
r_x,r_y=getRadiuss(A,B,C,D,E,F)
#determinant
Q=[A B/2;B/2 C]
det(Q)

qa=[A B D;B C E;D E F]


#################### functions ###################
function findLargestEllipse(poly_g,_solver)
    cntr, rad = chebyshevcenter(poly_g, _solver)
    inner_pnt = SetProg.InteriorPoint(cntr)
    m = Model(_solver)
    @variable(m, Srf, Ellipsoid(point=inner_pnt))
    @constraint(m, Srf ⊆ poly_g)
    @objective(m, Max, nth_root(volume(Srf)))
    optimize!(m)
    return Srf
end

function findABCDEF(S)
    A=convert(Float64,subs(value(S).set.p[4], value(S).set.x[1]=>1))
    B=convert(Float64,subs(value(S).set.p[5], value(S).set.x[1]=>1,value(S).set.x[2]=>1)/2)
    C=convert(Float64,subs(value(S).set.p[6], value(S).set.x[2]=>1))
    D=convert(Float64,subs(value(S).set.p[2], value(S).set.z=>1,value(S).set.x[1]=>1)/2)
    E=convert(Float64,subs(value(S).set.p[3], value(S).set.z=>1,value(S).set.x[2]=>1)/2)
    F=convert(Float64,subs(value(S).set.p[1], value(S).set.z=>1))
    return A,B,C,D,E,F
end

function getX0(A,B,C,D,E)
    x0=(2*C*D-B*E)/(B^2-4*A*C)
    return x0
end

function getY0(A,B,C,D,E)
    y0=(2*A*E-B*D)/(B^2-4*A*C)
    return y0
end

function getRadiuss(A,B,C,D,E,F)
    a=(-sqrt(complex(2*(A*E^2+C*D^2-B*D*E+(B^2-4*A*C)*F)*((A+C)+sqrt((A-C)^2+B^2)))))/(B^2-4*A*C)
    b=(-sqrt(complex(2*(A*E^2+C*D^2-B*D*E+(B^2-4*A*C)*F)*((A+C)-sqrt((A-C)^2+B^2)))))/(B^2-4*A*C)
    return a,b
end
