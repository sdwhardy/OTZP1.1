using Plots, JuMP, LinearAlgebra, Polyhedra, SetProg, Mosek, MosekTools, MathOptInterface
using TypedPolynomials
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
poly_gone = rec_1x2
#poly_gone = complex_shape
ngs=nogo()
mosek_solver = with_optimizer(Mosek.Optimizer, QUIET=false)
ngs.nbnd,ngs.ebnd,ngs.sbnd,ngs.wbnd,ngs.bndryPnts=lof_bndNESW_nogo(ocean.nogos[2].bndryPnts)
simplex=opt_makeHalfSpace(ngs)
poly_gone=simplex

#Find largest area ellipse within polygone
Elips=findLargestEllipse_testing(poly_gone,mosek_solver)

toPlot=@time processDual(value(Elips))
#plot(toPlot)
#plot([XY0.x XY1.x],[XY0.y XY1.y])
#gui()
dia_Major,XY0,XY1=get_MajorAxis(toPlot)
rad_Minor,rad_Major,cntrXY,xyMinor,xyMajor,offset=get_ElipseData(toPlot,XY0,XY1,dia_Major)
plot!(polyhedron(poly_gone), ratio=1)
plot!(toPlot,color=:red,linewidth=4, label="")
gui()
plot!([cntrXY.x,xyMajor.x],[cntrXY.y,xyMajor.y])
plot!([cntrXY.x,xyMinor.x],[cntrXY.y,xyMinor.y])


#plot polygone
plot(polyhedron(poly_gone), ratio=1)
#plot!(Elips.variable)
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
function findLargestEllipse_testing(poly_g,_solver)
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



function get_ElipseData(tp,maXY0,maXY1,maDia)
    str8minr=line()
    centre_xy=xy()
    mn_xy=xy()
    offset=0
    centre_xy.x=(maXY0.x+maXY1.x)/2
    centre_xy.y=(maXY0.y+maXY1.y)/2
    rad=Inf
    xy_mn=xy()
    for (indx,x) in enumerate(tp[1])
        xy_mn.x=x
        xy_mn.y=tp[2][indx]
        xy_dist_temp=lof_pnt2pnt_dist(centre_xy,xy_mn)
        if (xy_dist_temp<rad)
            rad=deepcopy(xy_dist_temp)
        end
    end
    #check for a circle
    if (!((rad <= maDia/2+0.001) && (rad >= maDia/2-0.001)))
        nd_tail=node()
        nd_head=node()
        nd_tail.xy=maXY0
        nd_head.xy=maXY1
        verln=lof_lineDirection(nd_tail,nd_head)
        str8majr=lof_getStr8line(nd_head,nd_tail)
        if verln==true
            str8minr.m_findy=-1*(1/str8majr.m_findy)
            str8minr.b_findy=centre_xy.y/(str8minr.m_findy*centre_xy.x)
            del_x=sqrt((rad^2)/(1+str8minr.m_findy^2))
            del_y=str8minr.m_findy*del_x
            mn_xy.x=centre_xy.x+del_x
            mn_xy.y=centre_xy.y+del_y
            if (maXY1.y<maXY0.y)
                temp=deepcopy(maXY1)
                maXY1=deepcopy(maXY0)
                maXY0=deepcopy(temp)
            end
            offset=atand((mn_xy.y-centre_xy.y)/(mn_xy.x-centre_xy.x))
        else
            str8minr.m_findx=-1*(1/str8majr.m_findx)
            str8minr.b_findx=centre_xy.x/(str8minr.m_findx*centre_xy.y)
            del_y=sqrt((rad^2)/(1+str8minr.m_findx^2))
            del_x=str8minr.m_findx*del_y
            mn_xy.x=centre_xy.x+del_x
            mn_xy.y=centre_xy.y+del_y
            if (maXY1.x<maXY0.x)
                temp=deepcopy(maXY1)
                maXY1=deepcopy(maXY0)
                maXY0=deepcopy(temp)
            end
            offset=atand((maXY1.y-centre_xy.y)/(maXY1.x-centre_xy.x))
        end
    #is a circle
    else
        mn_xy=deepcopy(centre_xy)
        maXY1=deepcopy(centre_xy)
        mn_xy.y=mn_xy.y+rad
        maXY1.x=maXY1.x+rad
        offset=0
    end
    return rad,maDia/2,centre_xy,mn_xy,maXY1,offset
end

function get_MajorAxis(toPlot)
    xy_ma0=xy()
    xy_ma1=xy()
    dist=0
    for (indx0,x0) in enumerate(toPlot[1][1:length(toPlot[1])-1])
        xy0=xy()
        xy0.x=x0
        xy0.y=toPlot[2][indx0]
        for indx1=indx0+1:1:length(toPlot[1])
            xy1=xy()
            xy1.x=toPlot[1][indx1]
            xy1.y=toPlot[2][indx1]
            xy_dist_temp=lof_pnt2pnt_dist(xy0,xy1)
            if (xy_dist_temp>dist)
                dist=deepcopy(xy_dist_temp)
                xy_ma0.x=deepcopy(xy0.x)
                xy_ma0.y=deepcopy(xy0.y)
                xy_ma1.x=deepcopy(xy1.x)
                xy_ma1.y=deepcopy(xy1.y)
            end
        end
    end
    return dist,xy_ma0,xy_ma1
end

function processDual(set)
    @assert dimension(set) == 2
    # z is a halfspace of the primal so a ray of the dual
    npoints=128
    z = [1.0, 0.0, 0.0]
    h1, h2 = set.set.h
    # a is a ray of the primal so a halfspace of the dual
    a = [1, h1, h2]
    b = [h1, -1, 0]
    @assert abs(dot(a, b)) < 1e-8
    c = [h2 / (1 + h1^2), h1*h2 / (1 + h1^2), -1]
    @assert abs(dot(b, c)) < 1e-8
    @assert abs(dot(a, c)) < 1e-8
    polyhedron = dual_contour(scaling_function(set), npoints, Float64,z, b, c, true)
    return Polyhedra.planar_contour(fixandeliminate(polyhedron, 1, 1.0))
end

function dual_contour(f::Function, nhalfspaces::Int, ::Type{T},
                      point::Vector{T} = [0.0, 0.0],
                      x_axis::Vector{T} = [1.0, 0.0],
                      y_axis::Vector{T} = [0.0, 1.0],
                      cone = false) where T
    h = hrep(Polyhedra.HyperPlane{T, Vector{T}}[],
             Polyhedra.HalfSpace{T, Vector{T}}[], d=length(x_axis))
    for α in range(0, stop=2π - 2π/nhalfspaces, length=nhalfspaces)
        ray = x_axis * cos(α) + y_axis * sin(α)
        λ = f(ray...)
        a = point + ray / λ
        intersect!(h, HalfSpace(cone ? -a : a, cone ? zero(T) : one(T)))
    end
    return polyhedron(h)
end


function scaling_function(set)
    @assert length(space_variables(set)) == 2
    vars = [perspective_variable(set); space_variables(set)]
    # z is a halfspace of the primal so a ray of the dual
    z = [1.0, 0.0, 0.0]
    in_set(Δ::Vector) = perspective_gauge0(set)(vars => z + Δ) < 0
    @assert in_set(zeros(3))
    return (Δz, Δx, Δy) -> begin
        Δ = [Δz, Δx, Δy]
        _in_set(λ::Real) = in_set(Δ * λ)
        λ = 1.0
        while _in_set(λ)
            if λ > 1e10
                error("Error in plotting : the `InteriorPoint` seems to be on the boundary")
            end
            λ *= 2
        end
        λmin = 0.0
        λmax = λ
        # Binary search. Invariant: in_set(λmin) and !in_set(λmax)
        while abs(λmin - λmax) > 1e-8
            λ = (λmin + λmax) / 2
            if _in_set(λ)
                λmin = λ
            else
                λmax = λ
            end
        end

        λ = (λmin + λmax) / 2
        return 1 / λ
    end
end


function dimension(set)
    return length(space_variables(set))
end

function space_variables(set)
    return set.set.x
end

function perspective_variable(set)
    return set.set.z
end
function perspective_gauge0(set)
    return set.set.p
end
