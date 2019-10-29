
#=
ocn=ocean
owpps=ocn.owpps
pcc=ocn.pccs[1]
paths=opt_mvOSSplacement(ocn,owpps,pcc)
=#



function opt_mvOSSplacement(ocn,owpps,pcc)
    #owpps is an array of owpps ordered closest to farthest from the designated pcc
    owpps_tbl=top_mvTopos(owpps)
    bus_dummies=Array{bus,1}()
    oss_system=Array{node,1}()
    oss_systems=Array{Array{node,1},1}()
    for indx0=1:length(owpps_tbl[:,1])
    #for indx0=1:5
        bus_dummies=bus[]
        mv_square=opt_mvConstraints(ocn,owpps_tbl[indx0,:])
        for (indx1,owpp) in enumerate(owpps_tbl[indx0,:])
            if (owpp==1)
                push!(bus_dummies,owpps[indx1])
            end
        end
        oss_system=opt_mvOssSystem(mv_square,bus_dummies,pcc,ocn)
        push!(oss_systems,deepcopy(oss_system))
    end
    return oss_systems
end

#buses=bus_dummies

function opt_mvOssSystem(mv_square,buses,pcc,ocn)
    #each owpp path
    paths=node[]
    for owpp in buses
        push!(paths,deepcopy(as_Astar(owpp.node,pcc.node,ocn.discretedom.nodes)))
    end

    #pick first node in each owpp path
    xys=xy[]
    push!(xys,pcc.node.xy)
    for path=1:length(paths)
        pathsParent=paths[path].parent
        while pathsParent.parent.num != buses[path].node.num
            pathsParent=pathsParent.parent
        end
        push!(xys,deepcopy(pathsParent.xy))
    end

    paths=node[]
    owp=opt_buses2nodes4MVoptLocal(mv_square,buses)
    oss_node=opt_mvOssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    push!(paths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    for owpp in buses
        push!(paths,deepcopy(as_Astar(domain_oss[owpp.node.num],oss_node,domain_oss)))
    end

    #pick first node in each owpp path
    xys=xy[]
    push!(xys,pcc.node.xy)
    #for path=2:length(pathsA)
    for path=2:length(paths)
        pathsParent=paths[path].parent
        while pathsParent.parent.num != buses[path-1].node.num
            pathsParent=pathsParent.parent
        end
        push!(xys,deepcopy(pathsParent.xy))
    end

    #re calculate the path
    paths=node[]
    oss_node=opt_mvOssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    push!(paths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    for owpp in buses
        push!(paths,deepcopy(as_Astar(domain_oss[owpp.node.num],oss_node,domain_oss)))
    end

    #pick first node in each owpp path
    xys=xy[]
    for path=1:length(paths)
        push!(xys,deepcopy(paths[path].parent.xy))
    end

    #re calculate the path
    paths=node[]
    oss_node=opt_mvOssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    push!(paths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    for owpp in buses
        push!(paths,deepcopy(as_Astar(domain_oss[owpp.node.num],oss_node,domain_oss)))
    end

    return paths
end

function opt_ifIn2Out(oss_node,ocn)
    bs,a=lof_test4pnt(oss_node.xy.x,oss_node.xy.y,ocn)
    dummy_dist=Float64
    bsf_dist=Inf
    bsf_indx=1
    if bs==false
        for (indx,pnt) in enumerate(a.nodes)
            dummy_dist=lof_pnt2pnt_dist(oss_node.xy,pnt.xy)
            if dummy_dist<bsf_dist
                bsf_dist=deepcopy(dummy_dist)
                bsf_indx=deepcopy(indx)
            end
        end
        oss_node=a.nodes[bsf_indx]
    else
    end
    return bs,oss_node
end

function opt_buses2nodes4MVoptLocal(mv_square,buses)
    target_owpp=1
    owp=bus()
    owp=buses[target_owpp]
    for (indx,bs) in enumerate(buses)
        if buses[target_owpp].node.xy.y<mv_square.ymn
            target_owpp=deepcopy(indx)
            owp=buses[target_owpp]
        end
    end
    return owp
end

function opt_mvOssOptimalLocal(xys,constrain,owp,num_owpp)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(m, x)
    @variable(m, y)
    @variable(m, lamda)
    eps=1e-9

    #@NLconstraint(m, sum((x-xys[i].x)/(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys))==2*x*lamda)
    #@NLconstraint(m, sum((y-xys[i].y)/(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys))==2*y*lamda)
    @NLobjective(m, Min,sum((sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys)))
    #@NLconstraint(m, (x-constrain.x0)^2/(constrain.rx)^2 + (y-constrain.y0)^2/(constrain.ry)^2 == 1)
    #@NLconstraint(m, (x-constrain.x0)^2/(constrain.rx)^2 + (y-constrain.y0)^2/(constrain.ry)^2 >= 1)
    #@NLconstraint(m, (x-constrain.x0)^2/(owpp.zone.pos_width)^2 + (y-constrain.y0)^2/(constrain.ry)^2 <= 1)
    for ellipse in constrain[1:length(constrain)-num_owpp]
        #println(ellipse)
        @NLconstraint(m, (x-ellipse.x0)^2/(ellipse.rx)^2 + (y-ellipse.y0)^2/(ellipse.ry)^2 >= 1)
    end
        @NLconstraint(m, (x-constrain[length(constrain)-num_owpp+owp.num].x0)^2/(constrain[length(constrain)-num_owpp+owp.num].rx)^2 + (y-constrain[length(constrain)-num_owpp+owp.num].y0)^2/(constrain[length(constrain)-num_owpp+owp.num].ry)^2 >= 1)
    #println(mv_square.ymx)
    #=println(mv_square.ymn)
    println(mv_square.xmx)
    println(mv_square.xmn)

    println(xys)=#
    #@constraint(m, x >= mv_square.xmn)
    @constraint(m, y == owp.node.xy.y)
    #@constraint(m, x <= mv_square.xmx)
    #@constraint(m, y <= mv_square.ymx)
    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy

    #node_oss=opt_findClosestNode(temp_xy,owpp.zone.nodes)
    return node_oss
end

function opt_findClosestNode(x_y,ndes)
    node_oss=node()
    distance=Inf
    oss_indx=Inf

    for (indx,nd) in enumerate(ndes)
        temp_dist=lof_pnt2pnt_dist(x_y,nd.xy)
        if(temp_dist<distance)
            distance=deepcopy(temp_dist)
            oss_indx=deepcopy(indx)
        end
    end
    node_oss=ndes[oss_indx]
    return node_oss
end

#MV constraints
function opt_mvConstraints(ocn, owpps)
    sqr_bnd=square()
    sqr_bnd.ymn=0
    sqr_bnd.ymx=Inf
    sqr_bnd.xmn=0
    sqr_bnd.xmx=Inf
    for (indx,owpp) in enumerate(ocn.owpps)

        #Check MV ranges
        if (owpps[indx]==1)
            ymn=owpp.node.xy.y-owpp.mv_zone.neg_height
            ymx=owpp.node.xy.y+owpp.mv_zone.pos_height
            xmn=owpp.node.xy.x-owpp.mv_zone.neg_width
            xmx=owpp.node.xy.x+owpp.mv_zone.pos_width

            #south boundary
            if (ymn>sqr_bnd.ymn)
                sqr_bnd.ymn=deepcopy(ymn)
            end

            #north boundary
            if (ymx<sqr_bnd.ymx)
                sqr_bnd.ymx=deepcopy(ymx)
            end

            #west boundary
            if (xmn>sqr_bnd.xmn)
                sqr_bnd.xmn=deepcopy(xmn)
            end

            #east boundary
            if (xmx<sqr_bnd.xmx)
                sqr_bnd.xmx=deepcopy(xmx)
            end
        else
        end
    end
    return sqr_bnd
end

#OWWP constraints
function opt_owppConstraints(ocn)
    #@NLconstraint(m, ((x-x1)/a)^2+((y-y1)/b)^2 > 1)
    for owpp in ocn.owpps
        ellipse_con=ellipse()
        ellipse_con.y0=owpp.node.xy.y-(owpp.zone.neg_height-owpp.zone.pos_height)/2
        ellipse_con.x0=owpp.node.xy.x-(owpp.zone.neg_width-owpp.zone.pos_width)/2
        ellipse_con.ry=owpp.node.xy.y+owpp.zone.pos_height-ellipse_con.y0
        ellipse_con.rx=owpp.node.xy.x+owpp.zone.pos_width-ellipse_con.x0
        push!(ocn.constrain.ellipses,deepcopy(ellipse_con))
    end
end


#OWWP constraints
function opt_nogoConstraints(ocn)
    #@NLconstraint(m, ((x-x1)/a)^2+((y-y1)/b)^2 > 1)
    for nogo in ocn.nogos
        ellipse_con=ellipse()
        cntr=opt_findShapeCentre(nogo)
        ellipse_con.y0=cntr.y
        ellipse_con.x0=cntr.x
        ellipse_con.ry=opt_findShapeRadius(nogo,cntr)
        ellipse_con.rx=ellipse_con.ry
        push!(ocn.constrain.ellipses,deepcopy(ellipse_con))
    end
end
#=
ng=ocean.nogos[1]
=#
#Finds closest edge
function opt_findShapeRadius(ng,cntr)
    rad=Inf
    for bnds in [ng.nbnd,ng.sbnd,ng.ebnd,ng.wbnd]
        for bnd in bnds
            x=(cntr.x+(cntr.y*bnd.m_findx+bnd.b_findx))/2
            y=(cntr.y+(cntr.x*bnd.m_findy+bnd.b_findy))/2
            rad_dummy=sqrt((cntr.x-x)^2+(cntr.y-y)^2)
            if (rad_dummy<rad)
                rad=deepcopy(rad_dummy)
            end
        end
    end
    return rad
end

#=
ocean.owpps[2].node
mean=opt_findShapeCentre(ocean.owpps[2].zone)
=#
function opt_findShapeCentre(shape_bnds)
    mn0=opt_shapeMean(shape_bnds)
    triangles=Array{xy,1}()
    for bnd in shape_bnds.nbnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.xmn
        xy0.y=bnd.xmn*bnd.m_findy+bnd.b_findy
        xy1.x=bnd.xmx
        xy1.y=bnd.xmx*bnd.m_findy+bnd.b_findy
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    for bnd in shape_bnds.sbnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.xmn
        xy0.y=bnd.xmn*bnd.m_findy+bnd.b_findy
        xy1.x=bnd.xmx
        xy1.y=bnd.xmx*bnd.m_findy+bnd.b_findy
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    for bnd in shape_bnds.ebnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.ymn*bnd.m_findx+bnd.b_findx
        xy0.y=bnd.ymn
        xy1.x=bnd.ymx*bnd.m_findx+bnd.b_findx
        xy1.y=bnd.ymx
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    for bnd in shape_bnds.wbnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.ymn*bnd.m_findx+bnd.b_findx
        xy0.y=bnd.ymn
        xy1.x=bnd.ymx*bnd.m_findx+bnd.b_findx
        xy1.y=bnd.ymx
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    cntr=xy()
    cntr.x=0
    cntr.y=0
    for x_y in triangles
        cntr.x=cntr.x+deepcopy(x_y.x)
        cntr.y=cntr.y+deepcopy(x_y.y)
    end
    cntr.x=cntr.x/length(triangles)
    cntr.y=cntr.y/length(triangles)
    return cntr
end

function opt_shapeMean(shape_bnds)
    mean=xy()
    mean.y=0
    mean.x=0

    #north bnd
    for bnd in shape_bnds.nbnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end

    #south bnd
    for bnd in shape_bnds.sbnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end
    #west bnd
    for bnd in shape_bnds.wbnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end
    #east bnd
    for bnd in shape_bnds.ebnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end
    l=2*(length(shape_bnds.ebnd)+length(shape_bnds.wbnd)+length(shape_bnds.sbnd)+length(shape_bnds.nbnd))
    mean.y=mean.y/(l)
    mean.x=mean.x/(l)
    return mean
end
#=
y0=ocean.nogos[1].wbnd[1].ymn
y1=ocean.nogos[1].wbnd[1].ymx
y2=ocean.nogos[1].ebnd[1].ymn
x0=ocean.nogos[1].wbnd[1].m_findx*y0+ocean.nogos[1].wbnd[1].b_findx
x1=ocean.nogos[1].wbnd[1].m_findx*y1+ocean.nogos[1].wbnd[1].b_findx
x2=ocean.nogos[1].ebnd[1].m_findx*y2+ocean.nogos[1].ebnd[1].b_findx
pnt0=xy()
pnt1=xy()
pnt2=xy()
pnt0.x=x0
pnt0.y=y0
pnt1.x=x1
pnt1.y=y1
pnt2.x=x2
pnt2.y=y2
cntr=opt_findTriAngCentre(pnt0, pnt1, pnt2)
=#
function opt_findTriAngCentre(p0, p1, p2)
    cntr=xy()
    cntr.x=(p0.x+p1.x+p2.x)/3
    cntr.y=(p0.y+p1.y+p2.y)/3
    return cntr
end
