#ocn=ocean
#pcc=ocn.pccs[2]
#**
function opt_hvOSSplacement_zero_size(ocn,pcc)
    hv_circs=ocn.circuits
    owpps_tbl_hv=top_hvTopos(ocn.owpps)
    owpps_tbl_hv=owpps_tbl_hv[end:-1:1,end:-1:1]
    oss_system=circuit()
    oss_systems=Array{circuit,1}()
    discreteE=deepcopy(ocn.discretedom.edges)
    for indx0=1:length(owpps_tbl_hv[:,1])
    #indx0=5
        bus_dummies=bus[]
        for (indx1,owp) in enumerate(owpps_tbl_hv[indx0,:])
            if (owp==1)
                push!(bus_dummies,ocn.owpps[indx1])
            end
        end

        if (length(bus_dummies)>1)
            oss_system, discrete_dom, discrete_edges=opt_hvOssSystem_zero_size(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0,:])
        else
            oss_system=opt_str8Connect(bus_dummies[1],pcc,ocn,owpps_tbl_hv[indx0,:])
            discrete_dom=ocn.discretedom.nodes
            discrete_edges=ocn.discretedom.edges
        end
        #rnk=floor(Int32,top_bin2dec(owpps_tbl_hv[indx0]))
        #if (oss_system.cost<hv_circs[rnk].cost)
            #hv_circs[rnk]=deepcopy(oss_system)
            push!(hv_circs,oss_system)
            if (length(ocn.discretedom.nodes)<length(discrete_dom))
                ocn.discretedom.nodes=deepcopy(discrete_dom)
            end
            le=deepcopy(length(discreteE))
            if (le<length(discrete_edges))
                discreteE=opt_addEdges2temps(discreteE,discrete_edges,le)
            end
        #end
    end
    ocn.discretedom.edges=deepcopy(discreteE)
    return hv_circs
end
#=
function opt_hvOSSplacement_zero_size(ocn,pcc)
    hv_circs=ocn.circuits
    owpps_tbl_hv=top_hvTopos(ocn.owpps)
    owpps_tbl_hv=owpps_tbl_hv[end:-1:1,end:-1:1]
    oss_system=circuit()
    oss_systems=Array{circuit,1}()

    for indx0=1:length(owpps_tbl_hv[:,1])
    #indx0=5
        bus_dummies=bus[]
        for (indx1,owp) in enumerate(owpps_tbl_hv[indx0,:])
            if (owp==1)
                push!(bus_dummies,ocn.owpps[indx1])
            end
        end

        if (length(bus_dummies)>1)
            oss_system, discrete_dom, discrete_edges=opt_hvOssSystem_zero_size(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0,:])
        else
            oss_system=opt_str8Connect(bus_dummies[1],pcc,ocn,owpps_tbl_hv[indx0,:])
            discrete_dom=ocn.discretedom.nodes
            discrete_edges=ocn.discretedom.edges
        end
        #rnk=floor(Int32,top_bin2dec(owpps_tbl_hv[indx0]))
        #if (oss_system.cost<hv_circs[rnk].cost)
            #hv_circs[rnk]=deepcopy(oss_system)
            push!(hv_circs,oss_system)
            ocn.discretedom.nodes=deepcopy(discrete_dom)
            ocn.discretedom.edges=deepcopy(discrete_edges)
        #end
    end
    return hv_circs
end=#

#buses=bus_dummies
#bn=owpps_tbl_hv[indx0,:]
#**
function opt_hvOssSystem_zero_size(buses,pcc,ocn,bn)
    circ=circuit()
    owp=buses[1]
    opNode,domain_oss,edges_oss = opt_findOptPointHV_zero_size(buses,pcc,ocn,owp)
    circ,domain_oss,edges_oss=opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode,domain_oss,edges_oss)
    return circ,domain_oss,edges_oss
end

function opt_findOptPointHV_zero_size(buses,pcc,ocn,owp)
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    opNode=opt_1stNode(pcc,ocn,owp)
    lngth=opNode.H_cost
    rmnxys=opt_rmnNodeHV_zero_size(opNode,buses,ocn,pcc,owp)
    #rmnxys=rmnXys
    opNode=opt_OssOptimal(rmnxys,ocn.constrain.ellipses,opNode,lngth,pcc,ocn.yaxisMajor)
    bs,opNode=opt_ifIn2Out(opNode,ocn,owp)
    if bs==true
        opNode.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,opNode)
        lof_edgeifyOss(opNode,ocn,domain_oss,edges_oss)
    end
    return opNode,domain_oss,edges_oss
end

#frstNode=opNode
#**
function opt_rmnNodeHV_zero_size(frstNode,buses,ocn,pcc,owp)

    rmnXys=Array{xy,1}()

    for bs in buses
        push!(rmnXys,deepcopy(bs.node.xy))
    end
    return rmnXys
end

function opt_updateMVocean(ocn_mv,ocn_hv,bl)
    hve=deepcopy(ocn_hv.discretedom.edges)
    hve_l=length(hve)
    for edg in ocn_mv.discretedom.edges[bl+1:length(ocn_mv.discretedom.edges)]
        push!(hve,edg)
    end
    ocn_mv.discretedom.edges=deepcopy(hve)
    for (i,circ) in enumerate(ocn_mv.circuits)
        if (circ.cost<1)
            ocn_mv.circuits[i]=deepcopy(ocn_hv.circuits[i])
        end
    end
    return ocn_mv
end
########################################## Start MV System #####################################
#=
ocn=ocean
owpps=ocn.owpps
pcc=ocn.pccs[2]
=#
#**
function opt_mvOSSplacement_zero_size(ocn,owpps,pcc,hv_nodes)
    owpps_tbl=top_hvTopos(ocn.owpps)
    owpps_tbl=owpps_tbl[end:-1:1,end:-1:1]
    bus_dummies=Array{bus,1}()
    oss_system=circuit()
    oss_systems=ocn.circuits
    mvrng=cstF_mVrng(5,owpps[1].mva,owpps[1].wnd,ocn.finance,ocn.sys)
    dummy_sys=circuit()
    dummy_sys.cost=0
    ocn.discretedom.nodes=hv_nodes
    discreteE=deepcopy(ocn.discretedom.edges)
    #frst=1
    for indx0=1:length(owpps_tbl[:,1])
    #for indx0=1:length(owpps_tbl[251:253,1])
    #indx0=48
        bus_dummies=bus[]
        validMV,mv_square=opt_mvConstraints(ocn,owpps_tbl[indx0,:])

        if (validMV==true)
            #println("Valid: "*string(owpps_tbl[indx0,:]))
            for (indx1,owpp) in enumerate(owpps_tbl[indx0,:])
                if (owpp==1)
                    push!(bus_dummies,owpps[indx1])
                end
            end
            oss_system, discrete_dom, discrete_edges=opt_mvOssSystem_zero_size(mv_square,bus_dummies,pcc,ocn,owpps_tbl[indx0,:],mvrng)
            rnk=floor(Int32,top_bin2dec(owpps_tbl[indx0,:]))
            #if (oss_system.cost<oss_systems[rnk].cost)
                #oss_systems[rnk]=deepcopy(oss_system)
                push!(oss_systems,deepcopy(oss_system))
                if (length(ocn.discretedom.nodes)<length(discrete_dom))
                    ocn.discretedom.nodes=deepcopy(discrete_dom)
                end
                le=deepcopy(length(discreteE))
                if (le<length(discrete_edges))
                    discreteE=opt_addEdges2temps(discreteE,discrete_edges,le)
                end
                #ocn.discretedom.edges=deepcopy(discrete_edges)
            #end
        else
            println("Not Valid: "*string(owpps_tbl[indx0,:]))
            push!(oss_systems,dummy_sys)
        end
    end
    ocn.discretedom.edges=deepcopy(discreteE)
    return oss_systems
end

function opt_addEdges2temps(discreteE,discrete_edges,le)
    for edg in discrete_edges[le+1:length(discrete_edges)]
        push!(discreteE,deepcopy(edg))
    end
    return discreteE
end
#=function opt_mvOSSplacement_zero_size(ocn,owpps,pcc)
    owpps_tbl=top_hvTopos(ocn.owpps)
    owpps_tbl=owpps_tbl[end:-1:1,end:-1:1]
    bus_dummies=Array{bus,1}()
    oss_system=circuit()
    oss_systems=ocn.circuits
    mvrng=cstF_mVrng(5,owpps[1].mva,owpps[1].wnd,ocn.finance,ocn.sys)
    dummy_sys=circuit()
    dummy_sys.cost=0
    #frst=1
    for indx0=1:length(owpps_tbl[:,1])
    #indx0=48
        bus_dummies=bus[]
        validMV,mv_square=opt_mvConstraints(ocn,owpps_tbl[indx0,:])

        if (validMV==true)
            #println("Valid: "*string(owpps_tbl[indx0,:]))
            for (indx1,owpp) in enumerate(owpps_tbl[indx0,:])
                if (owpp==1)
                    push!(bus_dummies,owpps[indx1])
                end
            end
            oss_system, discrete_dom, discrete_edges=opt_mvOssSystem_zero_size(mv_square,bus_dummies,pcc,ocn,owpps_tbl[indx0,:],mvrng)
            rnk=floor(Int32,top_bin2dec(owpps_tbl[indx0,:]))
            #if (oss_system.cost<oss_systems[rnk].cost)
                #oss_systems[rnk]=deepcopy(oss_system)
                push!(oss_systems,deepcopy(oss_system))
                ocn.discretedom.nodes=deepcopy(discrete_dom)
                ocn.discretedom.edges=deepcopy(discrete_edges)
            #end
        else
            println("Not Valid: "*string(owpps_tbl[indx0,:]))
            push!(oss_systems,dummy_sys)
        end
    end
    return oss_systems
end=#
#buses=bus_dummies
#bn=owpps_tbl[indx0,:]
function opt_mvOssSystem_zero_size(mv_square,buses,pcc,ocn,bn,mvrng)
    circ=circuit()
    owp=opt_buses2nodes4MV_zero_size(mv_square,ocn.owpps,buses)
    opNode,domain_oss,edges_oss = opt_findOptPointMV_zero_size(buses,pcc,ocn,owp,mv_square,mvrng)
    circ,domain_oss,edges_oss=opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode,domain_oss,edges_oss)
    return circ,domain_oss,edges_oss
end
#=function opt_mvOssSystem_zero_size(mv_square,buses,pcc,ocn,bn,mvrng)
    circ=circuit()
    owp=opt_buses2nodes4MV_zero_size(mv_square,ocn.owpps,buses)
    opNode,domain_oss,edges_oss = opt_findOptPointMV_zero_size(buses,pcc,ocn,owp,mv_square,mvrng)
    circ,domain_oss,edges_oss=opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode,domain_oss,edges_oss)
    return circ,domain_oss,edges_oss
end=#


#all_opps=ocn.owpps
#**
function opt_buses2nodes4MV_zero_size(mv_square,all_opps,buses)
    target_owpp=1
    owp=deepcopy(buses[1])
    #=for (indx,bs) in enumerate(all_opps)
        if ((bs.node.xy.y<mv_square.ymn) || (bs.node.xy.x < mv_square.xmn) || (bs.num < owp.num))
            target_owpp=deepcopy(indx+1)
        end
    end
    owp=deepcopy(all_opps[target_owpp])=#
    return owp
end

#sqr=mv_square
function opt_findOptPointMV_zero_size(buses,pcc,ocn,owp,sqr,mvrng)
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    #opNode=opt_1stNode(pcc,ocn,owp)
    #lngth=opNode.H_cost
    rmnxys=opt_rmnNodeMV_zero_size(buses)
    #rmnxys=rmnXys
    opNode=opt_OssOptimalMV_zero_size(rmnxys,ocn.constrain.ellipses,owp,sqr,pcc,ocn.yaxisMajor,mvrng)
    bs,opNode=opt_ifIn2Out(opNode,ocn,owp)
    if bs==true
        opNode.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,opNode)
        lof_edgeifyOss(opNode,ocn,domain_oss,edges_oss)
    end
    return opNode,domain_oss,edges_oss
end

#=function opt_findOptPointMV_zero_size(buses,pcc,ocn,owp,sqr,mvrng)
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    #opNode=opt_1stNode(pcc,ocn,owp)
    #lngth=opNode.H_cost
    rmnxys=opt_rmnNodeMV_zero_size(buses)
    #rmnxys=rmnXys
    opNode=opt_OssOptimalMV_zero_size(rmnxys,ocn.constrain.ellipses,owp,sqr,pcc,ocn.yaxisMajor,mvrng)
    bs,opNode=opt_ifIn2Out(opNode,ocn,owp)
    if bs==true
        opNode.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,opNode)
        lof_edgeifyOss(opNode,ocn,domain_oss,edges_oss)
    end
    return opNode,domain_oss,edges_oss
end=#

function opt_rmnNodeMV_zero_size(buses)

    rmnXys=Array{xy,1}()

    for bs in buses
        push!(rmnXys,deepcopy(bs.node.xy))
    end
    return rmnXys
end
#=
xys=rmnxys
constrain=ocn.constrain.ellipses
yaxis=ocn.yaxisMajor
=#
function opt_OssOptimalMV_zero_size(xys,constrain,owp,sqr,pcc,yaxis,mvrng)
    mX = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(mX, x)
    @variable(mX, y)
    #@variable(m, lamda)
    #eps=1e-9

    @NLobjective(mX, Min,sum(((xys[i].x-x)^2+(xys[i].y-y)^2) for i in 1:length(xys)))
    for (ind,ellipse) in enumerate(constrain[1:length(constrain)])
        A=(((cos(ellipse.alpha)^2)/(ellipse.rx)^2)+((sin(ellipse.alpha)^2)/(ellipse.ry)^2))*(x-ellipse.x0)^2
        B=(2*cos(ellipse.alpha)*sin(ellipse.alpha)*((1/((ellipse.rx)^2))-(1/((ellipse.ry)^2))))*(x-ellipse.x0)*(y-ellipse.y0)
        C=(((sin(ellipse.alpha)^2)/(ellipse.rx)^2)+((cos(ellipse.alpha)^2)/(ellipse.ry)^2))*(y-ellipse.y0)^2
        #@NLconstraint(m, (x-ellipse.x0)^2/(ellipse.rx)^2 + (y-ellipse.y0)^2/(ellipse.ry)^2 >= 1.1)
        @constraint(mX,A+B+C >= 1)
    end
    @NLconstraint(mX, (x-xys[1].x)^2+(y-xys[1].y)^2 <= ((mvrng)^2))
    #@NLconstraint(mX, (x-xys[2].x)^2+(y-xys[2].y)^2 <= ((mvrng)^2))

    @constraint(mX, y <= sqr.ymx)
    @constraint(mX, y >= sqr.ymn)
    @constraint(mX, x <= sqr.xmx)
    @constraint(mX, x >= sqr.xmn)


    optimize!(mX)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end

function opt_hv2mvCircuit(ocn)
    for (i,crc) in enumerate(ocn.circuits)
        if (sum(x->x>0,crc.binary)>1)
            ocn.circuits[i]=circuit()
        end
    end
    return ocn
end
