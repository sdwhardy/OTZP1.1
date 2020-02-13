
#**
function opt_makeHalfSpace(ng)
    hsAr=HalfSpace[]
    #southern
    #@constraint(m, l.m_findy*x-y <= -l.b_findy)
    for l in ng.sbnd
        xc=l.m_findy
        if abs(l.m_findy)<(1e-10)
            xc=0
        end
        push!(hsAr,HalfSpace([xc, -1], -l.b_findy))
    end
    #western
    #@constraint(m, -x+l.m_findx*y <= -l.b_findx)
    for l in ng.wbnd
        yc=l.m_findx
        if abs(l.m_findx)<(1e-10)
            yc=0
        end
        push!(hsAr,HalfSpace([-1, yc], -l.b_findx))
    end
    #northern
    #@constraint(m, -l.m_findy*x+y <= l.b_findy)
    for l in ng.nbnd
        xc=l.m_findy
        if abs(l.m_findy)<(1e-10)
            xc=0
        end
        push!(hsAr,HalfSpace([-xc, 1], l.b_findy))
    end
    #eastern
    #@constraint(m, x-l.m_findx*y <= l.b_findx)
    for l in ng.ebnd
        yc=l.m_findx
        if abs(l.m_findx)<(1e-10)
            yc=0
        end
        push!(hsAr,HalfSpace([1, -yc], l.b_findx))
    end
    simplex = hsAr[1]
    for hs in hsAr[2:length(hsAr)]
        simplex=simplex ∩ deepcopy(hs)
    end
    return simplex
end
#ocn=ocean
#**

function opt_rollUp(ocn)
    circs=ocn.circuits
    complete_systems=circuit[]
    push!(complete_systems,deepcopy(circs[length(circs)]))
    for (indx0,crc0) in enumerate(circs[1:floor(Int64,length(circs)/2)])
        bn0=findall(x->x==1,crc0.binary)
        for indx1=indx0+1:1:length(circs)
            bn1=findall(x->x==1,circs[indx1].binary)
            cmb, bn01=opt_willCombine(bn0,bn1,length(circs[indx1].binary))
            if (cmb)
                rnk=floor(Int64,top_bin2dec(bn01))
                if (rnk==length(circs))
                    push!(complete_systems,deepcopy(opt_replaceRnk(circs[rnk],crc0,circs[indx1])))
                elseif ((crc0.cost+circs[indx1].cost) < (circs[rnk].cost))
                    circs[rnk]=deepcopy(opt_replaceRnk(circs[rnk],crc0,circs[indx1]))
                end
            end
        end
        println("Circuit Number: "*string(indx0))
    end
    complete_systems=deepcopy(opt_sortFllSys(complete_systems))
    return complete_systems, circs
end

#**
function opt_sortFllSys(complete_systems)
    c_is=Array{Tuple{Float64,Int64},1}()
    sorted_circs=circuit[]
    for (i,sys) in enumerate(complete_systems)
        push!(c_is,deepcopy((sys.cost,i)))
    end
    sort!(c_is, by = x -> x[1])
    for cst in c_is
        push!(sorted_circs,deepcopy(complete_systems[cst[2]]))
    end
    return sorted_circs
end
#=
circ=circs[rnk]
crc1=circs[indx1]
=#
#**
function opt_replaceRnk(circ,crc0,crc1)
    crc01=deepcopy(circ)
    crc01.osss_owp=deepcopy(crc0.osss_owp)
    crc01.osss_mog=deepcopy(crc0.osss_mog)
    crc01.owp_MVcbls=deepcopy(crc0.owp_MVcbls)
    crc01.owp_HVcbls=deepcopy(crc0.owp_HVcbls)
    crc01.oss2oss_cbls=deepcopy(crc0.oss2oss_cbls)
    crc01.pcc_cbls=deepcopy(crc0.pcc_cbls)
    crc01.cost=deepcopy(crc0.cost+crc1.cost)
    for o_owp in crc1.osss_owp
        push!(crc01.osss_owp,deepcopy(o_owp))
    end
    for o_mog in crc1.osss_mog
        push!(crc01.osss_mog,deepcopy(o_mog))
    end
    for mv_c in crc1.owp_MVcbls
        push!(crc01.owp_MVcbls,deepcopy(mv_c))
    end
    for hv_c in crc1.owp_HVcbls
        push!(crc01.owp_HVcbls,deepcopy(hv_c))
    end
    for o2o_c in crc1.oss2oss_cbls
        push!(crc01.oss2oss_cbls,deepcopy(o2o_c))
    end
    for pcc_c in crc1.pcc_cbls
        push!(crc01.pcc_cbls,deepcopy(pcc_c))
    end
    return crc01
end
#**

#l=length(circs[indx1].binary)
function opt_willCombine(bn0,bn1,l)
    cmb=true
    bn01=zeros(Int8,l)
    for one0 in bn0
        for one1 in bn1
            if (one0==one1)
                cmb=false
            end
        end
    end
    if (cmb==true)
        for one0 in bn0
            bn01[one0]=1
        end
        for one1 in bn1
            bn01[one1]=1
        end
    end
    return cmb, bn01
end
########################################## Start Compound System ###############################
#ocn=ocean
#**
function opt_compoundOSS(ocn)
    owpps_dec_mhv=Array{Int32, 1}()
    owpps_tbl_mhv=Array{Array{Int8, 1}, 1}()
    owpp_tbl_mhv=Array{Int8, 1}()
    owpp_dec_mhv=Array{Int32, 1}()

    for crc in ocn.circuits
        push!(owpps_dec_mhv,deepcopy(crc.decimal))
        push!(owpps_tbl_mhv,deepcopy(crc.binary))
    end
#dec=54
#println("Testing mode!!!!!!")
    #for dec in owpps_dec_mhv[54:55]
    for dec in owpps_dec_mhv[7:length(owpps_dec_mhv)]
        owpp_tbl_mhv=owpps_tbl_mhv[dec]
        println("dec: "*string(owpp_tbl_mhv))
        min_Bin=findlast(x->x==1,owpp_tbl_mhv)
        osss_Pos=findfirst(x->x==1,owpp_tbl_mhv)
        zeros_mhv=findall(x->x==0,owpp_tbl_mhv[1:length(owpp_tbl_mhv)])
        min_dec=floor(Int32,2^(min_Bin-1))
        oss_dec=floor(Int32,2^(osss_Pos-1))
        for bin in reverse(owpps_tbl_mhv[min_dec+1:dec-1])

            #bin=reverse(owpps_tbl_mhv[min_dec+1:dec-1])[6]
            #println("bin: "*string(bin))
            owpp_dec_mhv=Int32[]
            if (findall(x->x==1,bin[1:osss_Pos]) == [])
                if (sum(bin[zeros_mhv]) < 1)
                    push!(owpp_dec_mhv,deepcopy(oss_dec))
                    dc=floor(Int32,top_bin2dec(bin))
                    push!(owpp_dec_mhv,deepcopy(dc))
                    if (sum(owpp_dec_mhv) == dec)
                        opt_bS(owpp_dec_mhv,ocn)
                    elseif (sum(owpp_dec_mhv) < dec)
                        binner_array=deepcopy(reverse(owpps_tbl_mhv[oss_dec+1:(dec-sum(owpp_dec_mhv))]))
                        dinner_array=Int32[]
                        for binner in binner_array
                            if (findall(x->x==1,binner[1:osss_Pos]) == [])
                                if (sum(binner[zeros_mhv]) < 1)
                                    push!(dinner_array,top_bin2dec(binner))
                                end
                            end
                        end
                        if (length(dinner_array) != 0)
                            opt_breakdownSystem(owpp_dec_mhv,dinner_array,dec,ocn)
                        else
                            println("Error! Decimal array is empty")
                        end

                    else
                        println("Ërror! Sum of binaries exceeds maximum.")
                    end
                end
            end
        end
    end
end
#**
#i=1
#dc=dinner_array[i]
function opt_breakdownSystem(owpp_dec_mhv,dinner_array,dec,ocn)
    for (i,dc) in enumerate(dinner_array)
        temp_dec_mhv=deepcopy(owpp_dec_mhv)
        push!(temp_dec_mhv,deepcopy(dc))
        if (sum(temp_dec_mhv) == dec)
            opt_bS(temp_dec_mhv,ocn)
        #elseif (sum(owpp_dec_mhv) < dec)# original
        elseif (sum(temp_dec_mhv) < dec)
            opt_breakdownSystem(temp_dec_mhv,dinner_array[1:i-1],dec,ocn)
        else
            #println("Error: decimal sum is over maximum!")
        end
    end
end
#**
#owpp_dec_mhv=temp_dec_mhv
function opt_bS(owpp_dec_mhv,ocn)
    oss_systems=Array{circuit,1}()
    rnk=sum(owpp_dec_mhv)

    ioPos=Int32[]
    #collect systems to sum
    for crc in owpp_dec_mhv
        #println(crc)
        push!(oss_systems,ocn.circuits[crc])
        ios=findall(x->x==1,oss_systems[length(oss_systems)].binary)
        for io in ios
            push!(ioPos,deepcopy(io))
        end
    end
    ioConcise=unique(ioPos)
    if (length(ioConcise) == length(ioPos))
        opt_buildSystem(oss_systems,rnk,ocn)
    else
    end
end
#**
function opt_compoundCbl(circ,oss_system,ocn)

    for cb in circ.owp_MVcbls
        push!(oss_system.owp_MVcbls,cb)
    end

    for cb in circ.owp_HVcbls
        push!(oss_system.owp_HVcbls,cb)
    end

    for cb in circ.oss2oss_cbls
        push!(oss_system.oss2oss_cbls,cb)
    end

    for xm in circ.osss_owp
        push!(oss_system.osss_owp,xm)
    end

    for xm in circ.osss_mog
        push!(oss_system.osss_mog,xm)
    end

    #find oss to mog connections
    pth=deepcopy(as_Astar(circ.osss_mog[1].node,oss_system.osss_mog[1].node,ocn.discretedom.nodes))
    oss2ossCbl=cstF_Compound_HvCblo2o(pth.G_cost,circ.oss_mva,circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt,circ.oss_wind,ocn.finance,circ.pcc_cbls[length(circ.pcc_cbls)].mva,circ.pcc_cbls[length(circ.pcc_cbls)].size,circ.pcc_cbls[length(circ.pcc_cbls)].num)
    return oss2ossCbl,pth
end

#**
function opt_combineScratch(crc,oss_system,ocn)
    circ=circuit()
    opNode,domain_oss,edges_oss = opt_findOptPoint(crc.owpps,oss_system.osss_mog[1],ocn,crc.base_owp)
    circ,domain_oss,edges_oss=opt_mvhvEquipment(circ,crc.base_owp,crc.owpps,oss_system.osss_mog[1],ocn,crc.binary,opNode,domain_oss,edges_oss)
    ocn.discretedom.nodes=deepcopy(domain_oss)
    ocn.discretedom.edges=deepcopy(edges_oss)
    oss2ossCbl,pth=opt_compoundCbl(circ,oss_system,ocn)
    return oss2ossCbl,pth
end

#**
function opt_combineCblRoute(pth,oss2ossCbl,oss_system)
    currentNode=deepcopy(pth)
    goalNode=deepcopy(pth.goal)
    while currentNode.num != goalNode
        push!(oss2ossCbl.pth,deepcopy(currentNode))
        currentNode=currentNode.parent
    end
    push!(oss2ossCbl.pth,deepcopy(currentNode))
    oss2ossCbl.pth=reverse!(oss2ossCbl.pth)
    push!(oss_system.oss2oss_cbls,deepcopy(oss2ossCbl))
    return oss2ossCbl,oss_system
end

#**
function opt_combineCblPW(oss_system,circ,pHv,wHv,p132,w132,p220,w220,p400,w400)
    #store wind for mog transformer calcs
    if (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
        push!(pHv,circ.oss_mva)
        push!(wHv,circ.oss_wind)
    elseif (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==132)
        push!(p132,circ.oss_mva)
        push!(w132,circ.oss_wind)
    elseif (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==220)
        push!(p220,circ.oss_mva)
        push!(w220,circ.oss_wind)
    elseif (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==400)
        push!(p400,circ.oss_mva)
        push!(w400,circ.oss_wind)
    else
    end
    return pHv,wHv,p132,w132,p220,w220,p400,w400
end

#**
function opt_fromScrtch(dest,crc,ax)
    Scrtch=false
    if (ax==true)
    #yaxis Major
    #println(crc.decimal)
        if ((dest.node.xy.x-crc.base_owp.node.xy.x)*(crc.owp_MVcbls[length(crc.owp_MVcbls)].pth[2].xy.x-crc.base_owp.node.xy.x)<0)
            Scrtch=true
        end
    else
    #xaxis Major
        if ((dest.node.xy.y-crc.base_owp.node.xy.y)*(crc.owp_MVcbls[length(crc.owp_MVcbls)].pth[2].xy.y-crc.base_owp.node.xy.y)<0)
            Scrtch=true
        end
    end
    return Scrtch
end

#**
function opt_buildSystem(oss_systems,rnk,ocn)
    oss_system=circuit()

    #take base information
    oss_system.binary=ocn.circuits[rnk].binary
    oss_system.decimal=ocn.circuits[rnk].decimal
    oss_system.owpps=ocn.circuits[rnk].owpps
    oss_system.base_owp=ocn.circuits[rnk].base_owp
    oss_system.oss_wind=ocn.circuits[rnk].oss_wind
    oss_system.oss_mva=ocn.circuits[rnk].oss_mva

    #Take base OWP,PCC connections and MOG without transformers
    push!(oss_system.owp_MVcbls,deepcopy(ocn.circuits[rnk].owp_MVcbls[1]))
    push!(oss_system.osss_mog,deepcopy(ocn.circuits[rnk].osss_mog[1]))
    oss_system.osss_mog[1].xfmrs=xfo[]
    push!(oss_system.pcc_cbls,deepcopy(ocn.circuits[rnk].pcc_cbls[1]))
    if (oss_system.owp_MVcbls[1].pth[length(oss_system.owp_MVcbls[1].pth)].num != oss_system.osss_mog[1].node.num)
        push!(oss_system.owp_HVcbls,deepcopy(ocn.circuits[rnk].owp_HVcbls[1]))
        push!(oss_system.osss_owp,deepcopy(ocn.circuits[rnk].osss_owp[1]))
    else
    end
    #=base_sys=0
    for crc in oss_systems[2:length(oss_systems)]
        base_sys=base_sys+crc.cost-crc.pcc_cbls[1].costs.ttl+(crc.pcc_cbls[1].costs.ttl/crc.pcc_cbls[1].length)*lof_pnt2pnt_dist(crc.pcc_cbls[1].pth[1].xy,oss_system.osss_mog[1].node.xy)
    end
    opt_ttlMvCirc(oss_system)
    #println(string(mox+oss_system.cost+base_sys)*" - "*string(ocn.circuits[rnk].cost))
    if (oss_system.cost+base_sys <= (ocn.circuits[rnk].cost)*1)=#
        #println("inside!!!")
        #store wind for transformer calculation at MOG
        pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_InitPW()
        if (oss_system.owp_MVcbls[1].pth[length(oss_system.owp_MVcbls[1].pth)].num == oss_system.osss_mog[1].node.num)
            push!(pMv,oss_systems[1].oss_mva)
            push!(wMv,oss_systems[1].oss_wind)
        elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
            push!(pHv,oss_systems[1].oss_mva)
            push!(wHv,oss_systems[1].oss_wind)
        elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==132)
            push!(p132,oss_systems[1].oss_mva)
            push!(w132,oss_systems[1].oss_wind)
        elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==220)
            push!(p220,oss_systems[1].oss_mva)
            push!(w220,oss_systems[1].oss_wind)
        elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==400)
            push!(p400,oss_systems[1].oss_mva)
            push!(w400,oss_systems[1].oss_wind)
        else
        end
        combinesSyss=Float64[]
        push!(combinesSyss,oss_systems[1].decimal)
    ############################# update this section with new set calculation ########################
        if (length(oss_systems)>1)
            connect_Type=""
            #crc=oss_systems[2:length(oss_systems)][2]
            for crc in oss_systems[2:length(oss_systems)]
                push!(combinesSyss,crc.decimal)
                #copy base cables and transformers
                if (length(crc.owpps)==1)
                    oss_system,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_str8Oss2Oss(crc,oss_system,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
                else

                    if (opt_fromScrtch(oss_system.osss_mog[1],crc,ocn.yaxisMajor))#check if full recalc required
                        connect_Type=connect_Type*"_Scratch"
                        oss2ossCbl,pth=opt_compoundCbl(crc,oss_system,ocn)
                        #oss2ossCbl,pth=opt_combineScratch(crc,oss_system,ocn)
                    else
                        connect_Type=connect_Type*"_Orig"
                        oss2ossCbl,pth=opt_compoundCbl(crc,oss_system,ocn)
                    end
                    #println(string(crc.pcc_cbls[1].num)*" - "*string(crc.pcc_cbls[1].size)*" || "*string(oss2ossCbl.num)*" - "*string(oss2ossCbl.size))
                    oss2ossCbl,oss_system=opt_combineCblRoute(pth,oss2ossCbl,oss_system)
                    pHv,wHv,p132,w132,p220,w220,p400,w400=opt_combineCblPW(oss_system,crc,pHv,wHv,p132,w132,p220,w220,p400,w400)
                end
            end

            oss_system.osss_mog[1]=opt_mogXfmrs(oss_system.osss_mog[1],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
            opt_ttlMvCirc(oss_system)
            #println(string(oss_system.cost)*" - "*string(ocn.circuits[rnk].cost)*" - "*string(combinesSyss))
            if (oss_system.cost < ocn.circuits[rnk].cost)
                #println(string(rnk)*" - "*connect_Type)
                ocn.circuits[rnk]=deepcopy(oss_system)
            end
        end
    #else
        #println("skip!!!!")
    #end
end

########################################## End Compound System #################################
########################################## Start HV System #####################################
#ocn=ocean
#pcc=ocn.pccs[2]
#**
function opt_hvOSSplacement(ocn,pcc)
    hv_circs=ocn.circuits
    #owpps_tbl_hv=Array{Array{Int8, 1}, 1}()
    owpps_tbl_hv=top_hvTopos(ocn.owpps)
    owpps_tbl_hv=owpps_tbl_hv[end:-1:1,end:-1:1]
    oss_system=circuit()
    oss_systems=Array{circuit,1}()
    #=for crc in mv_circs
        if ((length(findall(x->x==true,crc.binary)) != 1) && (crc.base_owp.num != crc.owpps[1].num))
            push!(owpps_tbl_hv,crc.binary)
        end
    end=#
    for indx0=1:length(owpps_tbl_hv[:,1])
    #indx0=5
        bus_dummies=bus[]
        for (indx1,owp) in enumerate(owpps_tbl_hv[indx0,:])
            if (owp==1)
                push!(bus_dummies,ocn.owpps[indx1])
            end
        end
        #oss_system, discrete_dom, discrete_edges=opt_hvOssSystem(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0,:],side)
        if (length(bus_dummies)>1)
            oss_system, discrete_dom, discrete_edges=opt_hvOssSystem2(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0,:])
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
end

#=
buses=bus_dummies
bn=owpps_tbl[indx0,:]
=#
#**
function opt_mvOssSystem2(mv_square,buses,pcc,ocn,bn)
    circ=circuit()
    owp=opt_buses2nodes4MVoptLocal(mv_square,ocn.owpps,buses)
    #opNode,domain_oss,edges_oss = opt_findOptPointMV2(buses,pcc,ocn,owp,mv_square)
    opNode,domain_oss,edges_oss = opt_findOptPointMV(buses,pcc,ocn,owp,mv_square)
    circ,domain_oss,edges_oss=opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode,domain_oss,edges_oss)
    return circ,domain_oss,edges_oss
end

#buses=bus_dummies
#bn=owpps_tbl_hv[indx0,:]
#**
function opt_hvOssSystem2(buses,pcc,ocn,bn)
    circ=circuit()
    owp=buses[1]
    opNode,domain_oss,edges_oss = opt_findOptPoint(buses,pcc,ocn,owp)
    circ,domain_oss,edges_oss=opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode,domain_oss,edges_oss)
    return circ,domain_oss,edges_oss
end
#**
function opt_findWindPths2(pcc,opt_node,domain_oss,buses)
    #each owpp path
    opt_paths=node[]
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float32,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float32,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0

    push!(opt_paths,deepcopy(as_Astar(pcc.node,opt_node,domain_oss)))
    #opp=buses[2]
    for opp in buses
        push!(opt_paths,deepcopy(as_Astar(opp.node,opt_node,domain_oss)))
        wind_sum.pu=(wind_sum.pu.+opp.wnd.pu)
        wind_sum.ce=(wind_sum.ce.+opp.wnd.ce)
        wind_sum.delta=(wind_sum.delta+opp.wnd.delta)
        wind_sum.lf=(wind_sum.lf+opp.wnd.lf)
        power_sum=(power_sum+opp.mva)
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)
    return power_sum,wind_sum,opt_paths
end
#**
function opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode,domain_oss,edges_oss)
    circ.base_owp=owp
    #opNode=owp
    as_pths=node[]
    #if (length(buses)>1)
    circ.oss_mva,circ.oss_wind,as_pths=opt_findWindPths2(pcc,opNode,domain_oss,buses)

    #Find PCC connection cable  and PCC transformer
    println("mva: "*string(circ.oss_mva))
    cbl_xfo=cstF_HvCblallKvo2p(as_pths[1].G_cost,circ.oss_mva,circ.oss_wind,ocn.finance,pcc)

    #Set cable path
    currentNode=deepcopy(as_pths[1])
    goalNode=as_pths[1].goal
    push!(cbl_xfo[1].pth,deepcopy(currentNode))
    while currentNode.num != goalNode
        push!(cbl_xfo[1].pth,deepcopy(currentNode.parent))
        currentNode=currentNode.parent
    end
    push!(circ.pcc_cbls,cbl_xfo[1])

    #pcc transformer
    circ.pcc=deepcopy(pcc)
    push!(circ.pcc.xfmrs,cbl_xfo[2])

    #MV cable sizes and paths
    pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_InitPW()
    for (i,owp) in enumerate(buses)
        cb,xs=cstF_MvHvCbloss(as_pths[i+1].G_cost,owp.mva,owp.wnd,ocn.finance,circ.pcc_cbls[1].elec.volt,ocn.sys)
        if (length(cb)==1)
            currentNode=deepcopy(as_pths[i+1])
            goalNode=currentNode.goal
            push!(cb[1].pth,deepcopy(currentNode))
            while (currentNode.num != goalNode)
                push!(cb[1].pth,deepcopy(currentNode.parent))
                currentNode=currentNode.parent
            end
            cb[1].pth=reverse!(cb[1].pth)
            push!(circ.owp_MVcbls,cb[1])
            #store winds and powers
            push!(pMv,owp.mva)
            push!(wMv,owp.wnd)
        elseif (length(cb)==2)
            #Cable
            currentNode=deepcopy(as_pths[i+1])
            goalNode=currentNode.goal
            push!(cb[2].pth,deepcopy(currentNode))
            while currentNode.parent.num != goalNode
                push!(cb[2].pth,deepcopy(currentNode.parent))
                currentNode=currentNode.parent
            end
            cb[2].pth=reverse!(cb[2].pth)
            push!(cb[1].pth,deepcopy(currentNode.parent))
            push!(cb[1].pth,deepcopy(currentNode))
            push!(circ.owp_MVcbls,deepcopy(cb[1]))
            push!(circ.owp_HVcbls,deepcopy(cb[2]))

            #OSS
            ossmv=bus()
            ossmv.node=currentNode
            ossmv.base_cost=ocn.finance.FC_bld
            push!(circ.osss_owp,ossmv)
            push!(circ.osss_owp[length(circ.osss_owp)].xfmrs,deepcopy(xs[1]))

            #store winds and powers
            if (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)
                push!(pHv,owp.mva)
                push!(wHv,owp.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 132)
                push!(p132,owp.mva)
                push!(w132,owp.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 220)
                push!(p220,owp.mva)
                push!(w220,owp.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 400)
                push!(p400,owp.mva)
                push!(w400,owp.wnd)
            end
        end
    end
#circ.osss_mog=[]
    #mog transformer
    ossMOG=bus()
    ossMOG.mva=circ.oss_mva
    ossMOG.wnd=circ.oss_wind
    ossMOG.node=as_pths[1]
    ossMOG.base_cost=ocn.finance.FC_bld
    push!(circ.osss_mog,ossMOG)
    circ.osss_mog[length(circ.osss_mog)]=opt_mogXfmrs(circ.osss_mog[length(circ.osss_mog)],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)

    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.owpps=buses
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    return circ,domain_oss,edges_oss
end

#**
function opt_findOptPoint(buses,pcc,ocn,owp)
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    opNode=opt_1stNode(pcc,ocn,owp)
    lngth=opNode.H_cost
    rmnxys=opt_rmnNode(opNode,buses,ocn,pcc,owp)
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
#sqr=mv_square
function opt_findOptPointMV(buses,pcc,ocn,owp,sqr)
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    opNode=opt_1stNode(pcc,ocn,owp)
    #lngth=opNode.H_cost
    rmnxys=opt_rmnNodeMV(opNode,buses,ocn,pcc,owp)
    #rmnxys=rmnXys
    opNode=opt_OssOptimalMV2(rmnxys,ocn.constrain.ellipses,opNode,sqr,pcc,ocn.yaxisMajor)
    bs,opNode=opt_ifIn2Out(opNode,ocn,owp)
    if bs==true
        opNode.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,opNode)
        lof_edgeifyOss(opNode,ocn,domain_oss,edges_oss)
    end
    return opNode,domain_oss,edges_oss
end

function opt_findOptPointMV2(buses,pcc,ocn,owp,sqr)
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    eNode,wNode=opt_1stNodeMV(owp)
    rmnxys,opNode=opt_rmnNodeMV(eNode,wNode,buses,ocn,pcc,owp)
    lngth=opNode.H_cost
    opNode=opt_OssOptimalMV2(rmnxys,ocn.constrain.ellipses,opNode,sqr,pcc,ocn.yaxisMajor)
    bs,opNode=opt_ifIn2Out(opNode,ocn,owp)
    if bs==true
        opNode.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,opNode)
        lof_edgeifyOss(opNode,ocn,domain_oss,edges_oss)
    end
    return opNode,domain_oss,edges_oss
end
#=
xys=rmnxys
constrain=ocn.constrain.ellipses
bsf_node=opNode
yaxis=ocn.yaxisMajor
=#
#**

function opt_OssOptimalMV(xys,constrain,bsf_node,zne,pcc,yaxis)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))
    @variable(m, x)
    @variable(m, y)
    #@variable(m, lamda)
    eps=1e-9

    @NLobjective(m, Min,sum(((xys[i].x-x)^2+(xys[i].y-y)^2+eps) for i in 1:length(xys)))
    for (ind,ellipse) in enumerate(constrain[1:length(constrain)])
        A=(((cos(ellipse.alpha)^2)/(ellipse.rx)^2)+((sin(ellipse.alpha)^2)/(ellipse.ry)^2))*(x-ellipse.x0)^2
        B=(2*cos(ellipse.alpha)*sin(ellipse.alpha)*((1/((ellipse.rx)^2))-(1/((ellipse.ry)^2))))*(x-ellipse.x0)*(y-ellipse.y0)
        C=(((sin(ellipse.alpha)^2)/(ellipse.rx)^2)+((cos(ellipse.alpha)^2)/(ellipse.ry)^2))*(y-ellipse.y0)^2
        #@NLconstraint(m, (x-ellipse.x0)^2/(ellipse.rx)^2 + (y-ellipse.y0)^2/(ellipse.ry)^2 >= 1.1)
        @constraint(m,A+B+C >= 1)
    end

    if (yaxis==true)
        if (bsf_node.xy.y>pcc.node.xy.y)
            @constraint(m, y <= zne.node.xy.y)
        else
            @constraint(m, y >= zne.node.xy.y)
        end
    else
        if (bsf_node.xy.x>pcc.node.xy.x)
            @constraint(m, x <= zne.node.xy.x)
        else
            @constraint(m, x >= zne.node.xy.x)
        end
    end


    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end
#=xys=rmnxys
constrain=ocn.constrain.ellipses
bsf_node=opNode
yaxis=ocn.yaxisMajor
=#
function opt_OssOptimalMV2(xys,constrain,bsf_node,sqr,pcc,yaxis)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))
    @variable(m, x)
    @variable(m, y)
    #@variable(m, lamda)
    eps=1e-9

    @NLobjective(m, Min,sum(((xys[i].x-x)^2+(xys[i].y-y)^2+eps) for i in 1:length(xys)))
    for (ind,ellipse) in enumerate(constrain[1:length(constrain)])
        A=(((cos(ellipse.alpha)^2)/(ellipse.rx)^2)+((sin(ellipse.alpha)^2)/(ellipse.ry)^2))*(x-ellipse.x0)^2
        B=(2*cos(ellipse.alpha)*sin(ellipse.alpha)*((1/((ellipse.rx)^2))-(1/((ellipse.ry)^2))))*(x-ellipse.x0)*(y-ellipse.y0)
        C=(((sin(ellipse.alpha)^2)/(ellipse.rx)^2)+((cos(ellipse.alpha)^2)/(ellipse.ry)^2))*(y-ellipse.y0)^2
        #@NLconstraint(m, (x-ellipse.x0)^2/(ellipse.rx)^2 + (y-ellipse.y0)^2/(ellipse.ry)^2 >= 1.1)
        @constraint(m,A+B+C >= 1)
    end

    @constraint(m, y <= sqr.ymx)
    @constraint(m, y >= sqr.ymn)
    @constraint(m, x <= sqr.xmx)
    @constraint(m, x >= sqr.xmn)


    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end

function opt_OssOptimal(xys,constrain,bsf_node,lngth,pcc,yaxis)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))
    @variable(m, x)
    @variable(m, y)
    #@variable(m, lamda)
    eps=1e-9

    @NLobjective(m, Min,sum(((xys[i].x-x)^2+(xys[i].y-y)^2+eps) for i in 1:length(xys)))
    for (ind,ellipse) in enumerate(constrain[1:length(constrain)])
        A=(((cos(ellipse.alpha)^2)/(ellipse.rx)^2)+((sin(ellipse.alpha)^2)/(ellipse.ry)^2))*(x-ellipse.x0)^2
        B=(2*cos(ellipse.alpha)*sin(ellipse.alpha)*((1/((ellipse.rx)^2))-(1/((ellipse.ry)^2))))*(x-ellipse.x0)*(y-ellipse.y0)
        C=(((sin(ellipse.alpha)^2)/(ellipse.rx)^2)+((cos(ellipse.alpha)^2)/(ellipse.ry)^2))*(y-ellipse.y0)^2
        #@NLconstraint(m, (x-ellipse.x0)^2/(ellipse.rx)^2 + (y-ellipse.y0)^2/(ellipse.ry)^2 >= 1.1)
        @constraint(m,A+B+C >= 1)
    end

    @NLconstraint(m, (x-pcc.node.xy.x)^2+(y-pcc.node.xy.y)^2 <= (lngth^2))
    if (yaxis==true)
        if (bsf_node.xy.y>pcc.node.xy.y)
            @constraint(m, y <= bsf_node.xy.y)
        else
            @constraint(m, y >= bsf_node.xy.y)
        end
    else
        if (bsf_node.xy.x>pcc.node.xy.x)
            @constraint(m, x <= bsf_node.xy.x)
        else
            @constraint(m, x >= bsf_node.xy.x)
        end
    end


    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end

#frstNode=opNode
#**
function opt_rmnNode(frstNode,buses,ocn,pcc,owp)

    rmnXys=Array{xy,1}()

    for bs in buses
        #bs=buses[2]
        if (bs.num != owp.num)
            nde=as_Astar(bs.node,frstNode,ocn.discretedom.nodes)
            nde=opt_backUpNode(nde,bs.node,ocn.yaxisMajor)
            push!(rmnXys,deepcopy(nde.xy))
        end
    end
    return rmnXys
end
function opt_rmnNodeMV(frstNode,buses,ocn,pcc,owp)

    rmnXys=Array{xy,1}()

    for bs in buses
        #bs=buses[2]
        #if (bs.num != owp.num)
            nde=as_Astar(bs.node,frstNode,ocn.discretedom.nodes)
            nde=opt_backUpNode(nde,bs.node,ocn.yaxisMajor)
            push!(rmnXys,deepcopy(nde.xy))
        #end
    end
    return rmnXys
end

function opt_rmnNodeMV_failed(eNode,wNode,buses,ocn,pcc,owp)
    rmnXysE=Array{xy,1}()
    rmnXysW=Array{xy,1}()
    le=0
    lw=0
    for bs in buses
        if (bs.num != owp.num)
            ndeE=as_Astar(bs.node,eNode,ocn.discretedom.nodes)
            le=le+ndeE.G_cost
            ndeE=opt_backUpNode(ndeE,bs.node,ocn.yaxisMajor)

            ndeW=as_Astar(bs.node,wNode,ocn.discretedom.nodes)
            lw=lw+ndeW.G_cost
            ndeW=opt_backUpNode(ndeW,bs.node,ocn.yaxisMajor)
            push!(rmnXysE,deepcopy(ndeE.xy))
            push!(rmnXysW,deepcopy(ndeW.xy))
        end
    end
    if (le<lw)
        rmnxys=deepcopy(rmnXysE)
        opNode=deepcopy(eNode)
    else
        rmnxys=deepcopy(rmnXysW)
        opNode=deepcopy(wNode)
    end
    return rmnxys,opNode
end


#**
function opt_backUpNode(nde,owp,axisY)
        xy_0=deepcopy(nde.xy)
        xy_1=deepcopy(nde.parent.xy)
        xy_2=deepcopy(nde.parent.parent.xy)
        nde=deepcopy(nde.parent)
        if (axisY == true)
            div0=xy_1.x-xy_0.x
            div1=xy_2.x-xy_1.x
            while ((nde.parent.num != owp.num) && (div0*div1 > 0))
                xy_0=deepcopy(nde.xy)
                xy_1=deepcopy(nde.parent.xy)
                xy_2=deepcopy(nde.parent.parent.xy)
                nde=deepcopy(nde.parent)
                div0=deepcopy(xy_1.x-xy_0.x)
                div1=deepcopy(xy_2.x-xy_1.x)
            end
        else
            div0=xy_1.y-xy_0.y
            div1=xy_2.y-xy_1.y
            while ((nde.parent.num != owp.num) && (div0*div1 > 0))
                xy_0=deepcopy(nde.xy)
                xy_1=deepcopy(nde.parent.xy)
                xy_2=deepcopy(nde.parent.parent.xy)
                nde=deepcopy(nde.parent)
                div0=deepcopy(xy_1.y-xy_0.y)
                div1=deepcopy(xy_2.y-xy_1.y)
            end
        end
    return nde
end

#**
function opt_1stNode(pcc,ocn,owp)
    first_pth=as_Astar(owp.node,pcc.node,ocn.discretedom.nodes)
    while first_pth.parent.num != owp.node.num
        first_pth=first_pth.parent
    end
    return first_pth
end

#**

function opt_1stNodeMV(owp)
    ePnt=node()
    wPnt=node()
    exy=deepcopy(owp.node.xy)
    wxy=deepcopy(owp.node.xy)
    exy.x=exy.x-owp.zone.neg_width
    wxy.x=wxy.x+owp.zone.neg_width
    el=Inf
    wl=Inf
    for nd in owp.zone.nodes
        l=lof_pnt2pnt_dist(nd.xy,exy)
        if (l<el)
            ePnt=deepcopy(nd)
            el=deepcopy(l)
        end
        l=lof_pnt2pnt_dist(nd.xy,wxy)
        if (l<wl)
            wPnt=deepcopy(nd)
            wl=deepcopy(l)
        end
    end
    return ePnt, wPnt
end

########################################## End HV System #####################################
########################################## Start MV System #####################################
#=
ocn=ocean
owpps=ocn.owpps
pcc=ocn.pccs[2]
=#
#**
function opt_mvOSSplacement(ocn,owpps,pcc)
    #owpps is an array of owpps ordered closest to farthest from the designated pcc
    #owpps_tbl=top_mvTopos(ocn.owpps)
    owpps_tbl=top_hvTopos(ocn.owpps)
    owpps_tbl=owpps_tbl[end:-1:1,end:-1:1]
    bus_dummies=Array{bus,1}()
    oss_system=circuit()
    oss_systems=ocn.circuits
    #frst=1
    for indx0=1:length(owpps_tbl[:,1])
    #indx0=14
        bus_dummies=bus[]
        mv_square=opt_mvConstraints(ocn,owpps_tbl[indx0,:])
        for (indx1,owpp) in enumerate(owpps_tbl[indx0,:])
            if (owpp==1)
                push!(bus_dummies,owpps[indx1])
            end
        end

        #=if (length(oss_systems)==savedIndx)
            println(sizeof(oss_systems))
            oss_systems=ppf_saveSystem(oss_systems,frst,indx0-1)
            frst=deepcopy(indx0)
            println(sizeof(oss_systems))
        end=#

        if (sum(owpps_tbl[indx0,:])>1)
            oss_system, discrete_dom, discrete_edges=opt_mvOssSystem2(mv_square,bus_dummies,pcc,ocn,owpps_tbl[indx0,:])
            rnk=floor(Int32,top_bin2dec(owpps_tbl[indx0,:]))
            #if (oss_system.cost<oss_systems[rnk].cost)
                oss_systems[rnk]=deepcopy(oss_system)
                #push!(oss_systems,deepcopy(oss_system))
                ocn.discretedom.nodes=deepcopy(discrete_dom)
                ocn.discretedom.edges=deepcopy(discrete_edges)
            #end
        else
        end
    end
    return oss_systems
end

#owp=bus_dummies[1]
#bn=owpps_tbl[indx0,:]
#**

function opt_str8Oss2Oss(crc,oss_system,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
    #Calculate length and type of cables
    dumb_pth=as_Astar(crc.base_owp.node,oss_system.osss_mog[1].node,ocn.discretedom.nodes)
    cb,xs=cstF_MvHvCbloss(dumb_pth.G_cost,crc.base_owp.mva,crc.base_owp.wnd,ocn.finance,oss_system.pcc_cbls[1].elec.volt,ocn.sys)

    #store cable data
    currentNode=deepcopy(dumb_pth)
    goalNode=deepcopy(dumb_pth.goal)
    if length(cb) == 1
        while currentNode.num != goalNode
            push!(cb[1].pth,deepcopy(currentNode))
            currentNode=currentNode.parent
        end
        push!(cb[1].pth,deepcopy(currentNode))
        cb[1].pth=reverse!(cb[1].pth)
        push!(oss_system.owp_MVcbls,deepcopy(cb[1]))
        #wind/power
        push!(pMv,crc.oss_mva)
        push!(wMv,crc.oss_wind)
    elseif length(cb) == 2
        #Cable
        push!(cb[2].pth,deepcopy(currentNode))
        while currentNode.parent.num != goalNode
            push!(cb[2].pth,deepcopy(currentNode.parent))
            currentNode=currentNode.parent
        end
        cb[2].pth=reverse!(cb[2].pth)
        push!(cb[1].pth,deepcopy(currentNode.parent))
        push!(cb[1].pth,deepcopy(currentNode))
        push!(oss_system.owp_MVcbls,deepcopy(cb[1]))
        push!(oss_system.owp_HVcbls,deepcopy(cb[2]))

        #OSS
        ossmv=bus()
        ossmv.node=currentNode.parent
        ossmv.base_cost=ocn.finance.FC_bld
        push!(ossmv.xfmrs,deepcopy(xs[1]))
        push!(oss_system.osss_owp,ossmv)

        #wind/power
        #store wind for mog transformer calcs
        if (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt==oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
            push!(pHv,crc.oss_mva)
            push!(wHv,crc.oss_wind)
        elseif (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt==132)
            push!(p132,crc.oss_mva)
            push!(w132,crc.oss_wind)
        elseif (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt==220)
            push!(p220,crc.oss_mva)
            push!(w220,crc.oss_wind)
        elseif (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt==400)
            push!(p400,crc.oss_mva)
            push!(w400,crc.oss_wind)
        else
        end
    end

    return oss_system,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400
end
#owp=bus_dummies[1]
#bn=owpps_tbl_hv[indx0,:]
#**
function opt_str8Connect(owp,pcc,ocn,bn)
    circ=circuit()
    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.pcc=deepcopy(pcc)
    circ.base_owp=owp
    push!(circ.owpps,owp)

    #Calculate length and type of cables
    dumb_pth=as_Astar(owp.node,pcc.node,ocn.discretedom.nodes)
    cs=cstF_MvHvCblpcc(dumb_pth.G_cost,owp.mva,owp.wnd,ocn.finance,pcc,ocn.sys)

    #store cable data
    currentNode=dumb_pth
    goalNode=dumb_pth.goal
    if length(cs) == 2
        while currentNode.parent.num != goalNode
            push!(cs[1].pth,deepcopy(currentNode))
            currentNode=currentNode.parent
        end
        cs[1].pth=reverse!(cs[1].pth)
        push!(circ.owp_MVcbls,deepcopy(cs[1]))
        #PCC
        push!(circ.pcc.xfmrs,deepcopy(cs[2]))
    elseif length(cs) == 4
        #Cable
        push!(cs[2].pth,deepcopy(currentNode))
        while currentNode.parent.num != goalNode
            push!(cs[2].pth,deepcopy(currentNode.parent))
            currentNode=currentNode.parent
        end
        cs[2].pth=reverse!(cs[2].pth)
        push!(cs[1].pth,deepcopy(currentNode.parent))
        push!(cs[1].pth,deepcopy(currentNode))
        push!(circ.owp_MVcbls,deepcopy(cs[1]))
        push!(circ.pcc_cbls,deepcopy(cs[2]))

        #OSS
        ossmv=bus()
        ossmv.node=currentNode
        ossmv.base_cost=ocn.finance.FC_bld
        ossmv.mva=deepcopy(owp.mva)
        push!(circ.osss_mog,ossmv)
        push!(circ.osss_mog[1].xfmrs,deepcopy(cs[3]))
        #PCC
        push!(circ.pcc.xfmrs,deepcopy(cs[4]))
    end
    circ.oss_mva=deepcopy(owp.mva)
    circ.oss_wind=deepcopy(owp.wnd)
    opt_ttlMvCirc(circ)
    return circ
end

########################################## End MV System #####################################
#**
function opt_mogXfmrs(mog,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ks,kV)
    wo=wind()
    wo.pu=zeros(Float32,8759)
    wo.ce=zeros(Float32,8759)
    wo.delta=0
    wo.lf=0
    if length(pMv) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(pMv),opt_Wsum(deepcopy(wo),wMv),ks,66,kV))
    end
    if length(p132) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(p132),opt_Wsum(deepcopy(wo),w132),ks,132,kV))
    end
    if length(p220) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(p220),opt_Wsum(deepcopy(wo),w220),ks,220,kV))
    end
    if length(p400) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(p400),opt_Wsum(deepcopy(wo),w400),ks,400,kV))
    end
    return mog
end

#**
function opt_Wsum(wo,wA)
    for w in wA
        wo.pu=(wo.pu.+w.pu)
        wo.ce=(wo.ce.+w.ce)
        wo.delta=(wo.delta+w.delta)
        wo.lf=(wo.lf+w.lf)
    end
    wo.pu=wo.pu./length(wA)
    wo.ce=wo.ce./length(wA)
    wo.delta=wo.delta/length(wA)
    wo.lf=wo.lf/length(wA)
    return wo
end

#**
function opt_InitPW()
    pSum_mv=Float32[]
    pSum_hv=Float32[]
    pSum_hv132=Float32[]
    pSum_hv220=Float32[]
    pSum_hv400=Float32[]

    wind_sumMv=wind[]
    wind_sumHv=wind[]
    wind_sum132=wind[]
    wind_sum220=wind[]
    wind_sum400=wind[]
    return pSum_mv,pSum_hv,pSum_hv132,pSum_hv220,pSum_hv400,wind_sumMv,wind_sumHv,wind_sum132,wind_sum220,wind_sum400
end

#**
function opt_findClosestNode(pos_xy,ndes)
    lwr_node=node()
    lwr_node_mag=Inf
    for nd in (ndes)
        tl_mag=lof_pnt2pnt_dist(nd.xy,pos_xy)
        if (tl_mag<lwr_node_mag)
            lwr_node_mag=deepcopy(tl_mag)
            lwr_node=deepcopy(nd)
        end
    end
    return lwr_node
end

#**
function opt_ttlMvCirc(circ)
    circ.cost=0
    for cb in circ.pcc_cbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for cb in circ.oss2oss_cbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for cb in circ.owp_MVcbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for cb in circ.owp_HVcbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for oss in circ.osss_owp
        circ.cost=circ.cost+oss.base_cost
        for xf in oss.xfmrs
            circ.cost=circ.cost+xf.costs.ttl
        end
    end

    for oss in circ.osss_mog
        circ.cost=circ.cost+oss.base_cost
        for xf in oss.xfmrs
            circ.cost=circ.cost+xf.costs.ttl
        end
    end

    for xf in circ.pcc.xfmrs
        circ.cost=circ.cost+xf.costs.ttl
    end

end
#oss_node=opNode
#**
function opt_ifIn2Out(oss_node,ocn,wp)
    bs,a=lof_test4pnt(oss_node.xy.x,oss_node.xy.y,ocn)
    dummy_dist=Float32
    bsf_dist=Inf
    bsf_indx=1
    if bs==false
        verln=lof_lineDirection(wp.node,oss_node)
        vect=lof_getStr8line(wp.node,oss_node)
        for (indx,pnt) in enumerate(a.nodes)
            if (((pnt.xy.x >= wp.node.xy.x && pnt.xy.x <= oss_node.xy.x)||(pnt.xy.x <= wp.node.xy.x && pnt.xy.x >= oss_node.xy.x))&&((pnt.xy.y >= wp.node.xy.y && pnt.xy.y <= oss_node.xy.y)||(pnt.xy.y <= wp.node.xy.y && pnt.xy.y >= oss_node.xy.y)))
                dummy_node=node()
                if (verln==true)
                    dummy_node.xy.x=pnt.xy.y*vect.m_findx+vect.b_findx
                    dummy_node.xy.y=pnt.xy.y
                    dummy_dist=lof_pnt2pnt_dist(oss_node.xy,pnt.xy)
                    if dummy_dist<bsf_dist
                        bsf_dist=deepcopy(dummy_dist)
                        bsf_indx=deepcopy(indx)
                    end
                else
                    dummy_node.xy.y=pnt.xy.x*vect.m_findy+vect.b_findy
                    dummy_node.xy.x=pnt.xy.x
                    dummy_dist=lof_pnt2pnt_dist(oss_node.xy,pnt.xy)
                    if dummy_dist<bsf_dist
                        bsf_dist=deepcopy(dummy_dist)
                        bsf_indx=deepcopy(indx)
                    end
                end
            end
        end
        oss_node=a.nodes[bsf_indx]
    else
    end
    return bs,oss_node
end
#all_opps=ocn.owpps
#**
function opt_buses2nodes4MVoptLocal(mv_square,all_opps,buses)
    target_owpp=1
    owp=deepcopy(buses[1])
    for (indx,bs) in enumerate(all_opps)
        if ((bs.node.xy.y<mv_square.ymn) || (bs.node.xy.x < mv_square.xmn) || (bs.num < owp.num))
            target_owpp=deepcopy(indx+1)
        end
    end
    owp=deepcopy(all_opps[target_owpp])
    return owp
end

#MV constraints
#**
#owpps=owpps_tbl[indx0,:]

function opt_mvConstraints(ocn, owpps)
    #last_owp=findlast(x->x==1,owpps)
    sqr_bnd=square()
    valid=false
    if (sum(owpps)>1)
        first_owp=findfirst(x->x==1,owpps)
        sqr_bnd.ymn=ocn.owpps[first_owp].node.xy.y-ocn.owpps[first_owp].mv_zone.neg_height
        sqr_bnd.ymx=ocn.owpps[first_owp].node.xy.y+ocn.owpps[first_owp].mv_zone.pos_height
        sqr_bnd.xmn=ocn.owpps[first_owp].node.xy.x-ocn.owpps[first_owp].mv_zone.neg_width
        sqr_bnd.xmx=ocn.owpps[first_owp].node.xy.x+ocn.owpps[first_owp].mv_zone.pos_width
        for (i,owp) in enumerate(owpps[first_owp+1:length(owpps)])
            if (owp==1)
                ymn=ocn.owpps[first_owp+i].node.xy.y-ocn.owpps[first_owp+i].mv_zone.neg_height
                ymx=ocn.owpps[first_owp+i].node.xy.y+ocn.owpps[first_owp+i].mv_zone.pos_height
                xmn=ocn.owpps[first_owp+i].node.xy.x-ocn.owpps[first_owp+i].mv_zone.neg_width
                xmx=ocn.owpps[first_owp+i].node.xy.x+ocn.owpps[first_owp+i].mv_zone.pos_width
                if (xmx>=sqr_bnd.xmn && xmn<=sqr_bnd.xmx && ymx>=sqr_bnd.ymn && ymn<=sqr_bnd.ymx)
                    valid=true
                end
            end
        end
    else
        valid=false
    end
    return valid,sqr_bnd
end
#=
function opt_mvConstraints_old(ocn, owpps)
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
end=#

#OWWP constraints
#**
function opt_owppConstraints(ocn)
    #@NLconstraint(m, ((x-x1)/a)^2+((y-y1)/b)^2 > 1)
    for owpp in ocn.owpps
        ellipse_con=ellipse()
        ellipse_con.y0=owpp.node.xy.y-(owpp.zone.neg_height-owpp.zone.pos_height)/2
        ellipse_con.x0=owpp.node.xy.x-(owpp.zone.neg_width-owpp.zone.pos_width)/2
        ellipse_con.ry=owpp.node.xy.y+owpp.zone.pos_height-ellipse_con.y0
        ellipse_con.rx=owpp.node.xy.x+owpp.zone.pos_width-ellipse_con.x0
        ellipse_con.alpha=0.0
        push!(ocn.constrain.ellipses,deepcopy(ellipse_con))
    end
end


#OWWP constraints
#**
function opt_nogoConstraints(ocn)
    #@NLconstraint(m, ((x-x1)/a)^2+((y-y1)/b)^2 > 1)
    for nogo in ocn.nogos
        ellipse_con=ellipse()
        simplex=opt_makeHalfSpace(nogo)
        dual_surf=findLargestEllipse(simplex)
        r_x,r_y,cntr,offset=form_EllipseConstraint(dual_surf)
        ellipse_con.y0=cntr.y
        ellipse_con.x0=cntr.x
        ellipse_con.ry=r_y
        ellipse_con.rx=r_x
        ellipse_con.alpha=offset
        push!(ocn.constrain.ellipses,deepcopy(ellipse_con))
    end
end

function optSysCompare(MVC,HVC)
    for (i,hc) in enumerate(HVC)
        if (MVC[i].cost < hc.cost)
            HVC[i]=deepcopy(MVC[i])
        end
    end
    return HVC
end

################################## depricated ##################################

#NOGO constraints - original
#=function opt_nogoConstraints(ocn)
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
end=#
#=
ng=ocean.nogos[1]
=#
#Finds closest edge
#=
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
=#
#=
ocean.owpps[2].node
mean=opt_findShapeCentre(ocean.owpps[2].zone)
=#
#=
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
=#
#=
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
=#
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
#=
function opt_findTriAngCentre(p0, p1, p2)
    cntr=xy()
    cntr.x=(p0.x+p1.x+p2.x)/3
    cntr.y=(p0.y+p1.y+p2.y)/3
    return cntr
end
=#
#=
function opt_buildSystem_orig(oss_systems,rnk,ocn)
    oss_system=circuit()

    #take base information
    oss_system.binary=ocn.circuits[rnk].binary
    oss_system.decimal=ocn.circuits[rnk].decimal
    oss_system.owpps=ocn.circuits[rnk].owpps
    oss_system.base_owp=ocn.circuits[rnk].base_owp
    oss_system.oss_wind=ocn.circuits[rnk].oss_wind
    oss_system.oss_mva=ocn.circuits[rnk].oss_mva

    #Take base OWP,PCC connections and MOG without transformers
    push!(oss_system.owp_MVcbls,deepcopy(ocn.circuits[rnk].owp_MVcbls[1]))
    oss_system.pcc=deepcopy(ocn.circuits[rnk].pcc)
    push!(oss_system.osss_mog,deepcopy(ocn.circuits[rnk].osss_mog[1]))
    oss_system.osss_mog[1].xfmrs=xfo[]
    push!(oss_system.pcc_cbls,deepcopy(ocn.circuits[rnk].pcc_cbls[1]))
    if (oss_system.owp_MVcbls[1].pth[length(oss_system.owp_MVcbls[1].pth)].num != oss_system.osss_mog[1].node.num)
        push!(oss_system.owp_HVcbls,deepcopy(ocn.circuits[rnk].owp_HVcbls[1]))
        push!(oss_system.osss_owp,deepcopy(ocn.circuits[rnk].osss_owp[1]))
    else
    end

    #store wind for transformer calculation at MOG
    pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_InitPW()
    if (oss_system.owp_MVcbls[1].pth[length(oss_system.owp_MVcbls[1].pth)].num == oss_system.osss_mog[1].node.num)
        push!(pMv,oss_systems[1].oss_mva)
        push!(wMv,oss_systems[1].oss_wind)
    elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
        push!(pHv,oss_systems[1].oss_mva)
        push!(wHv,oss_systems[1].oss_wind)
    elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==132)
        push!(p132,oss_systems[1].oss_mva)
        push!(w132,oss_systems[1].oss_wind)
    elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==220)
        push!(p220,oss_systems[1].oss_mva)
        push!(w220,oss_systems[1].oss_wind)
    elseif (oss_systems[1].pcc_cbls[length(oss_systems[1].pcc_cbls)].elec.volt==400)
        push!(p400,oss_systems[1].oss_mva)
        push!(w400,oss_systems[1].oss_wind)
    else
    end

    for crc in oss_systems[2:length(oss_systems)]
        #copy base cables and transformers
        if (length(crc.owpps)==1)
            oss_system,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_str8Oss2Oss(crc,oss_system,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
        else
            for cb in crc.owp_MVcbls
                push!(oss_system.owp_MVcbls,cb)
            end

            for cb in crc.owp_HVcbls
                push!(oss_system.owp_HVcbls,cb)
            end

            for cb in crc.oss2oss_cbls
                push!(oss_system.oss2oss_cbls,cb)
            end

            for xm in crc.osss_owp
                push!(oss_system.osss_owp,xm)
            end

            for xm in crc.osss_mog
                push!(oss_system.osss_mog,xm)
            end

            #find oss to mog connections
            pth=deepcopy(as_Astar(crc.osss_mog[1].node,oss_system.osss_mog[1].node,ocn.discretedom.nodes))
            oss2ossCbl=cstF_HvCblo2o(pth.G_cost,crc.oss_mva,crc.pcc_cbls[length(crc.pcc_cbls)].elec.volt,crc.oss_wind,ocn.finance)

            currentNode=deepcopy(pth)
            goalNode=deepcopy(pth.goal)
            while currentNode.num != goalNode
                push!(oss2ossCbl.pth,deepcopy(currentNode))
                currentNode=currentNode.parent
            end
            push!(oss2ossCbl.pth,deepcopy(currentNode))
            oss2ossCbl.pth=reverse!(oss2ossCbl.pth)
            push!(oss_system.oss2oss_cbls,deepcopy(oss2ossCbl))

            #store wind for mog transformer calcs
            if (crc.pcc_cbls[length(crc.pcc_cbls)].elec.volt==oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
                push!(pHv,crc.oss_mva)
                push!(wHv,crc.oss_wind)
            elseif (crc.pcc_cbls[length(crc.pcc_cbls)].elec.volt==132)
                push!(p132,crc.oss_mva)
                push!(w132,crc.oss_wind)
            elseif (crc.pcc_cbls[length(crc.pcc_cbls)].elec.volt==220)
                push!(p220,crc.oss_mva)
                push!(w220,crc.oss_wind)
            elseif (crc.pcc_cbls[length(crc.pcc_cbls)].elec.volt==400)
                push!(p400,crc.oss_mva)
                push!(w400,crc.oss_wind)
            else
            end
        end
    end

    oss_system.osss_mog[1]=opt_mogXfmrs(oss_system.osss_mog[1],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
    opt_ttlMvCirc(oss_system)
    println(string(oss_system.cost)*" - "*string(ocn.circuits[rnk].cost))
    if (oss_system.cost < ocn.circuits[rnk].cost)
        ocn.circuits[rnk]=deepcopy(oss_system)
    end
end
=#

#=
function opt_hvOssSystem(buses,pcc,ocn,bn,side)
    owp=buses[1]
    circ,domain_oss,edges_oss=opt_mvHvMain(buses,pcc,ocn,bn,side,owp)
    return circ,domain_oss,edges_oss
end=#

#=
function opt_mvOssSystem(mv_square,buses,pcc,ocn,bn,side)
    owp=opt_buses2nodes4MVoptLocal(mv_square,ocn.owpps,buses)
    circ,domain_oss,edges_oss=opt_mvHvMain(buses,pcc,ocn,bn,side,owp)
    #ocn.discretedom.nodes=deepcopy(domain_oss)
    #ocn.discretedom.edges=deepcopy(edges_oss)
    return circ, domain_oss, edges_oss
end
=#

#=
bn=owpps_tbl[indx0,:]
buses=bus_dummies
=#

#=
function opt_mvHvMain(buses,pcc,ocn,bn,side,owp)
    circ=circuit()
    circ.base_owp=owp
    oss_node,domain_oss,edges_oss,circ.oss_mva,circ.oss_wind,as_pths=opt_mvhvOss1stLocal(owp,buses,pcc,ocn,side)

    #Find PCC connection cable  and PCC transformer
    println("mva: "*string(circ.oss_mva))
    cbl_xfo=cstF_HvCblallKvo2p(as_pths[1].G_cost,circ.oss_mva,circ.oss_wind,ocn.finance,pcc)

    #Set cable path
    currentNode=deepcopy(as_pths[1])
    goalNode=as_pths[1].goal
    push!(cbl_xfo[1].pth,deepcopy(currentNode))
    while currentNode.num != goalNode
        push!(cbl_xfo[1].pth,deepcopy(currentNode.parent))
        currentNode=currentNode.parent
    end
    push!(circ.pcc_cbls,cbl_xfo[1])

    #pcc transformer
    circ.pcc=deepcopy(pcc)
    push!(circ.pcc.xfmrs,cbl_xfo[2])

    #MV cable sizes and paths
    pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_InitPW()
    for (i,owp) in enumerate(buses)
        cb,xs=cstF_MvHvCbloss(as_pths[i+1].G_cost,owp.mva,owp.wnd,ocn.finance,circ.pcc_cbls[1].elec.volt,ocn.sys)
        if (length(cb)==1)
            currentNode=deepcopy(as_pths[i+1])
            goalNode=currentNode.goal
            push!(cb[1].pth,deepcopy(currentNode))
            while (currentNode.num != goalNode)
                push!(cb[1].pth,deepcopy(currentNode.parent))
                currentNode=currentNode.parent
            end
            cb[1].pth=reverse!(cb[1].pth)
            push!(circ.owp_MVcbls,cb[1])
            #store winds and powers
            push!(pMv,owp.mva)
            push!(wMv,owp.wnd)
        elseif (length(cb)==2)
            #Cable
            currentNode=deepcopy(as_pths[i+1])
            goalNode=currentNode.goal
            push!(cb[2].pth,deepcopy(currentNode))
            while currentNode.parent.num != goalNode
                push!(cb[2].pth,deepcopy(currentNode.parent))
                currentNode=currentNode.parent
            end
            cb[2].pth=reverse!(cb[2].pth)
            push!(cb[1].pth,deepcopy(currentNode.parent))
            push!(cb[1].pth,deepcopy(currentNode))
            push!(circ.owp_MVcbls,deepcopy(cb[1]))
            push!(circ.owp_HVcbls,deepcopy(cb[2]))

            #OSS
            ossmv=bus()
            ossmv.node=currentNode
            ossmv.base_cost=ocn.finance.FC_bld
            push!(circ.osss_owp,ossmv)
            push!(circ.osss_owp[length(circ.osss_owp)].xfmrs,deepcopy(xs[1]))

            #store winds and powers
            if (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)
                push!(pHv,owp.mva)
                push!(wHv,owp.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 132)
                push!(p132,owp.mva)
                push!(w132,owp.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 220)
                push!(p220,owp.mva)
                push!(w220,owp.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 400)
                push!(p400,owp.mva)
                push!(w400,owp.wnd)
            end
        end
    end

    #mog transformer
    ossMOG=bus()
    ossMOG.mva=circ.oss_mva
    ossMOG.wnd=circ.oss_wind
    ossMOG.node=as_pths[1]
    ossMOG.base_cost=ocn.finance.FC_bld
    push!(circ.osss_mog,ossMOG)
    circ.osss_mog[length(circ.osss_mog)]=opt_mogXfmrs(circ.osss_mog[length(circ.osss_mog)],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)

    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.owpps=buses
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    return circ,domain_oss,edges_oss
end
=#

#=
function opt_adjustPath(owp,ocn,xys,buses,pcc)
    paths=node[]
    opt_MakeConstraints(xys,owp,ocn.owpps,pcc)
    oss_node=opt_OssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
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
    pcc_xy=xy()
    #if (paths[1].parent.num != pcc.node.num)
    #    pcc_xy=paths[1].parent.xy
    #else
        pcc_xy=pcc.node.xy
    #end
    push!(xys,pcc_xy)
    for path=2:length(paths)
        pathsParent=paths[path]
        while pathsParent.parent.num != buses[path-1].node.num
            pathsParent=pathsParent.parent
        end
        push!(xys,deepcopy(pathsParent.xy))
    end

    #re calculate the path
    paths=node[]
    #error
    oss_node=opt_OssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    push!(paths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float32,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float32,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0
    for (i,owpp) in enumerate(buses)
        wind_sum.pu=(wind_sum.pu.+deepcopy(owpp.wnd.pu))
        wind_sum.ce=(wind_sum.ce.+deepcopy(owpp.wnd.ce))
        wind_sum.delta=(wind_sum.delta+deepcopy(owpp.wnd.delta))
        wind_sum.lf=(wind_sum.lf+deepcopy(owpp.wnd.lf))
        power_sum=power_sum+deepcopy(owpp.mva)
        push!(paths,deepcopy(as_Astar(domain_oss[owpp.node.num],oss_node,domain_oss)))
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)
    #pick first node in each owpp path
    xys=xy[]
    for path=1:length(paths)
        if (paths[path].parent.num != paths[path].parent.goal)
            push!(xys,deepcopy(paths[path].parent.xy))
        else
            push!(xys,deepcopy(paths[path].xy))
        end
    end

    #re calculate the path
    paths=node[]
    oss_node=opt_OssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    return domain_oss,oss_node,power_sum,wind_sum,edges_oss
end
=#
#=
function opt_mvhvOss1stLocal_c(owp,buses,pcc,ocn)
    #find PCC major axis
    ymajorPCC=true
    if (abs(pcc.node.xy.y-owp.node.xy.y)>=abs(pcc.node.xy.x-owp.node.xy.x))
        ymajorPCC=true
    else
        ymajorPCC=false
    end

    ymajorOWPP=true
    if (abs(buses[1].node.xy.y-buses[length(buses)].node.xy.y)>=abs(buses[1].node.xy.x-buses[length(buses)].node.xy.x))
        ymajorOWPP=true
    else
        ymajorOWPP=false
    end


    #find min and max connection points
    cnt_xy=xy()
    if (ymajorPCC==true)
        if (pcc.node.xy.y>=owp.node.xy.y)
            cnt_xy.y=owp.node.xy.y+owp.zone.pos_height
        else
            cnt_xy.y=owp.node.xy.y-owp.zone.neg_height
        end
        if (owp.node.xy.x-owp.zone.neg_width<pcc.node.xy.x && pcc.node.xy.x<owp.node.xy.x+owp.zone.pos_width)
            cnt_xy.x=pcc.node.xy.x
        else
            cnt_xy.x=owp.node.xy.x
        end
    else
        if (pcc.node.xy.x>=owp.node.xy.y)
            cnt_xy.x=owp.node.xy.x+owp.zone.pos_width
        else
            cnt_xy.x=owp.node.xy.x-owp.zone.neg_width
        end
        if (owp.node.xy.y-owp.zone.neg_height<pcc.node.xy.y && pcc.node.xy.y<owp.node.xy.y+owp.zone.pos_height)
            cnt_xy.y=pcc.node.xy.y
        else
            cnt_xy.y=owp.node.xy.y
        end
    end

    #find nodes
    cnt_node=node()
    cnt_node_mag=Inf

    for nd in (owp.zone.nodes)
        tc_mag=lof_pnt2pnt_dist(nd.xy,cnt_xy)
        if (tc_mag<cnt_node_mag)
            cnt_node_mag=deepcopy(tc_mag)
            cnt_node=deepcopy(nd)
        end
    end


    #each owpp path
    cnt_paths=node[]
    cnt_lths=0
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float32,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float32,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0
    push!(cnt_paths,deepcopy(as_Astar(pcc.node,cnt_node,ocn.discretedom.nodes)))
    for opp in buses
        push!(cnt_paths,deepcopy(as_Astar(opp.node,cnt_node,ocn.discretedom.nodes)))
        wind_sum.pu=(wind_sum.pu.+opp.wnd.pu)
        wind_sum.ce=(wind_sum.ce.+opp.wnd.ce)
        wind_sum.delta=(wind_sum.delta+opp.wnd.delta)
        wind_sum.lf=(wind_sum.lf+opp.wnd.lf)
        power_sum=(power_sum+opp.mva)
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)



        best_path=cnt_paths
        best_node=deepcopy(best_path[1])

        ojama_paths=node[]
        ojama_xys=xy[]
        for nd in ocn.owpps[owp.num+1:buses[length(buses)].num]
            dummy_path=as_Astar(nd,best_node,ocn.discretedom.nodes)
            push!(ojama_paths,deepcopy(dummy_path))
            goal=deepcopy(dummy_path.goal)
            while (dummy_path.parent.num != goal)
                dummy_path=dummy_path.parent
            end
            push!(ojama_xys,deepcopy(dummy_path.xy))
        end
        #adjust best connection point
        path2adjust=deepcopy(best_node)
        if (ymajorOWPP==true)
            room2adjust=deepcopy(best_node.xy.x)
            if (best_node.xy.x<=owp.node.xy.x)
                goal=deepcopy(path2adjust.goal)
                while path2adjust.num != goal
                    path2adjust=path2adjust.parent
                    if (path2adjust.xy.x<room2adjust)
                        room2adjust=deepcopy(path2adjust.xy.x)
                    end
                end
                for xy in ojama_xys
                    if (xy.x<best_node.xy.x && xy.x>=room2adjust)
                        best_node.xy.x=deepcopy(best_node.xy.x)
                    end
                end
            else
                goal=deepcopy(path2adjust.goal)
                while path2adjust.num != goal
                    path2adjust=path2adjust.parent
                    if (path2adjust.xy.x>room2adjust)
                        room2adjust=deepcopy(path2adjust.xy.x)
                    end
                end
                for xy in ojama_xys
                    if (xy.x>best_node.xy.x && xy.x<=room2adjust)
                        best_node.xy.x=deepcopy(best_node.xy.x)
                    end
                end
            end
        else
            room2adjust=deepcopy(best_node.xy.y)
            goal=deepcopy(path2adjust.goal)
            if (best_node.xy.y<=owp.node.xy.y)
                while path2adjust.num != goal
                    path2adjust=path2adjust.parent
                    if (path2adjust.xy.y<room2adjust)
                        room2adjust=deepcopy(path2adjust.xy.y)
                    end
                end
                for xy in ojama_xys
                    if (xy.y<best_node.xy.y && xy.y>=room2adjust)
                        best_node.xy.y=deepcopy(best_node.xy.y)
                    end
                end
            else
                while path2adjust.num != goal
                    path2adjust=path2adjust.parent
                    if (path2adjust.xy.y>room2adjust)
                        room2adjust=deepcopy(path2adjust.xy.y)
                    end
                end
                for xy in ojama_xys
                    if (xy.y>best_node.xy.y && xy.y<=room2adjust)
                        best_node.xy.y=deepcopy(best_node.xy.y)
                    end
                end
            end
        end

        domain_oss=ocn.discretedom.nodes
        edges_oss=ocn.discretedom.edges
        bs,best_node=opt_ifIn2Out(best_node,ocn)
        if bs==true
            best_node.num=deepcopy(length(domain_oss)+1)
            push!(domain_oss,best_node)
            lof_edgeifyOss(best_node,ocn,domain_oss,edges_oss)
        end

        if ((best_node.xy.x != best_path[1].xy.x) || (best_node.xy.y != best_path[1].xy.y))
            best_path=node[]
            push!(best_path,as_Astar(pcc.node,best_node,ocn.discretedom.nodes))
            for opp in buses
                push!(best_path,as_Astar(opp.node,best_node,ocn.discretedom.nodes))
            end
        end
        return best_node,domain_oss,edges_oss,power_sum,wind_sum,best_path
end
=#
#=
function opt_mvhvOss1stLocal_h(owp,buses,pcc,ocn)
    #find PCC major axis
    ymajorPCC=true
    if (abs(pcc.node.xy.y-owp.node.xy.y)>=abs(pcc.node.xy.x-owp.node.xy.x))
        ymajorPCC=true
    else
        ymajorPCC=false
    end

    ymajorOWPP=true
    if (abs(buses[1].node.xy.y-buses[length(buses)].node.xy.y)>=abs(buses[1].node.xy.x-buses[length(buses)].node.xy.x))
        ymajorOWPP=true
    else
        ymajorOWPP=false
    end

    #find min and max connection points
    upr_xy=xy()
    if (ymajorPCC==true)
        if (pcc.node.xy.y>=owp.node.xy.y)
            upr_xy.y=owp.node.xy.y+owp.zone.pos_height
        else
            upr_xy.y=owp.node.xy.y-owp.zone.neg_height
        end
        upr_xy.x=owp.node.xy.x+owp.zone.pos_width
    else
        if (pcc.node.xy.x>=owp.node.xy.y)
            upr_xy.x=owp.node.xy.x+owp.zone.pos_width
        else
            upr_xy.x=owp.node.xy.x-owp.zone.neg_width
        end
        upr_xy.y=owp.node.xy.y+owp.zone.pos_height

        if (owp.node.xy.y-owp.zone.neg_height<pcc.node.xy.y && pcc.node.xy.y<owp.node.xy.y+owp.zone.pos_height)
            cnt_xy.y=pcc.node.xy.y
        else
            cnt_xy.y=owp.node.xy.y
        end
    end

    #find nodes
    upr_node=node()
    upr_node_mag=Inf

    for nd in (owp.zone.nodes)
        tu_mag=lof_pnt2pnt_dist(nd.xy,upr_xy)
        if (tu_mag<upr_node_mag)
            upr_node_mag=deepcopy(tu_mag)
            upr_node=deepcopy(nd)
        end
    end


    #each owpp path
    upr_paths=node[]
    upr_lths=0
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float32,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float32,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0

    push!(upr_paths,deepcopy(as_Astar(pcc.node,upr_node,ocn.discretedom.nodes)))
    for opp in buses
        push!(upr_paths,deepcopy(as_Astar(opp.node,upr_node,ocn.discretedom.nodes)))
        wind_sum.pu=(wind_sum.pu.+opp.wnd.pu)
        wind_sum.ce=(wind_sum.ce.+opp.wnd.ce)
        wind_sum.delta=(wind_sum.delta+opp.wnd.delta)
        wind_sum.lf=(wind_sum.lf+opp.wnd.lf)
        power_sum=(power_sum+opp.mva)
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)

    best_path=upr_paths
    best_node=deepcopy(best_path[1])

    ojama_paths=node[]
    ojama_xys=xy[]
    for nd in ocn.owpps[owp.num+1:buses[length(buses)].num]
        dummy_path=as_Astar(nd,best_node,ocn.discretedom.nodes)
        push!(ojama_paths,deepcopy(dummy_path))
        goal=deepcopy(dummy_path.goal)
        while (dummy_path.parent.num != goal)
            dummy_path=dummy_path.parent
        end
        push!(ojama_xys,deepcopy(dummy_path.xy))
    end
    #adjust best connection point
    path2adjust=deepcopy(best_node)
    if (ymajorOWPP==true)
        room2adjust=deepcopy(best_node.xy.x)
        if (best_node.xy.x<=owp.node.xy.x)
            goal=deepcopy(path2adjust.goal)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.x<room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.x)
                end
            end
            for xy in ojama_xys
                if (xy.x<best_node.xy.x && xy.x>=room2adjust)
                    best_node.xy.x=deepcopy(best_node.xy.x)
                end
            end
        else
            goal=deepcopy(path2adjust.goal)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.x>room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.x)
                end
            end
            for xy in ojama_xys
                if (xy.x>best_node.xy.x && xy.x<=room2adjust)
                    best_node.xy.x=deepcopy(best_node.xy.x)
                end
            end
        end
    else
        room2adjust=deepcopy(best_node.xy.y)
        goal=deepcopy(path2adjust.goal)
        if (best_node.xy.y<=owp.node.xy.y)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.y<room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.y)
                end
            end
            for xy in ojama_xys
                if (xy.y<best_node.xy.y && xy.y>=room2adjust)
                    best_node.xy.y=deepcopy(best_node.xy.y)
                end
            end
        else
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.y>room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.y)
                end
            end
            for xy in ojama_xys
                if (xy.y>best_node.xy.y && xy.y<=room2adjust)
                    best_node.xy.y=deepcopy(best_node.xy.y)
                end
            end
        end
    end

    domain_oss=ocn.discretedom.nodes
    edges_oss=ocn.discretedom.edges
    bs,best_node=opt_ifIn2Out(best_node,ocn)
    if bs==true
        best_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,best_node)
        lof_edgeifyOss(best_node,ocn,domain_oss,edges_oss)
    end

    if ((best_node.xy.x != best_path[1].xy.x) || (best_node.xy.y != best_path[1].xy.y))
        best_path=node[]
        push!(best_path,as_Astar(pcc.node,best_node,ocn.discretedom.nodes))
        for opp in buses
            push!(best_path,as_Astar(opp.node,best_node,ocn.discretedom.nodes))
        end
    end
    return best_node,domain_oss,edges_oss,power_sum,wind_sum,best_path
end
=#
#=
function opt_OssOptimalLocal(xys,constrain,owp,num_owpp)

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
=#
#=
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
=#
#=
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
=#
#=
function opt_findCorner_l(ymajorPCC,pcc,owp)
    lwr_xy=xy()
    if (ymajorPCC==true)
        if (pcc.node.xy.y>=owp.node.xy.y)
            lwr_xy.y=owp.node.xy.y+owp.zone.pos_height
        else
            lwr_xy.y=owp.node.xy.y-owp.zone.neg_height
        end
        lwr_xy.x=owp.node.xy.x-owp.zone.neg_width
    else
        if (pcc.node.xy.x>=owp.node.xy.y)
            lwr_xy.x=owp.node.xy.x+owp.zone.pos_width
        else
            lwr_xy.x=owp.node.xy.x-owp.zone.neg_width
        end
        lwr_xy.y=owp.node.xy.y-owp.zone.neg_height
    end
    return lwr_xy
end
=#
#=
function opt_findCorner_h(ymajorPCC,pcc,owp)
    upr_xy=xy()
    if (ymajorPCC==true)
        if (pcc.node.xy.y>=owp.node.xy.y)
            upr_xy.y=owp.node.xy.y+owp.zone.pos_height
        else
            upr_xy.y=owp.node.xy.y-owp.zone.neg_height
        end
        upr_xy.x=owp.node.xy.x+owp.zone.pos_width
    else
        if (pcc.node.xy.x>=owp.node.xy.y)
            upr_xy.x=owp.node.xy.x+owp.zone.pos_width
        else
            upr_xy.x=owp.node.xy.x-owp.zone.neg_width
        end
        upr_xy.y=owp.node.xy.y+owp.zone.pos_height
    end

    return upr_xy
end
=#
#=
function opt_findCorner_c(ymajorPCC,pcc,owp)
    #find min and max connection points
    cnt_xy=xy()
    if (ymajorPCC==true)
        if (pcc.node.xy.y>=owp.node.xy.y)
            cnt_xy.y=owp.node.xy.y+owp.zone.pos_height
        else
            cnt_xy.y=owp.node.xy.y-owp.zone.neg_height
        end
        if (owp.node.xy.x-owp.zone.neg_width<pcc.node.xy.x && pcc.node.xy.x<owp.node.xy.x+owp.zone.pos_width)
            cnt_xy.x=pcc.node.xy.x
        else
            cnt_xy.x=owp.node.xy.x
        end
    else
        if (pcc.node.xy.x>=owp.node.xy.y)
            cnt_xy.x=owp.node.xy.x+owp.zone.pos_width
        else
            cnt_xy.x=owp.node.xy.x-owp.zone.neg_width
        end
        if (owp.node.xy.y-owp.zone.neg_height<pcc.node.xy.y && pcc.node.xy.y<owp.node.xy.y+owp.zone.pos_height)
            cnt_xy.y=pcc.node.xy.y
        else
            cnt_xy.y=owp.node.xy.y
        end
    end

    return cnt_xy
end
=#
#=
function opt_findWindPths(pcc,lwr_node,ocn,buses)
    #each owpp path
    lwr_paths=node[]
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float32,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float32,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0

    push!(lwr_paths,deepcopy(as_Astar(pcc.node,lwr_node,ocn.discretedom.nodes)))
    for opp in buses
        push!(lwr_paths,deepcopy(as_Astar(opp.node,lwr_node,ocn.discretedom.nodes)))
        wind_sum.pu=(wind_sum.pu.+opp.wnd.pu)
        wind_sum.ce=(wind_sum.ce.+opp.wnd.ce)
        wind_sum.delta=(wind_sum.delta+opp.wnd.delta)
        wind_sum.lf=(wind_sum.lf+opp.wnd.lf)
        power_sum=(power_sum+opp.mva)
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)
    return power_sum,wind_sum,lwr_paths
end
=#
#opt_node=opNode


#=
function opt_findOjamaXY(ocn,owp,buses,side)
    ojama_xys=xy[]
    upr_lim=buses[length(buses)].num
    lwr_lim=buses[1].num+1
    for nd in ocn.owpps[lwr_lim:upr_lim]
        dummy_xy=deepcopy(nd.node.xy)
        if (side=="low")
            dummy_xy.x=dummy_xy.x-nd.zone.neg_width
            push!(ojama_xys,deepcopy(dummy_xy))
        end
        if (side=="high")
            dummy_xy.x=dummy_xy.x+nd.zone.pos_width
            push!(ojama_xys,deepcopy(dummy_xy))
        end
    end
    return ojama_xys
end
=#
#=
function opt_adjustConnectionPnt(best_node,ymajorOWPP,owp,ojama_xys)
    path2adjust=deepcopy(best_node)
    if (ymajorOWPP==true)
        room2adjust=deepcopy(best_node.xy.x)
        if (best_node.xy.x<=owp.node.xy.x)
            goal=deepcopy(path2adjust.goal)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.x<room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.x)
                end
            end
            for xy in ojama_xys
                if (xy.x<best_node.xy.x && xy.x>=room2adjust)
                    best_node.xy.x=deepcopy(xy.x)
                end
            end
        else
            goal=deepcopy(path2adjust.goal)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.x>room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.x)
                end
            end
            for xy in ojama_xys
                if (xy.x>best_node.xy.x && xy.x<=room2adjust)
                    best_node.xy.x=deepcopy(xy.x)
                end
            end
        end
    else
        room2adjust=deepcopy(best_node.xy.y)
        goal=deepcopy(path2adjust.goal)
        if (best_node.xy.y<=owp.node.xy.y)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.y<room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.y)
                end
            end
            for xy in ojama_xys
                if (xy.y<best_node.xy.y && xy.y>=room2adjust)
                    best_node.xy.y=deepcopy(xy.y)
                end
            end
        else
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.y>room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.y)
                end
            end
            for xy in ojama_xys
                if (xy.y>best_node.xy.y && xy.y<=room2adjust)
                    best_node.xy.y=deepcopy(xy.y)
                end
            end
        end
    end
    return best_node
end
=#
#=
function opt_mvhvOss1stLocal(owp,buses,pcc,ocn,side)
    #find PCC major axis
    ymajorPCC=true
    if (abs(pcc.node.xy.y-owp.node.xy.y)>=abs(pcc.node.xy.x-owp.node.xy.x))
        ymajorPCC=true
    else
        ymajorPCC=false
    end

    ymajorOWPP=true
    if (abs(buses[1].node.xy.y-buses[length(buses)].node.xy.y)>=abs(buses[1].node.xy.x-buses[length(buses)].node.xy.x))
        ymajorOWPP=true
    else
        ymajorOWPP=false
    end

    #find min and max connection points
    if (side=="low")
        lwr_xy=opt_findCorner_l(ymajorPCC,pcc,owp)
    elseif (side=="high")
        lwr_xy=opt_findCorner_h(ymajorPCC,pcc,owp)
    else
        lwr_xy=opt_findCorner_c(ymajorPCC,pcc,owp)
    end
    #find nodes
    lwr_node=opt_findClosestNode(lwr_xy,owp.zone.nodes)

    power_sum,wind_sum,lwr_paths=opt_findWindPths(pcc,lwr_node,ocn,buses)

    best_path=lwr_paths
    best_node=deepcopy(best_path[1])
    if (side=="low" || side=="high")
        ojama_xys=opt_findOjamaXY(ocn,owp,buses,side)
        #adjust best connection point
        best_node=opt_adjustConnectionPnt(best_node,ymajorOWPP,owp,ojama_xys)
    end


    domain_oss=ocn.discretedom.nodes
    edges_oss=ocn.discretedom.edges
    bs,best_node=opt_ifIn2Out(best_node,ocn)
    if bs==true
        best_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,best_node)
        lof_edgeifyOss(best_node,ocn,domain_oss,edges_oss)
    end

    if ((best_node.xy.x != best_path[1].xy.x) || (best_node.xy.y != best_path[1].xy.y))
        best_path=node[]
        push!(best_path,deepcopy(as_Astar(pcc.node,best_node,ocn.discretedom.nodes)))
        for opp in buses
            push!(best_path,deepcopy(as_Astar(opp.node,best_node,ocn.discretedom.nodes)))
        end
    end
    return best_node,domain_oss,edges_oss,power_sum,wind_sum,best_path
end
=#

#=
ocn=ocean
indx0=1
indx1=indx0+1

crc0 = circs[indx0]
=#
#=
function opt_getX0(A,B,C,D,E)
    x0=(2*C*D-B*E)/(B^2-4*A*C)
    return x0
end

function opt_getY0(A,B,C,D,E)
    y0=(2*A*E-B*D)/(B^2-4*A*C)
    return y0
end

function opt_getRadiuss(A,B,C,D,E,F)
    a=(-sqrt(complex(2*(A*E^2+C*D^2-B*D*E+(B^2-4*A*C)*F)*((A+C)+sqrt((A-C)^2+B^2)))))/(B^2-4*A*C)
    b=(-sqrt(complex(2*(A*E^2+C*D^2-B*D*E+(B^2-4*A*C)*F)*((A+C)-sqrt((A-C)^2+B^2)))))/(B^2-4*A*C)
    return a,b
end
=#
