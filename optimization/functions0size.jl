#=
circ=ocean_mv.circuits[255]
co_ords=new_coords
ocn=ocean_mv
=#
function opt_reAdjust_oss(system,mv_rng,mv_cl)
    mog_xys=Array{Tuple{xy,Int64},1}()
    oss_xys=Array{Tuple{xy,Int64},1}()
    eps=10^-6
    for mog in system.osss_mog
        push!(mog_xys,(mog.node.xy,mog.node.num))
    end
    for os in system.osss_owp
        push!(oss_xys,(os.node.xy,os.node.num))
    end
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))

    @variables(m, begin
               MOG_y[i=1:length(mog_xys)]>=0, (start = mog_xys[i][1].y, base_name = "MOGY_"*string(mog_xys[i][2]))
               MOG_x[i=1:length(mog_xys)]>=0, (start = mog_xys[i][1].x, base_name = "MOGX_"*string(mog_xys[i][2]))
           end)

       @variables(m, begin
                  OSS_y[i=1:length(oss_xys)]>=0, (start = oss_xys[i][1].y, base_name = "OSSY_"*string(oss_xys[i][2]))
                  OSS_x[i=1:length(oss_xys)]>=0, (start = oss_xys[i][1].x, base_name = "OSSX_"*string(oss_xys[i][2]))
              end)

    connections=Array{Tuple{Any,Any,Any,Any,Float64},1}()
    for mv_c in system.owp_MVcbls
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==mv_c.pth[length(mv_c.pth)].num)
                push!(connections,(mv_c.pth[1].xy.x, mv_c.pth[1].xy.y, MOG_x[i], MOG_y[i], mv_c.costs.ttl/mv_c.length))
                @constraint(m,((connections[length(connections)][1]-connections[length(connections)][3])^2+(connections[length(connections)][2]-connections[length(connections)][4])^2) <= mv_rng^2)
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==mv_c.pth[length(mv_c.pth)].num)
                push!(connections,(mv_c.pth[1].xy.x, mv_c.pth[1].xy.y, OSS_x[i], OSS_y[i], mv_c.costs.ttl/mv_c.length))
                @constraint(m,((connections[length(connections)][1]-connections[length(connections)][3])^2+(connections[length(connections)][2]-connections[length(connections)][4])^2) <= mv_cl^2)
            end
        end
    end

    #PCC_connections=Array{Tuple{VariableRef,VariableRef,Float64,Float64,Float64},1}()
    for pcc_c in system.pcc_cbls
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==pcc_c.pth[1].num)
                push!(connections,(MOG_x[i], MOG_y[i], pcc_c.pth[length(pcc_c.pth)].xy.x, pcc_c.pth[length(pcc_c.pth)].xy.y, pcc_c.costs.ttl/pcc_c.length))
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==pcc_c.pth[length(pcc_c.pth)].num)
                push!(connections,(OSS_x[i], OSS_y[i], pcc_c.pth[length(pcc_c.pth)].xy.x, pcc_c.pth[length(pcc_c.pth)].xy.y, pcc_c.costs.ttl/pcc_c.length))
            end
        end
    end

    #HV_connections=Array{Tuple{VariableRef,VariableRef,VariableRef,VariableRef,Float64},1}()
    for hv_c in system.owp_HVcbls
        tail_x=VariableRef
        tail_y=VariableRef
        head_x=VariableRef
        head_y=VariableRef
        for (i,mg) in enumerate(mog_xys)

            if (mg[2]==hv_c.pth[length(hv_c.pth)].num)
                head_x=MOG_x[i]
                head_y=MOG_y[i]
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==hv_c.pth[1].num)
                tail_x=OSS_x[i]
                tail_y=OSS_y[i]
            end
        end
        push!(connections,(tail_x, tail_y, head_x, head_y, hv_c.costs.ttl/hv_c.length))
    end
    #o2o_connections=Array{Tuple{VariableRef,VariableRef,VariableRef,VariableRef,Float64},1}()
    for o2o_c in system.oss2oss_cbls
        tail_x=VariableRef
        tail_y=VariableRef
        head_x=VariableRef
        head_y=VariableRef
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==o2o_c.pth[1].num)
                tail_x=MOG_x[i]
                tail_y=MOG_y[i]
            end
            if (mg[2]==o2o_c.pth[length(o2o_c.pth)].num)
                head_x=MOG_x[i]
                head_y=MOG_y[i]
            end
        end

        push!(connections,(tail_x, tail_y, head_x, head_y, o2o_c.costs.ttl/o2o_c.length))
    end

    @NLobjective(m, Min,sum(sqrt((((connections[i][1]-connections[i][3])^2+(connections[i][2]-connections[i][4])^2)+eps))*connections[i][5] for i in 1:length(connections)))

    optimize!(m)
    mxy=(JuMP.value.(MOG_x),JuMP.value.(MOG_y),JuMP.value.(OSS_x),JuMP.value.(OSS_y))
    return mxy
end


function opt_readjust_circuits(ocn)
    mvrng=cstF_mVrng(5,ocn.owpps[1].mva,ocn.owpps[1].wnd,ocn.finance,ocn.sys)
    for i=1:length(ocn.circuits)
        print(string(i)*" Initial: "*string(ocn.circuits[i].cost))
        new_coords=opt_reAdjust_oss(ocn.circuits[i],mvrng,ocn.sys.mvCl)
        ocn.circuits[i]=opt_reAdjust_cbls(ocn.circuits[i],new_coords,ocn)
        println(" - Adjusted: "*string(ocn.circuits[i].cost))
    end
    return ocn
end


function opt_combineCblRoute_0(pth,oss2ossCbl,oss_system)
    push!(oss_system.oss2oss_cbls,deepcopy(oss2ossCbl))
    return oss2ossCbl,oss_system
end

function opt_compoundCbl_0(circ,oss_system,ocn)

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
    pth=deepcopy(lof_pnt2pnt_dist(circ.osss_mog[1].node.xy,oss_system.osss_mog[1].node.xy))
    oss2ossCbl=cstF_nextSizeDown(pth,circ.oss_mva,circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt,circ.oss_wind,ocn.finance,circ.pcc_cbls[length(circ.pcc_cbls)].size,circ.pcc_cbls[length(circ.pcc_cbls)].num)
    push!(oss2ossCbl.pth,circ.osss_mog[1].node)
    push!(oss2ossCbl.pth,oss_system.osss_mog[1].node)
    return oss2ossCbl,pth
end


function opt_str8Oss2Oss_0(crc,oss_system,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
    #Calculate length and type of cables
    dumb_pth=lof_pnt2pnt_dist(crc.base_owp.node.xy,oss_system.osss_mog[1].node.xy)
    cb,xs=cstF_MvHvCbloss(dumb_pth,crc.base_owp.mva,crc.base_owp.wnd,ocn.finance,oss_system.pcc_cbls[1].elec.volt,ocn.sys)

    #store cable data
    goalNode=deepcopy(oss_system.osss_mog[1].node)
    currentNode=deepcopy(crc.base_owp.node)
    if length(cb) == 1
        push!(cb[1].pth,deepcopy(currentNode))
        push!(cb[1].pth,deepcopy(goalNode))
        push!(oss_system.owp_MVcbls,deepcopy(cb[1]))
        #wind/power
        push!(pMv,crc.oss_mva)
        push!(wMv,crc.oss_wind)
    elseif length(cb) == 2
        verln=lof_lineDirection(currentNode,goalNode)
        vect=lof_getStr8line(currentNode,goalNode)
        dummy_node=node()
        if (verln==true)
            dummy_node.xy.y=currentNode.xy.y+cb[1].length*(abs(currentNode.xy.y-goalNode.xy.y)/(currentNode.xy.y-goalNode.xy.y))
            dummy_node.xy.x=dummy_node.xy.y*vect.m_findx+vect.b_findx
        else
            dummy_node.xy.x=currentNode.xy.x+cb[1].length*(abs(currentNode.xy.x-goalNode.xy.x)/(currentNode.xy.x-goalNode.xy.x))
            dummy_node.xy.y=dummy_node.xy.x*vect.m_findy+vect.b_findy
        end
        dummy_node.num=deepcopy(length(ocn.discretedom.nodes)+1)
        #Cable
        push!(cb[1].pth,deepcopy(currentNode))
        push!(cb[1].pth,deepcopy(dummy_node))
        push!(cb[2].pth,deepcopy(dummy_node))
        push!(cb[2].pth,deepcopy(goalNode))
        push!(oss_system.owp_MVcbls,deepcopy(cb[1]))
        push!(oss_system.owp_HVcbls,deepcopy(cb[2]))

        #OSS
        ossmv=bus()
        ossmv.node=deepcopy(dummy_node)
        ossmv.mva=crc.base_owp.mva
        ossmv.wnd=crc.base_owp.wnd
        #ossmv.node=currentNode.parent
        ossmv.base_cost=ocn.finance.FC_bld
        push!(ossmv.xfmrs,deepcopy(xs[1]))
        push!(oss_system.osss_owp,ossmv)
        push!(ocn.discretedom.nodes,deepcopy(dummy_node))
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

function opt_buildSystem_0(oss_systems_all,rnk,ocn,owpp_dec_mhv_orig,mvrng)
    owpp_dec_mhv_xpnd,oss_systems_xpnd=opt_setBreakDown(oss_systems_all,owpp_dec_mhv_orig)
#i_set=3
    for i_set=1:1:length(owpp_dec_mhv_xpnd)
        oss_systems=oss_systems_xpnd[i_set]
        owpp_dec_mhv=owpp_dec_mhv_xpnd[i_set]

        #take base information
        oss_system=opt_copyBaseEquipment(ocn.circuits[rnk],owpp_dec_mhv[1],ocn)
        pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_copyWindProfiles(ocn.owpps,oss_system)
        println(owpp_dec_mhv)
        #ppf_equipment_OSS_MOG_test(ocn,oss_system,owpp_dec_mhv)

        #=base_sys=0
        for crc in oss_systems[2:length(oss_systems)]
            base_sys=base_sys+crc.cost-crc.pcc_cbls[1].costs.ttl+(crc.pcc_cbls[1].costs.ttl/crc.pcc_cbls[1].length)*lof_pnt2pnt_dist(crc.pcc_cbls[1].pth[1].xy,oss_system.osss_mog[1].node.xy)
        end
        opt_ttlMvCirc(oss_system)
        #println(string(mox+oss_system.cost+base_sys)*" - "*string(ocn.circuits[rnk].cost))
        if (oss_system.cost+base_sys <= (ocn.circuits[rnk].cost)*1)=#
            #println("inside!!!")
            #store wind for transformer calculation at MOG
            #combinesSyss=Float64[]
            #push!(combinesSyss,oss_systems[1].decimal)
        ############################# update this section with new set calculation ########################
            if (length(oss_systems)>1)
                connect_Type=""
                #crc=oss_systems[2:length(oss_systems)][1]
                for crc in oss_systems[2:length(oss_systems)]
                    #push!(combinesSyss,crc.decimal)
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
                #new_coords=opt_reAdjust_oss(oss_system,mvrng,ocn.sys.mvCl)
                #oss_system=opt_reAdjust_cbls(oss_system,new_coords,ocn)
                #println(string(oss_system.cost)*" - "*string(ocn.circuits[rnk].cost)*" - "*string(combinesSyss))
                mvrng=6.6
                new_coords=opt_reAdjust_oss(oss_system,mvrng,ocn.sys.mvCl)
                oss_system=opt_reAdjust_cbls(oss_system,new_coords,ocn)
                if (oss_system.cost < ocn.circuits[rnk].cost)
                    #println(string(rnk)*" - "*connect_Type)
                    ocn.circuits[rnk]=deepcopy(oss_system)
                end
            end
        #else
            #println("skip!!!!")
        #end
    end
end

function opt_buildSystem_0_old(oss_systems_all,rnk,ocn,owpp_dec_mhv_orig,mvrng)
    owpp_dec_mhv_xpnd,oss_systems_xpnd=opt_setBreakDown(oss_systems_all,owpp_dec_mhv_orig)
#i_set=3
    for i_set=1:1:length(owpp_dec_mhv_xpnd)
        oss_systems=oss_systems_xpnd[i_set]
        owpp_dec_mhv=owpp_dec_mhv_xpnd[i_set]

        #take base information
        oss_system=opt_copyBaseEquipment(ocn.circuits[rnk],owpp_dec_mhv[1],ocn)
        pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_copyWindProfiles(ocn.owpps,oss_system)
        println(owpp_dec_mhv)
        #ppf_equipment_OSS_MOG_test(ocn,oss_system,owpp_dec_mhv)

        #=base_sys=0
        for crc in oss_systems[2:length(oss_systems)]
            base_sys=base_sys+crc.cost-crc.pcc_cbls[1].costs.ttl+(crc.pcc_cbls[1].costs.ttl/crc.pcc_cbls[1].length)*lof_pnt2pnt_dist(crc.pcc_cbls[1].pth[1].xy,oss_system.osss_mog[1].node.xy)
        end
        opt_ttlMvCirc(oss_system)
        #println(string(mox+oss_system.cost+base_sys)*" - "*string(ocn.circuits[rnk].cost))
        if (oss_system.cost+base_sys <= (ocn.circuits[rnk].cost)*1)=#
            #println("inside!!!")
            #store wind for transformer calculation at MOG
            #combinesSyss=Float64[]
            #push!(combinesSyss,oss_systems[1].decimal)
        ############################# update this section with new set calculation ########################
            if (length(oss_systems)>1)
                connect_Type=""
                #crc=oss_systems[2:length(oss_systems)][1]
                for crc in oss_systems[2:length(oss_systems)]
                    #push!(combinesSyss,crc.decimal)
                    #copy base cables and transformers
                    if (length(crc.owpps)==1)
                        println("inside crc.owpps")
                        oss_system,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_str8Oss2Oss_0(crc,oss_system,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
                    else

                        if (opt_fromScrtch(oss_system.osss_mog[1],crc,ocn.yaxisMajor))#check if full recalc required
                            connect_Type=connect_Type*"_Scratch"
                            oss2ossCbl,pth=opt_compoundCbl_0(crc,oss_system,ocn)
                            #oss2ossCbl,pth=opt_combineScratch(crc,oss_system,ocn)
                        else
                            connect_Type=connect_Type*"_Orig"
                            oss2ossCbl,pth=opt_compoundCbl_0(crc,oss_system,ocn)
                        end
                        #println(string(crc.pcc_cbls[1].num)*" - "*string(crc.pcc_cbls[1].size)*" || "*string(oss2ossCbl.num)*" - "*string(oss2ossCbl.size))
                        oss2ossCbl,oss_system=opt_combineCblRoute_0(pth,oss2ossCbl,oss_system)
                        pHv,wHv,p132,w132,p220,w220,p400,w400=opt_combineCblPW(oss_system,crc,pHv,wHv,p132,w132,p220,w220,p400,w400)
                    end
                end

                oss_system.osss_mog[1]=opt_mogXfmrs(oss_system.osss_mog[1],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)

                new_coords=opt_reAdjust_oss(oss_system,mvrng,ocn.sys.mvCl)
                oss_system=opt_reAdjust_cbls(oss_system,new_coords,ocn)
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
end


function opt_breakdownSystem_0(owpp_dec_mhv,dinner_array,dec,ocn,mvrng)
    for (i,dc) in enumerate(dinner_array)
        temp_dec_mhv=deepcopy(owpp_dec_mhv)
        push!(temp_dec_mhv,deepcopy(dc))
        if (sum(temp_dec_mhv) == dec)
            opt_bS_0(temp_dec_mhv,ocn,mvrng)
        #elseif (sum(owpp_dec_mhv) < dec)# original
        elseif (sum(temp_dec_mhv) < dec)
            opt_breakdownSystem_0(temp_dec_mhv,dinner_array[1:i-1],dec,ocn,mvrng)
        else
            #println("Error: decimal sum is over maximum!")
        end
    end
end

function opt_bS_0(owpp_dec_mhv,ocn,mvrng)
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
        opt_buildSystem_0(oss_systems,rnk,ocn,owpp_dec_mhv,mvrng)
    else
    end
end




function opt_compoundOSS_zeroSize(ocn)
    mvrng=cstF_mVrng(5,ocn.owpps[1].mva,ocn.owpps[1].wnd,ocn.finance,ocn.sys)
    owpps_dec_mhv=Array{Int32, 1}()
    owpps_tbl_mhv=Array{Array{Int8, 1}, 1}()
    owpp_tbl_mhv=Array{Int8, 1}()
    owpp_dec_mhv=Array{Int32, 1}()
    for crc in ocn.circuits
        push!(owpps_dec_mhv,deepcopy(crc.decimal))
        push!(owpps_tbl_mhv,deepcopy(crc.binary))
    end
#=
dec=23
println("Testing mode!!!!!!")
=#
    #for dec in owpps_dec_mhv[22:23]
    for dec in owpps_dec_mhv[7:length(owpps_dec_mhv)]

        owpp_tbl_mhv=owpps_tbl_mhv[dec]
        println("dec: "*string(owpp_tbl_mhv))
        min_Bin=findlast(x->x==1,owpp_tbl_mhv)
        osss_Pos=findfirst(x->x==1,owpp_tbl_mhv)
        zeros_mhv=findall(x->x==0,owpp_tbl_mhv[1:length(owpp_tbl_mhv)])
        min_dec=floor(Int32,2^(min_Bin-1))
        oss_dec=floor(Int32,2^(osss_Pos-1))
        for bin in reverse(owpps_tbl_mhv[min_dec+1:dec-1])

            #bin=reverse(owpps_tbl_mhv[min_dec+1:dec-1])[7]
            #println("bin: "*string(bin))
            owpp_dec_mhv=Int32[]
            if (findall(x->x==1,bin[1:osss_Pos]) == [])
                if (sum(bin[zeros_mhv]) < 1)
                    push!(owpp_dec_mhv,deepcopy(oss_dec))
                    dc=floor(Int32,top_bin2dec(bin))
                    push!(owpp_dec_mhv,deepcopy(dc))
                    if (sum(owpp_dec_mhv) == dec)
                        opt_bS_0(owpp_dec_mhv,ocn,mvrng)
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
                            opt_breakdownSystem_0(owpp_dec_mhv,dinner_array,dec,ocn,mvrng)
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

#circ=ocn.circuits[1]
#co_ords=new_coords
function opt_reAdjust_cbls(circ,co_ords,ocn)
    for mg_i=1:1:length(circ.osss_mog)
        circ.osss_mog[mg_i].node.xy.x=deepcopy(co_ords[1][mg_i])
        circ.osss_mog[mg_i].node.xy.y=deepcopy(co_ords[2][mg_i])
        circ.osss_mog[mg_i].node.goal=deepcopy(circ.osss_mog[mg_i].node.num)
        circ.osss_mog[mg_i].node.num=deepcopy(length(ocn.discretedom.nodes)+1)
        push!(ocn.discretedom.nodes,circ.osss_mog[mg_i].node)
    end
    for oss_i=1:1:length(circ.osss_owp)
        circ.osss_owp[oss_i].node.xy.x=deepcopy(co_ords[3][oss_i])
        circ.osss_owp[oss_i].node.xy.y=deepcopy(co_ords[4][oss_i])
        circ.osss_owp[oss_i].node.goal=deepcopy(circ.osss_owp[oss_i].node.num)
        circ.osss_owp[oss_i].node.num=deepcopy(length(ocn.discretedom.nodes)+1)
        push!(ocn.discretedom.nodes,circ.osss_owp[oss_i].node)
    end
    #set
    circ.owp_MVcbls=opt_updateMVC(circ,ocn)
    circ.owp_HVcbls=opt_updateHVC(circ,ocn)
    circ.oss2oss_cbls=opt_updateo2o(circ,ocn)
    circ.pcc_cbls=opt_updatePcc(circ,ocn)
    opt_ttlMvCirc(circ)
    return circ
end

function opt_updatePcc(circ,ocn)
    S=0
    wp=wind()
    for pcc_i=1:1:length(circ.pcc_cbls)
        pcc_mogNodeNum=deepcopy(circ.pcc_cbls[pcc_i].pth[1].num)
        circ.pcc_cbls[pcc_i].pth=node[]
        tale=node()
        hed=ocn.pccs[length(ocn.pccs)].node
        kble=cbl()
        for mog in circ.osss_mog
            if (mog.node.goal==pcc_mogNodeNum)
                tale=deepcopy(mog.node)
                S=deepcopy(mog.mva)
                wp=deepcopy(mog.wnd)
            else
            end
        end
        for oss in circ.osss_owp
            if (oss.node.goal==pcc_mogNodeNum)
                tale=deepcopy(oss.node)
                S=deepcopy(oss.mva)
                wp=deepcopy(oss.wnd)
            else
            end
        end
        pcc_l=lof_pnt2pnt_dist(tale.xy,hed.xy)
        if (pcc_l<circ.pcc_cbls[pcc_i].length)
            kble=cstF_nextSizeDownPcc(pcc_l,S,circ.pcc_cbls[pcc_i].elec.volt,wp,ocn.finance,circ.pcc_cbls[pcc_i].size,circ.pcc_cbls[pcc_i].num)
            circ.pcc_cbls[pcc_i]=kble
        elseif (pcc_l>circ.pcc_cbls[pcc_i].length)
            kble=cstF_nextSizeUpPcc(pcc_l,S,circ.pcc_cbls[pcc_i].elec.volt,wp,ocn.finance,circ.pcc_cbls[pcc_i].size,circ.pcc_cbls[pcc_i].num)
            circ.pcc_cbls[pcc_i]=kble
        else
        end
        push!(circ.pcc_cbls[pcc_i].pth,tale)
        push!(circ.pcc_cbls[pcc_i].pth,hed)
        new_edge=edge()
        new_edge.head=hed.num
        new_edge.tail=tale.num
        new_edge.lngth=pcc_l
        push!(ocn.discretedom.edges,new_edge)
    end
    return circ.pcc_cbls
end

function opt_updateo2o(circ,ocn)
    S=0
    wp=wind()
    for o2o_i=1:1:length(circ.oss2oss_cbls)
        o2o_ossNodeNum=deepcopy(circ.oss2oss_cbls[o2o_i].pth[1].num)
        o2o_mogNodeNum=deepcopy(circ.oss2oss_cbls[o2o_i].pth[length(circ.oss2oss_cbls[o2o_i].pth)].num)
        circ.oss2oss_cbls[o2o_i].pth=node[]
        tale=node()
        hed=node()
        kble=cbl()
        for mog in circ.osss_mog
            if (mog.node.goal==o2o_ossNodeNum)
                tale=deepcopy(mog.node)
                S=deepcopy(mog.mva)
                wp=deepcopy(mog.wnd)
            elseif (mog.node.goal==o2o_mogNodeNum)
                hed=deepcopy(mog.node)
            else
            end
        end
        o2o_l=lof_pnt2pnt_dist(tale.xy,hed.xy)
        if (o2o_l<circ.oss2oss_cbls[o2o_i].length)
            kble=cstF_nextSizeDown(o2o_l,S,circ.oss2oss_cbls[o2o_i].elec.volt,wp,ocn.finance,circ.oss2oss_cbls[o2o_i].size,circ.oss2oss_cbls[o2o_i].num)
            circ.oss2oss_cbls[o2o_i]=kble
        elseif (o2o_l>circ.oss2oss_cbls[o2o_i].length)
            kble=cstF_nextSizeUp(o2o_l,S,circ.oss2oss_cbls[o2o_i].elec.volt,wp,ocn.finance,circ.oss2oss_cbls[o2o_i].size,circ.oss2oss_cbls[o2o_i].num)
            circ.oss2oss_cbls[o2o_i]=kble
        else
        end
        push!(circ.oss2oss_cbls[o2o_i].pth,tale)
        push!(circ.oss2oss_cbls[o2o_i].pth,hed)
        new_edge=edge()
        new_edge.head=hed.num
        new_edge.tail=tale.num
        new_edge.lngth=o2o_l
        push!(ocn.discretedom.edges,new_edge)
    end
    return circ.oss2oss_cbls
end


function opt_updateHVC(circ,ocn)
    S=0
    wp=wind()
    for hvc_i=1:1:length(circ.owp_HVcbls)
        hvc_ossNodeNum=deepcopy(circ.owp_HVcbls[hvc_i].pth[1].num)
        hvc_mogNodeNum=deepcopy(circ.owp_HVcbls[hvc_i].pth[length(circ.owp_HVcbls[hvc_i].pth)].num)
        circ.owp_HVcbls[hvc_i].pth=node[]
        tale=node()
        hed=node()
        kble=cbl()
        #println(hvc_ossNodeNum)
        for oss in circ.osss_owp
            #println(oss.node.goal)
            if (oss.node.goal==hvc_ossNodeNum)
                tale=oss.node
                S=oss.mva
                #println("inside")
                #println(S)
                wp=oss.wnd
                break
            end
        end
        for mog in circ.osss_mog
            if (mog.node.goal==hvc_mogNodeNum)
                hed=mog.node
                break
            end
        end
        hv_l=lof_pnt2pnt_dist(tale.xy,hed.xy)
        if (hv_l<circ.owp_HVcbls[hvc_i].length)
            kble=cstF_nextSizeDown(hv_l,S,circ.owp_HVcbls[hvc_i].elec.volt,wp,ocn.finance,circ.owp_HVcbls[hvc_i].size,circ.owp_HVcbls[hvc_i].num)
            circ.owp_HVcbls[hvc_i]=kble
        elseif (hv_l>circ.owp_HVcbls[hvc_i].length)
            kble=cstF_nextSizeUp(hv_l,S,circ.owp_HVcbls[hvc_i].elec.volt,wp,ocn.finance,circ.owp_HVcbls[hvc_i].size,circ.owp_HVcbls[hvc_i].num)
            circ.owp_HVcbls[hvc_i]=kble
        else
        end
        push!(circ.owp_HVcbls[hvc_i].pth,tale)
        push!(circ.owp_HVcbls[hvc_i].pth,hed)
        new_edge=edge()
        new_edge.head=hed.num
        new_edge.tail=tale.num
        new_edge.lngth=hv_l
        push!(ocn.discretedom.edges,new_edge)
    end
    return circ.owp_HVcbls
end

#mvc_i=1
function opt_updateMVC(circ,ocn)
    for mvc_i=1:1:length(circ.owp_MVcbls)
        mvc_owpNodeNum=deepcopy(circ.owp_MVcbls[mvc_i].pth[1].num)
        mvc_ossNodeNum=deepcopy(circ.owp_MVcbls[mvc_i].pth[length(circ.owp_MVcbls[mvc_i].pth)].num)
        circ.owp_MVcbls[mvc_i].pth=node[]
        tale=bus()
        hed=bus()
        kble=cbl()
        for opp in circ.owpps
            if (opp.node.num==mvc_owpNodeNum)
                #opp.node.parent=node()
                tale=deepcopy(opp)
                break
            end
        end
        #oss=circ.osss_owp[1]
        for oss in circ.osss_owp
            if (oss.node.goal==mvc_ossNodeNum)
                #oss.node.parent=node()
                hed=deepcopy(oss)
                break
            end
        end
        for mog in circ.osss_mog
            if (mog.node.goal==mvc_ossNodeNum)
                #mog.node.parent=node()
                hed=deepcopy(mog)
                break
            end
        end
        #hed=circ.osss_mog[1].node
        #println(tail)
        #println(head)
        mv_l=lof_pnt2pnt_dist(tale.node.xy,hed.node.xy)
        if (mv_l<circ.owp_MVcbls[mvc_i].length)
            if (mv_l<ocn.sys.mvCl)
                mv_l=ocn.sys.mvCl
            end
            kble=cstF_MvCbl_nextSizeDown(mv_l,tale.mva,circ.owp_MVcbls[mvc_i].elec.volt,tale.wnd,ocn.finance,circ.owp_MVcbls[mvc_i].size,circ.owp_MVcbls[mvc_i].num)
            circ.owp_MVcbls[mvc_i]=kble
        elseif (mv_l>circ.owp_MVcbls[mvc_i].length)
            kble=cstF_MvCbl_nextSizeUp(mv_l,tale.mva,circ.owp_MVcbls[mvc_i].elec.volt,tale.wnd,ocn.finance,circ.owp_MVcbls[mvc_i].size,circ.owp_MVcbls[mvc_i].num)
            circ.owp_MVcbls[mvc_i]=kble
        else
        end
        push!(circ.owp_MVcbls[mvc_i].pth,tale.node)
        push!(circ.owp_MVcbls[mvc_i].pth,hed.node)
        new_edge=edge()
        new_edge.head=hed.node.num
        new_edge.tail=tale.node.num
        new_edge.lngth=mv_l
        push!(ocn.discretedom.edges,new_edge)
    end
    return circ.owp_MVcbls
end


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
#=
bfsmv=best_full_syss_cmv
bfshv=best_full_syss_chv
bfsmhv=best_full_syss_mvhv
cmv=bfsmv[1]
chv=bfshv[1]
=#
function combineAndrank(bfsmv,bfshv,bfsmhv)
    bfsmhv=comboCompare(bfsmhv,bfshv)
    bfsmhv=comboCompare(bfsmhv,bfsmv)
    bfsmhv=deepcopy(opt_sortFllSys(bfsmhv))
    return bfsmhv
end


function comboCompare(scs,bfs)
    eps=0.00001
    for cv in bfs
        crcCopy=false
        for sc in scs
            if (cv.cost<sc.cost+eps && cv.cost>sc.cost-eps)
                if (length(cv.osss_mog)==length(sc.osss_mog) && length(cv.osss_owp)==length(sc.osss_owp) && length(cv.owp_MVcbls)==length(sc.owp_MVcbls) && length(cv.owp_HVcbls)==length(sc.owp_HVcbls) && length(cv.oss2oss_cbls)==length(sc.oss2oss_cbls) && length(cv.pcc_cbls)==length(sc.pcc_cbls))
                    crcCopy=true
                end
            end
        end
        if (crcCopy==false)
            push!(scs,deepcopy(cv))
        end
    end
    return scs
end
