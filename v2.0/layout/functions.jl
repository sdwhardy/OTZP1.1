using Geodesy, Combinatorics, LinearAlgebra, Ipopt, JuMP, Cbc
include("struct.jl")
include("input_data/functions.jl")
include("circuits/functions.jl")

#orders OWPPS by distance to PCC - note mv_zone is only a convenient place holder and is nothing to do with the mv_zone size
function order_owpps(owpps,pcc)
    ordered_owpps=bus[]
    for wp in owpps
        wp.mv_zone=euclidian_distance(wp.node.xy,pcc.node.xy)
        ordered_owpps=(splice!(ordered_owpps, searchsorted(ordered_owpps,wp,by = v -> v.mv_zone), [wp]);ordered_owpps)
    end
    return ordered_owpps
end

#gives the nodes numbers for the owpps and pcc
function number_buses(ocean)
    for wp in ocean.owpps
        ocean.num=ocean.num+1
        wp.node.num=ocean.num
    end
    ocean.num=ocean.num+1
    ocean.pcc.node.num=ocean.num
    return ocean
end


#returns the euclidian distance between 2 points
function euclidian_distance(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
    return hyp
end

#sets PCC to (0,0) and OWPPs to cartesian coordinates of scale 1 tic=1km
#north_south=true for north
function utm_gps2xy(pcc,owpps,zone_utm,north_south)
    offset=xy()
    utm_desired = UTMfromLLA(zone_utm, north_south, wgs84)#sets UTM zone
    pcc_utm = utm_desired(LLA(pcc.node.gps.lat,pcc.node.gps.lng))#coverts to cartesian
    pcc.node.xy.x,pcc.node.xy.y=0,0#zeroes PCC
    offset.x=pcc_utm.x#saves offset for reverrse convertion
    offset.y=pcc_utm.y
    for wf in owpps
        _utm = utm_desired(LLA(wf.node.gps.lat,wf.node.gps.lng))
        wf.node.xy.x=(_utm.x-offset.x)/1000#converts to km
        wf.node.xy.y=(_utm.y-offset.y)/1000
    end
    return pcc,owpps,offset
end

#makes both Tbmv and Tbhv
function make_set_Tb(owpps,pcc,B,ocn)
    A=range(1,stop=length(owpps))
    A_squared=[2^(x-1) for x in A]
    Tb_basis_circuits=make_basis_of_Tb(owpps,pcc,B,ocn)
    for circ in Tb_basis_circuits
        if (circ[1].decimal in A_squared)
            push!(ocn.hv_circuits,circ)
            push!(ocn.mv_circuits,circ)
        else
            println("topology "*string(circ[1].decimal))
            push!(ocn.hv_circuits,[make_set_TbHV(ocn,deepcopy(circ[1]))])
            ocn.hv_circuits[length(ocn.hv_circuits)][1].id=string(ocn.hv_circuits[length(ocn.hv_circuits)][1].decimal)*"#h"
            push!(ocn.mv_circuits,[make_set_TbMV(ocn,deepcopy(circ[1]))])
            ocn.mv_circuits[length(ocn.mv_circuits)][1].id=string(ocn.mv_circuits[length(ocn.mv_circuits)][1].decimal)*"#m"
        end
    end
    return ocn
end

#Builds a tbmv topology
function make_set_TbMV(ocn,circ)
    mv_regionOWPP=owpps_within_mv_range(circ.owpps)
    #Get the coordinates of all OWPPS to be connected
    xys=Array{xy,1}()
    for owp in mv_regionOWPP
        push!(xys,deepcopy(owp.node.xy))
    end
    #find the location that minimizes connection distance
    oss_locationMV = oss_location_MV(xys,circ.base.mv_zone)
    mv_connections=mog2pcc_possibilities_noMOG(oss_locationMV,circ,ocn)
    #circ=owpps2oss(mv_connections,oss_locationMV,circ,ocn,false)
    circ=owpps2oss(mv_connections,oss_locationMV,circ,ocn,false)
    return circ
end



function oss_location_MV(xys,mv_max)
    m = Model(optimizer_with_attributes(Ipopt.Optimizer,"print_level"=>0))
    @variable(m, x)
    @variable(m, y)
    epsilon=1e-2


    @NLobjective(m, Min,sum(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2+epsilon) for i in 1:length(xys)))
    @NLconstraint(m, sqrt((xys[1].x-x)^2+(xys[1].y-y)^2) <= (mv_max-epsilon))

    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end

#finds which owpps can be connected via MV cables
function owpps_within_mv_range(owpps)
    mv_regionOWPP=bus[]
    push!(mv_regionOWPP,owpps[1])
    for owp in owpps[2:length(owpps)]
        mv_rng=true
        for mv_owp in mv_regionOWPP
            o2oL=euclidian_distance(mv_owp.node.xy,owp.node.xy)
            if (o2oL>=mv_owp.mv_zone+owp.mv_zone)
                mv_rng=false
            end
        end
        if (mv_rng==true)
            push!(mv_regionOWPP,owp)
        end
    end
    return mv_regionOWPP
end
#ocn.database["bits"]["hvdc"]=true
#Builds a tbhv topology
function make_set_TbHV(ocn,circ)
    #Get the coordinates of all OWPPS to be connected
    xys=Array{xy,1}()
    for owp in circ.owpps
        push!(xys,deepcopy(owp.node.xy))
    end
    #find the location that minimizes connection distance
    oss_locationHV = oss_location_HV(xys,circ.base,circ.pcc)
    hv_connections=mog2pcc_possibilities_noMOG(oss_locationHV,circ,ocn)
    #circ=owpps2oss(hv_connections,oss_locationHV,circ,ocn,true)
    circ=owpps2oss(hv_connections,oss_locationHV,circ,ocn,true)
    return circ
end


function oss_location_HV(xys,owp,pcc)
    m = Model(optimizer_with_attributes(Ipopt.Optimizer,"print_level"=>0))
    @variable(m, x)
    @variable(m, y)
    x_mn=findmin([owp.node.xy.x,pcc.node.xy.x])[1]    #x_mn=18.215405
    x_mx=findmax([owp.node.xy.x,pcc.node.xy.x])[1]
    y_mn=findmin([owp.node.xy.y,pcc.node.xy.y])[1]    #x_mn=18.215405
    y_mx=findmax([owp.node.xy.y,pcc.node.xy.y])[1]
    #@variable(m, lamda)
    epsilon=1e-2
    push!(xys,pcc.node.xy)
    @NLobjective(m, Min,sum(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2+epsilon) for i in 1:length(xys)))

    @constraint(m, y <= y_mx)
    @constraint(m, y >= y_mn)
    @constraint(m, x <= x_mx)
    @constraint(m, x >= x_mn)

    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end


#Constructs the common aspects of Set Tbmv and Tbhv
function make_basis_of_Tb(owpps,pcc,B,ocn)
    hmv_circuits=Vector{Array{circuit,1}}()
    A=range(1,stop=length(owpps))
    for b in B
        hmv_circuit=circuit()
        hmv_circuit.binary=b[1]
        hmv_circuit.decimal=b[2]
        first_owp=findfirst(x->x==1,b[1])[2]
        hmv_circuit.base=owpps[first_owp]
        attached_owpps=findall(x->x==1,b[1])
        net_mva=0.0
        net_wnd=Array{wind,1}()
        for owp in attached_owpps
            push!(net_wnd,owpps[owp[2]].wnd)
            net_mva=net_mva+owpps[owp[2]].mva
            push!(hmv_circuit.owpps,deepcopy(owpps[owp[2]]))
        end
        hmv_circuit.mva=net_mva
        hmv_circuit.wnd=find_netWind(net_wnd)
        hmv_circuit.pcc=deepcopy(pcc)
        if (length(attached_owpps)==1)
            hmv_circuit.id="mvhv"*string(hmv_circuit.decimal)
            hmv_circuit=owpp2pcc(hmv_circuit,ocn)
        end
        push!(hmv_circuits,[deepcopy(hmv_circuit)])
    end
    return hmv_circuits
end

#Builds set B and removes any unwanted connections
function make_set_B(owpps,pcc)
    A=range(1,stop=length(owpps))
    base_b=zeros(Int8,1,length(owpps))
    B=Array{Tuple{Array{Int8,2},Int64},1}()
    combos=collect(combinations(A))
    combos=filter_set_B(owpps,pcc,combos)
    for combo in combos
        b=deepcopy(base_b)
        dec=0
        for bits in combo
            b[bits]=1
            dec=dec+2^(bits-1)
        end
        push!(B,(b,deepcopy(dec)))
    end
    sort!(B, by = x -> x[2])
    return B
end

#removes any combinations where an upstream OWPP and the forward most OWPP are outside the inclusion angle
function filter_set_B(owpps,pcc,combos)
    bad_cons=bad_connections(owpps,pcc)
    #finds the set of indexes of all unwanted connections
    dont_wants=Array{Int32,1}()
    for (index,combo) in enumerate(combos)
        a=findmin(combo)
        for bad_con in bad_cons
            if (bad_con[2]==a[1])
                for bit in combo
                    if (bit==bad_con[1])
                        push!(dont_wants,deepcopy(index))
                        @goto next_connect
                    end
                end
            end
        end
        @label next_connect
    end
    #keeps all wanted connections
    to_keep=Vector{Array{Int64,1}}()
    for (index,combo) in enumerate(combos)
        keep=true
        for dont_want in dont_wants
            if (dont_want==index)
                keep=false
                break
            end
        end
        if (keep)
            push!(to_keep,deepcopy(combo))
        end
    end
    return to_keep
end

#calculates the angle between the 2 vectors from the upstream OWPP
function bad_connections(owpps,pcc)
    bad_cons=Array{Tuple,1}()
    for i=length(owpps):-1:2
        for j=i-1:-1:1
            theta=dot_product(owpps[i].node.xy,owpps[j].node.xy,pcc.node.xy)
            if (theta<1/sqrt(2))
                push!(bad_cons,(deepcopy(i),deepcopy(j)))
            end
        end
    end
    return bad_cons
end

#finds the dot product
function dot_product(owpp0,owpp1,pcc)
    x_ia=owpp0.x-owpp1.x
    y_ia=owpp0.y-owpp1.y
    x_ipcc=owpp0.x-pcc.x
    y_ipcc=owpp0.y-pcc.y
    v_ia=[x_ia;y_ia]./(sqrt(x_ia^2+y_ia^2))
    v_ipcc=[x_ipcc;y_ipcc]./(sqrt(x_ipcc^2+y_ipcc^2))
    theta=dot(v_ia,v_ipcc)
    return theta
end

#cs1=ocn.mv_circuits[255][1]
#cs1=finalize_cables(cs1,new_coords,ocn)
#=tst=deepcopy(ocn.hv_circuits[i][1])
new_coords=finalize_mog_location(tst,10e-2)
tst,new_size=finalize_cables(tst,new_coords,ocn,ks)

new_cable=deepcopy(mvac_cable(250.0,1,tst.owpps[1].wnd,ocn.database["cables"][string(66.0)][string(250.0)],ks))=#
####################### readjust circuit ###############
function finalize_circuits_layout(circs,ocn)
    ks=get_Cost_Data()
    for (i0,cs0) in enumerate(circs)
        for (i1,cs1) in enumerate(cs0)
            circs[i0][i1]=finalize_circuit_layout(cs1,ocn,ks)
        end
    end
    return circs
end

function finalize_circuit_layout(cs1,ocn,ks)
    orig=deepcopy(cs1)
    if ((length(cs1.mog)>0) && (length(cs1.MVcbls)>1))
        new_coords=finalize_mog_location(cs1,10e-2)
        cs1,new_size=finalize_cables(cs1,new_coords,ocn,ks)
        if (new_size==true)
            laps=0
            bsf=deepcopy(cs1)
            while (laps<5 && new_size==true)
                new_coords=finalize_mog_location(cs1,10e-2)
                cs1,new_size=finalize_cables(cs1,new_coords,ocn,ks)
                if (cs1.cost<bsf.cost)
                    println(string(cs1.decimal)*" new best: "*string(cs1.cost))
                    bsf=deepcopy(cs1)
                end
                laps=laps+1
            end
            cs1=deepcopy(bsf)
        end
    end
    if (orig.cost<cs1.cost)
        println(string(cs1.decimal)*" falling back on original layout.")
        cs1=deepcopy(orig)
    end
    return cs1
end

function finalize_cables(circ,co_ords,ocn,ks)
    new_size=false
    for mg_i=1:1:length(circ.mog)
        circ.mog[mg_i].node.xy.x=deepcopy(co_ords[1][mg_i])
        circ.mog[mg_i].node.xy.y=deepcopy(co_ords[2][mg_i])
    end
    for oss_i=1:1:length(circ.oss)
        circ.oss[oss_i].node.xy.x=deepcopy(co_ords[3][oss_i])
        circ.oss[oss_i].node.xy.y=deepcopy(co_ords[4][oss_i])
    end
    for mv_c in circ.MVcbls
        mv_c,new_size=update_MVcable(mv_c,ocn,ks,new_size)
    end

    for hv_c in circ.HVcbls
        if (hv_c.elec.volt==300.0)
            hv_c,new_size=update_HVDCcable(hv_c,ocn,ks,new_size)
        elseif (hv_c.mpc_ac==true)
            hv_c,new_size=update_MPCACcable(hv_c,ocn,ks,new_size)
        else
            hv_c,new_size=update_HVACcable(hv_c,ocn,ks,new_size)
        end
    end

    for hv_c in circ.PCCcbls
        if (hv_c.elec.volt==300.0)
            hv_c,new_size=update_HVDCcable(hv_c,ocn,ks,new_size)
        elseif (hv_c.mpc_ac==true)
            hv_c,new_size=update_MPCACcable(hv_c,ocn,ks,new_size)
        else
            hv_c,new_size=update_HVACcable(hv_c,ocn,ks,new_size)
        end
    end
    for hv_c in circ.O2Ocbls
        if (hv_c.elec.volt==300.0)
            hv_c,new_size=update_HVDCcable(hv_c,ocn,ks,new_size)
        elseif (hv_c.mpc_ac==true)
            hv_c,new_size=update_MPCACcable(hv_c,ocn,ks,new_size)
        else
            hv_c,new_size=update_HVACcable(hv_c,ocn,ks,new_size)
        end
    end
    circ=total_circuit_cost(circ)
    return circ,new_size
end

function update_MVcable(mv_c,ocn,ks,new_size)
    mv_c.length=euclidian_distance(mv_c.path[1].xy,mv_c.path[2].xy)
    new_cable=deepcopy(mvac_cable(mv_c.mva,mv_c.length,mv_c.wnd,ocn.database["cables"][string(mv_c.elec.volt)][string(mv_c.mva)],ks))
    if (new_cable.num != mv_c.num || new_cable.size != mv_c.size)
        new_size=true
    end
    new_cable.costs.grand_ttl=new_cable.costs.ttl
    mv_c.relia=deepcopy(new_cable.relia)
    mv_c.costs=deepcopy(new_cable.costs)
    mv_c.elec=deepcopy(new_cable.elec)
    mv_c.size=deepcopy(new_cable.size)
    mv_c.num=deepcopy(new_cable.num)
    return mv_c,new_size
end

function update_HVACcable(hv_c,ocn,ks,new_size)
    hv_c.length=euclidian_distance(hv_c.path[1].xy,hv_c.path[2].xy)
    new_cable=deepcopy(hvac_cable(hv_c.mva,hv_c.length,hv_c.wnd,ocn.database["cables"][string(hv_c.elec.volt)][string(hv_c.mva)],ks))
    if (new_cable.num != hv_c.num || new_cable.size != hv_c.size)
        new_size=true
    end
    hv_c.relia=deepcopy(new_cable.relia)
    new_cable.costs.grand_ttl=new_cable.costs.ttl
    hv_c.costs=deepcopy(new_cable.costs)
    hv_c.elec=deepcopy(new_cable.elec)
    hv_c.size=deepcopy(new_cable.size)
    hv_c.num=deepcopy(new_cable.num)
    hv_c.reactors=deepcopy(new_cable.reactors)
    return hv_c,new_size
end

function update_HVDCcable(hv_c,ocn,ks,new_size)
    hv_c.length=euclidian_distance(hv_c.path[1].xy,hv_c.path[2].xy)
    new_cable=deepcopy(hvdc_cable(hv_c.mva,hv_c.length,hv_c.wnd,ocn.database["cables"][string(hv_c.elec.volt)][string(hv_c.mva)],ks))
    if (new_cable.num != hv_c.num || new_cable.size != hv_c.size)
        new_size=true
    end
    hv_c.relia=deepcopy(new_cable.relia)
    new_cable.costs.grand_ttl=new_cable.costs.ttl
    hv_c.costs=deepcopy(new_cable.costs)
    hv_c.elec=deepcopy(new_cable.elec)
    hv_c.size=deepcopy(new_cable.size)
    hv_c.num=deepcopy(new_cable.num)
    return hv_c,new_size
end

function update_MPCACcable(hv_c,ocn,ks,new_size)
    hv_c.length=euclidian_distance(hv_c.path[1].xy,hv_c.path[2].xy)/2
    new_cable=deepcopy(hvac_cable(hv_c.mva,hv_c.length,hv_c.wnd,ocn.database["cables"][string(hv_c.elec.volt)][string(hv_c.mva)],ks))
    if (new_cable.num != hv_c.num || new_cable.size != hv_c.size)
        new_size=true
    end
    hv_c.relia=deepcopy(new_cable.relia)
    new_cable.costs.grand_ttl=new_cable.costs.ttl*2+hv_c.plat.costs.ttl
    hv_c.costs=deepcopy(new_cable.costs)
    hv_c.elec=deepcopy(new_cable.elec)
    hv_c.size=deepcopy(new_cable.size)
    hv_c.num=deepcopy(new_cable.num)
    hv_c.reactors=deepcopy(new_cable.reactors)
    return hv_c,new_size
end
#=NOTE ORIGINAL AND FUNCTIONAL b4 new size
function finalize_cables(circ,co_ords,ocn,ks)
    for mg_i=1:1:length(circ.mog)
        circ.mog[mg_i].node.xy.x=deepcopy(co_ords[1][mg_i])
        circ.mog[mg_i].node.xy.y=deepcopy(co_ords[2][mg_i])
    end
    for oss_i=1:1:length(circ.oss)
        circ.oss[oss_i].node.xy.x=deepcopy(co_ords[3][oss_i])
        circ.oss[oss_i].node.xy.y=deepcopy(co_ords[4][oss_i])
    end
    for mv_c in circ.MVcbls
        mv_c=update_MVcable(mv_c,ocn,ks)
    end

    for hv_c in circ.HVcbls
        if (hv_c.elec.volt==300.0)
            hv_c=update_HVDCcable(hv_c,ocn,ks)
        elseif (hv_c.mpc_ac==true)
            hv_c=update_MPCACcable(hv_c,ocn,ks)
        else
            hv_c=update_HVACcable(hv_c,ocn,ks)
        end
    end

    for hv_c in circ.PCCcbls
        if (hv_c.elec.volt==300.0)
            hv_c=update_HVDCcable(hv_c,ocn,ks)
        elseif (hv_c.mpc_ac==true)
            hv_c=update_MPCACcable(hv_c,ocn,ks)
        else
            hv_c=update_HVACcable(hv_c,ocn,ks)
        end
    end
    for hv_c in circ.O2Ocbls
        if (hv_c.elec.volt==300.0)
            hv_c=update_HVDCcable(hv_c,ocn,ks)
        elseif (hv_c.mpc_ac==true)
            hv_c=update_MPCACcable(hv_c,ocn,ks)
        else
            hv_c=update_HVACcable(hv_c,ocn,ks)
        end
    end
    circ=total_circuit_cost(circ)
    return circ
end

function update_MVcable(mv_c,ocn,ks)
    mv_c.length=euclidian_distance(mv_c.path[1].xy,mv_c.path[2].xy)
    new_cable=deepcopy(mvac_cable(mv_c.mva,mv_c.length,mv_c.wnd,ocn.database["cables"][string(mv_c.elec.volt)][string(mv_c.mva)],ks))
    mv_c.costs=deepcopy(new_cable.costs)
    mv_c.elec=deepcopy(new_cable.elec)
    mv_c.size=deepcopy(new_cable.size)
    mv_c.num=deepcopy(new_cable.num)
    return mv_c
end

function update_HVACcable(hv_c,ocn,ks)
    hv_c.length=euclidian_distance(hv_c.path[1].xy,hv_c.path[2].xy)
    new_cable=deepcopy(hvac_cable(hv_c.mva,hv_c.length,hv_c.wnd,ocn.database["cables"][string(hv_c.elec.volt)][string(hv_c.mva)],ks))
    new_cable.costs.grand_ttl=new_cable.costs.ttl
    hv_c.costs=deepcopy(new_cable.costs)
    hv_c.elec=deepcopy(new_cable.elec)
    hv_c.size=deepcopy(new_cable.size)
    hv_c.num=deepcopy(new_cable.num)
    hv_c.reactors=deepcopy(new_cable.reactors)
    return hv_c
end

function update_HVDCcable(hv_c,ocn,ks)
    hv_c.length=euclidian_distance(hv_c.path[1].xy,hv_c.path[2].xy)
    new_cable=deepcopy(hvdc_cable(hv_c.mva,hv_c.length,hv_c.wnd,ocn.database["cables"][string(hv_c.elec.volt)][string(hv_c.mva)],ks))
    new_cable.costs.grand_ttl=new_cable.costs.ttl
    hv_c.costs=deepcopy(new_cable.costs)
    hv_c.elec=deepcopy(new_cable.elec)
    hv_c.size=deepcopy(new_cable.size)
    hv_c.num=deepcopy(new_cable.num)
    return hv_c
end

function update_MPCACcable(hv_c,ocn,ks)
    hv_c.length=euclidian_distance(hv_c.path[1].xy,hv_c.path[2].xy)/2
    new_cable=deepcopy(hvac_cable(hv_c.mva,hv_c.length,hv_c.wnd,ocn.database["cables"][string(hv_c.elec.volt)][string(hv_c.mva)],ks))
    new_cable.costs.grand_ttl=new_cable.costs.ttl*2+hv_c.plat.costs.ttl
    hv_c.costs=deepcopy(new_cable.costs)
    hv_c.elec=deepcopy(new_cable.elec)
    hv_c.size=deepcopy(new_cable.size)
    hv_c.num=deepcopy(new_cable.num)
    hv_c.reactors=deepcopy(new_cable.reactors)
    return hv_c
end=#

#=ks=get_Cost_Data()
tst=deepcopy(ocn.hv_circuits[45][1])
new_coords=finalize_mog_location(ocn.hv_circuits[45][1],10e-2)
plot_circuit(tst)
for os in tst.mog
    println(os.node.xy)
end
new_coords=finalize_mog_location(tst,10e-2)
tst=finalize_cables(tst,new_coords,ocn,ks)=#
function finalize_mog_location(topology,epsilon)

    mog_xys=Array{Tuple{xy,Int64},1}()
    oss_xys=Array{Tuple{xy,Int64},1}()
    for mog in topology.mog
        push!(mog_xys,(mog.node.xy,mog.node.num))
    end
    for os in topology.oss
        push!(oss_xys,(os.node.xy,os.node.num))
    end
    m = Model(optimizer_with_attributes(Ipopt.Optimizer,"print_level"=>0))



    @variables(m, begin
               MOG_y[i=1:length(mog_xys)], (start = mog_xys[i][1].y, base_name = "MOGY_"*string(mog_xys[i][2]))
               MOG_x[i=1:length(mog_xys)], (start = mog_xys[i][1].x, base_name = "MOGX_"*string(mog_xys[i][2]))
           end)

       @variables(m, begin
                  OSS_y[i=1:length(oss_xys)], (start = oss_xys[i][1].y, base_name = "OSSY_"*string(oss_xys[i][2]))
                  OSS_x[i=1:length(oss_xys)], (start = oss_xys[i][1].x, base_name = "OSSX_"*string(oss_xys[i][2]))
              end)
#limits
      x_mn=topology.pcc.node.xy.x    #x_mn=18.215405
      x_mx=topology.pcc.node.xy.x#x_mx=22.73574
      y_mn=topology.pcc.node.xy.y
      y_mx=topology.pcc.node.xy.y
      for co_ord in topology.owpps
          if (co_ord.node.xy.x<x_mn)
              x_mn=deepcopy(co_ord.node.xy.x)
          elseif (co_ord.node.xy.x>x_mx)
              x_mx=deepcopy(co_ord.node.xy.x)
          end
          if (co_ord.node.xy.y<y_mn)
              y_mn=deepcopy(co_ord.node.xy.y)
          elseif (co_ord.node.xy.y>y_mx)
              y_mx=deepcopy(co_ord.node.xy.y)
          end
      end
      for i=1:1:length(mog_xys)
          @constraint(m, MOG_y[i] <= y_mx)
          @constraint(m, MOG_y[i] >= y_mn)
          @constraint(m, MOG_x[i] <= x_mx)
          @constraint(m, MOG_x[i] >= x_mn)
      end
      for i=1:1:length(oss_xys)
          @constraint(m, OSS_y[i] <= y_mx)
          @constraint(m, OSS_y[i] >= y_mn)
          @constraint(m, OSS_x[i] <= x_mx)
          @constraint(m, OSS_x[i] >= x_mn)
      end

    connections=Array{Tuple{Any,Any,Any,Any,Float64},1}()
    mv_rng=0
    for mv_c in topology.MVcbls
        for owp in topology.owpps
            if (owp.node.num==mv_c.path[1].num)
                mv_rng=owp.mv_zone
                break
            end
        end

        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==mv_c.path[length(mv_c.path)].num)
                push!(connections,(mv_c.path[1].xy.x, mv_c.path[1].xy.y, MOG_x[i], MOG_y[i], mv_c.costs.perkm_ttl))
                @NLconstraint(m,sqrt(abs(connections[length(connections)][1]-connections[length(connections)][3])^2+abs(connections[length(connections)][2]-connections[length(connections)][4])^2+epsilon) <= (mv_rng))
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==mv_c.path[length(mv_c.path)].num)
                push!(connections,(mv_c.path[1].xy.x, mv_c.path[1].xy.y, OSS_x[i], OSS_y[i], mv_c.costs.perkm_ttl))
                @NLconstraint(m,sqrt(abs(connections[length(connections)][1]-connections[length(connections)][3])^2+abs(connections[length(connections)][2]-connections[length(connections)][4])^2+epsilon) <= (mv_c.length+epsilon))#(mv_c.length+epsilon))
            end
        end
    end

    #PCC_connections=Array{Tuple{VariableRef,VariableRef,Float64,Float64,Float64},1}()
    for pcc_c in topology.PCCcbls
        if (pcc_c.mpc_ac==true)
            mp_ac=2
        else
            mp_ac=1
        end
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==pcc_c.path[1].num)
                push!(connections,(MOG_x[i], MOG_y[i], pcc_c.path[length(pcc_c.path)].xy.x, pcc_c.path[length(pcc_c.path)].xy.y, pcc_c.costs.perkm_ttl))
                @NLconstraint(m,sqrt(abs(connections[length(connections)][1]-connections[length(connections)][3])^2+abs(connections[length(connections)][2]-connections[length(connections)][4])^2+epsilon) <= (mp_ac*pcc_c.mx_rng))
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==pcc_c.path[length(pcc_c.path)].num)
                push!(connections,(OSS_x[i], OSS_y[i], pcc_c.path[length(pcc_c.path)].xy.x, pcc_c.path[length(pcc_c.path)].xy.y, pcc_c.costs.perkm_ttl))
                @NLconstraint(m,sqrt(abs(connections[length(connections)][1]-connections[length(connections)][3])^2+abs(connections[length(connections)][2]-connections[length(connections)][4])^2+epsilon) <= (mp_ac*pcc_c.mx_rng))
            end
        end
    end

    #HV_connections=Array{Tuple{VariableRef,VariableRef,VariableRef,VariableRef,Float64},1}()
    for hv_c in topology.HVcbls
        if (hv_c.mpc_ac==true)
            mp_ac=2
        else
            mp_ac=1
        end
        tail_x=VariableRef
        tail_y=VariableRef
        head_x=VariableRef
        head_y=VariableRef
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==hv_c.path[length(hv_c.path)].num)
                head_x=MOG_x[i]
                head_y=MOG_y[i]
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==hv_c.path[1].num)
                tail_x=OSS_x[i]
                tail_y=OSS_y[i]
            end
        end
        push!(connections,(tail_x, tail_y, head_x, head_y, hv_c.costs.perkm_ttl))
        @NLconstraint(m,sqrt(abs(connections[length(connections)][1]-connections[length(connections)][3])^2+abs(connections[length(connections)][2]-connections[length(connections)][4])^2+epsilon) <= (mp_ac*hv_c.mx_rng))
    end
    #o2o_connections=Array{Tuple{VariableRef,VariableRef,VariableRef,VariableRef,Float64},1}()
    for o2o_c in topology.O2Ocbls
        if (o2o_c.mpc_ac==true)
            mp_ac=2
        else
            mp_ac=1
        end
        tail_x=VariableRef
        tail_y=VariableRef
        head_x=VariableRef
        head_y=VariableRef
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==o2o_c.path[1].num)
                tail_x=MOG_x[i]
                tail_y=MOG_y[i]
            end
            if (mg[2]==o2o_c.path[length(o2o_c.path)].num)
                head_x=MOG_x[i]
                head_y=MOG_y[i]
            end
        end

        push!(connections,(tail_x, tail_y, head_x, head_y, o2o_c.costs.perkm_ttl))
        @NLconstraint(m,sqrt(abs(connections[length(connections)][1]-connections[length(connections)][3])^2+abs(connections[length(connections)][2]-connections[length(connections)][4])^2+epsilon) <= (mp_ac*o2o_c.mx_rng))
    end

    @NLobjective(m, Min,sum(sqrt(((abs(connections[i][1]-connections[i][3])^2+abs(connections[i][2]-connections[i][4])^2)+epsilon))*connections[i][5] for i in 1:length(connections)))
    optimize!(m)
    mxy=(JuMP.value.(MOG_x),JuMP.value.(MOG_y),JuMP.value.(OSS_x),JuMP.value.(OSS_y))
    return mxy
end
