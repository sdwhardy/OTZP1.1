
function lof_layoutEez_basis()
        #reading data from excel input files
        ocean=eez()
        ocean.pccs=lof_getPccData()
        ocean.owpps=lof_getOwppData()
        ocean.sys=lof_getSysData()
        ocean.finance=lof_getCfsData()
        ocean.eqp_data=lof_getEqpData()
        #transformation from gps to cartesian
        ocean.base=lof_bseCrd(ocean)#find base coordinates
        print("Base coordinates are: ")
        println(ocean.base)
        lof_gps2cartesian4bus(ocean.owpps,ocean.base)#projects owpps onto cartesian plane
        lof_gps2cartesian4bus(ocean.pccs,ocean.base)#projects pccs onto cartesian plane
        lof_transformAxis(ocean)
        return ocean
end
#**
function lof_layoutEez_expand(ocn,pcc)
    #set area of owpp
    lof_setAreaOwpp(ocn)
    #set range of MV
    lof_mVrng(ocn)

    ##########Printing
    for value in ocn.pccs
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    for value in ocn.owpps
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    println("GPS coordinates projected onto cartesian plane.")
    println("Axis transformed.")
    return ocn
end

function lof_setAreaOwpp(ocn)
    for op in ocn.owpps
        #radius
        op.zone=ocn.sys.mvCl
    end
end
#**
function lof_mVrng(ocn)
    for owpp in ocn.owpps
        mvrng=cstF_mVrng(6,owpp.mva,owpp.wnd,ocn.finance,ocn.sys,ocn.eqp_data)
        owpp.mv_zone=deepcopy(mvrng)
    end
end
#**
#**
############################ Order the OWPP with chosen PCC ####################
################################################################################
function lof_order2Pcc(ocn,pcc)
    lngth_owpps=Array{Tuple{Float64,Int64},1}()
    ordrd_owpps=Array{bus,1}()
    for (indx,owp) in enumerate(ocn.owpps)
        owp.num=deepcopy(indx)
        push!(lngth_owpps,(deepcopy(lof_pnt2pnt_dist(owp.node.xy,pcc.node.xy)),deepcopy(owp.num)))
    end
    sort!(lngth_owpps, by = x -> x[1])
    for indx in lngth_owpps
        for owp in ocn.owpps
            if (indx[2] == owp.num)
                push!(ordrd_owpps,owp)
            end
        end
    end
    for (i,owp) in enumerate(ordrd_owpps)
        owp.num=deepcopy(i)
        owp.node.num=deepcopy(i)
    end
    for (i,pc) in enumerate(ocn.pccs)
        pc.node.num=deepcopy(i+length(ocn.owpps))
    end
    ocn.buses=length(ocn.owpps)+length(ocn.pccs)
    return ordrd_owpps
end

############################ GPS to cartesian transform ########################
################################################################################
#sets the gps coords that used as reference coords
#**
function lof_bseCrd(ocean)
    base=gps()
    #for india (type layouts)
    if ocean.owpps[length(ocean.owpps)].node.gps.lat < ocean.pccs[length(ocean.pccs)].node.gps.lat
        base.lat=ocean.owpps[length(ocean.owpps)].node.gps.lat#base lat
        base.lng=ocean.owpps[length(ocean.owpps)].node.gps.lng#base long
    #for belgium (type layouts)
elseif ocean.pccs[length(ocean.pccs)].node.gps.lat < ocean.owpps[length(ocean.owpps)].node.gps.lat
        base.lat=ocean.pccs[length(ocean.pccs)].node.gps.lat#base lat
        base.lng=ocean.pccs[length(ocean.pccs)].node.gps.lng#base long
    else
        error("No proper base coordinates system established!")
    end
    return base
end

#as lattitude changes number of km should be updated
#**
function lof_gps2cartesian4bus(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.node.xy.x=lof_deg2lgth(value.node.gps.lng-base.lng,lof_lg1deg(value.node.gps.lat,lnthLT))
        value.node.xy.y=lof_deg2lgth(value.node.gps.lat-base.lat,lnthLT)
    end
end
#changes angle to an arc length
#**
function lof_deg2lgth(d,dPl)
    return d*dPl
end
#calculates length of 1 deg of longitude at given lattitude
#**
function lof_lg1deg(lat,lngth)
    return cos(lof_d2r(lat))*lngth
end

#Change degrees to radians
#**
function lof_d2r(deg)
    return deg*pi/180
end
#rotates and slides cartesian axis
function lof_transformAxis(ocn)
    offset=lof_rotateAxis(ocn)
    lof_slideAxis(ocn,offset)
end
#finds angle to rotate and applies to owpps and pccs
#rotates axis to align n-s with y
#**
function lof_rotateAxis(ocn)
    theta=atan((ocn.pccs[length(ocn.pccs)].node.xy.x-ocn.owpps[length(ocn.owpps)].node.xy.x)/(ocn.owpps[length(ocn.owpps)].node.xy.y-ocn.pccs[length(ocn.pccs)].node.xy.y))
    ocn.theta=theta
    offset=0.0
    offset=lof_rotateGroup4bus(ocn.owpps,theta,offset)
    offset=lof_rotateGroup4bus(ocn.pccs,theta,offset)
    ocn.offset=deepcopy(offset)
    return ocn.offset
end
#loops through to apply rotations for a specified group
#**
function lof_rotateGroup4bus(locations,theta,os)
    for value in locations
        xy=lof_rotatePnt(value.node.xy.x,value.node.xy.y,theta)
        value.node.xy.x=xy[1]
        value.node.xy.y=xy[2]
        if value.node.xy.x<os
            os=value.node.xy.x
        end
    end
    return os
end
#applies rotational matrix individual coordinates
#**
function lof_rotatePnt(x,y,theta)
    co_od=[x y]
    rotated=co_od*[cos(theta) -1*sin(theta);sin(theta) cos(theta)]
    return rotated
end
#translates the entire region by specified offset
#sets unique IDs for owpps and pccs
#**
function lof_slideAxis(ocn,os)
    for value in ocn.owpps
        value.node.xy.x=value.node.xy.x-os
    end
    for value in ocn.pccs
        value.node.xy.x=value.node.xy.x-os
    end
end
#returns the hypotenuse distance between 2 cartesian points
#minimum distance for a path is 1km
#**
function lof_pnt2pnt_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
    return hyp
end
############################ Reads Data from Excel ########################
################################################################################
#**#reads the data on PCCs
function lof_getPccData()
    df = DataFrame(XLSX.readtable("Zero_size//layout//data.xlsx", "pcc_data")...)
    pccs=Array{bus,1}()
    for index=1:length(df[:, 1])
        dummy_bus=bus()
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.kV=df.kv[index]
        push!(pccs,deepcopy(dummy_bus))
        pccs[length(pccs)].num=length(pccs)
    end
    return pccs
end
#**
#**#reads OWPP data in excel file
function lof_getOwppData()
    df = DataFrame(XLSX.readtable("Zero_size//layout//data.xlsx", "owpp_data")...)
    owpps=Array{bus,1}()
    for index=1:length(df[:, 1])
        dummy_bus=bus()
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.name=df.name[index]
        dummy_bus.mva=df.mva[index]
        push!(owpps,deepcopy(dummy_bus))
        owpps[length(owpps)].num=length(owpps)
    end
    lof_getOwppWindData(owpps)
    return owpps
end
#**Reads Wind data from excel file
function lof_getOwppWindData(owpps)
    df = DataFrame(XLSX.readtable("Zero_size//layout//data.xlsx", "wind_data")...)
    for owpp in owpps
        wnd=wndF_wndPrf([getproperty(df,Symbol(owpp.name))])
        owpp.wnd=deepcopy(wnd)
    end
end
#**
#**reads system data from Excel file
function lof_getSysData()
    sht = XLSX.readxlsx("Zero_size//layout//data.xlsx")["sys_data"]
    sys=system()
    sys.nogoNum=sht["B1"]
    sys.prec=sht["B2"]
    sys.mwPerKm=sht["B3"]
    sys.mvCl=sht["B4"]
    return sys
end
#Sets values of all cost factors discribed below
#the structure of object ks is described in file cst_structure.jl
function lof_getCfsData()
    sht = XLSX.readxlsx("Zero_size//layout//data.xlsx")["cost_data"]
    ks=cstS_ks()#create an instance of object ks
    ks.FC_ac=sht["B1"]#fixed AC cost
    ks.FC_dc=sht["B2"]#fixed DC cost
    ks.dc=sht["B3"]#penalization factor for different than 2 xfrms
    ks.f_ct=sht["B4"]#generating plant variable cost
    ks.p_ct=sht["B5"]#substructure variable cost
    ks.c_ct=sht["B6"]#hvdc converter variable cost
    ks.Qc_oss=sht["B7"]#M£/MVAr
    ks.Qc_pcc=sht["B8"]#M£/MVAr
    ks.life=sht["B9"]#lifetime of wind farm
    ks.T_op=sht["B10"]#Operational lifetime in hours
    ks.E_op=sht["B11"]#Energy price £/Wh
    ks.cf=sht["B12"]#Capitalization factor
    ks.FC_bld=sht["B13"]#Build cost
    ks.p2e=sht["B14"]#pounds/euro exchange
    return ks
end
#Gets the data for equipment to be used
function lof_getEqpData()
    df = DataFrame(XLSX.readtable("Zero_size//layout//data.xlsx", "eqp_data")...)
    eqps=eqp_data()#create an instance of object ks
    for index=1:length(df[:, 1])
        #get 33kv cable data
        if (typeof(df.kv_33[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_33[index])
            push!(cbd,df.mm2_33[index])
            push!(cbd,df.ohms_per_km_33[index])
            push!(cbd,df.nf_per_km_33[index])
            push!(cbd,df.Amps_33[index])
            push!(cbd,df.euro_per_m_33[index])
            push!(cbd,df.mh_per_km_33[index])
            push!(eqps.cbls33kV,deepcopy(cbd))
        end

        #get 66kv cable data
        if (typeof(df.kv_66[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_66[index])
            push!(cbd,df.mm2_66[index])
            push!(cbd,df.ohms_per_km_66[index])
            push!(cbd,df.nf_per_km_66[index])
            push!(cbd,df.Amps_66[index])
            push!(cbd,df.euro_per_m_66[index])
            push!(cbd,df.mh_per_km_66[index])
            push!(eqps.cbls66kV,deepcopy(cbd))
        end

        #get 132kv cable data
        if (typeof(df.kv_132[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_132[index])
            push!(cbd,df.mm2_132[index])
            push!(cbd,df.ohms_per_km_132[index])
            push!(cbd,df.nf_per_km_132[index])
            push!(cbd,df.Amps_132[index])
            push!(cbd,df.euro_per_m_132[index])
            push!(cbd,df.mh_per_km_132[index])
            push!(eqps.cbls132kV,deepcopy(cbd))
        end

        #get 220kv cable data
        if (typeof(df.kv_220[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_220[index])
            push!(cbd,df.mm2_220[index])
            push!(cbd,df.ohms_per_km_220[index])
            push!(cbd,df.nf_per_km_220[index])
            push!(cbd,df.Amps_220[index])
            push!(cbd,df.euro_per_m_220[index])
            push!(cbd,df.mh_per_km_220[index])
            push!(eqps.cbls220kV,deepcopy(cbd))
        end


        #get 400kv cable data
        if (typeof(df.kv_400[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_400[index])
            push!(cbd,df.mm2_400[index])
            push!(cbd,df.ohms_per_km_400[index])
            push!(cbd,df.nf_per_km_400[index])
            push!(cbd,df.Amps_400[index])
            push!(cbd,df.euro_per_m_400[index])
            push!(cbd,df.mh_per_km_400[index])
            push!(eqps.cbls400kV,deepcopy(cbd))
        end
    end
    return eqps
end


#**
function lof_lineDirection(nd_tail,nd_head)
    if (abs(nd_tail.xy.y-nd_head.xy.y) <  abs(nd_tail.xy.x-nd_head.xy.x))
        vertLn=false
    else
        vertLn=true
    end
    return vertLn
end

#**
function lof_getStr8line(nd_hd,nd_tl)
    dummy_line=line()
    dummy_line.xmx=max(nd_hd.xy.x,nd_tl.xy.x)
    dummy_line.xmn=min(nd_hd.xy.x,nd_tl.xy.x)
    dummy_line.ymx=max(nd_hd.xy.y,nd_tl.xy.y)
    dummy_line.ymn=min(nd_hd.xy.y,nd_tl.xy.y)

    if (nd_hd.xy.x==nd_tl.xy.x)
        dummy_line.b_findy,dummy_line.m_findy=reverse([[nd_hd.xy.x,nd_tl.xy.x+(10^-5)] ones(2)]\[nd_hd.xy.y,nd_tl.xy.y])
    else
        dummy_line.b_findy,dummy_line.m_findy=reverse([[nd_hd.xy.x,nd_tl.xy.x] ones(2)]\[nd_hd.xy.y,nd_tl.xy.y])
    end

    if (nd_hd.xy.y==nd_tl.xy.y)
        dummy_line.b_findx,dummy_line.m_findx=reverse([[nd_hd.xy.y,nd_tl.xy.y+(10^-5)] ones(2)]\[nd_hd.xy.x,nd_tl.xy.x])
    else
        dummy_line.b_findx,dummy_line.m_findx=reverse([[nd_hd.xy.y,nd_tl.xy.y] ones(2)]\[nd_hd.xy.x,nd_tl.xy.x])
    end
    return dummy_line
end
