function lof_layoutEez()
    ocean=eez()
    ocean.id_count=0
    ocean.pccs,ocean.id_count=lof_getPccData(ocean.id_count)
    ocean.owpps,ocean.id_count=lof_getOwppData(ocean.id_count)
    ocean.bndryPnts=lof_getBndData()
    ocean.sys=lof_getSysData()
    for i=1:ocean.sys.nogoNum
        nogoT=nogo()
        nogoT.nodes=lof_getNoGoData(i)
        push!(ocean.nogos,deepcopy(nogoT))
    end

    #transformation to cartesian
    ocean.base=lof_bseCrd(ocean)#find base coordinates
    print("Base coordinates are: ")
    println(ocean.base)
    lof_gps2cartesian(ocean.owpps,ocean.base)#projects owpps onto cartesian plane
    lof_gps2cartesian(ocean.pccs,ocean.base)#projects pccs onto cartesian plane
    lof_gps2cartesian(ocean.bndryPnts,ocean.base)#projects boundary points onto cartesian plane
    for i=1:ocean.sys.nogoNum
        lof_gps2cartesian(ocean.nogos[i].nodes,ocean.base)#projects boundary points onto cartesian plane
    end
    lof_transformAxis(ocean)

    #build domain "mesh"
    lof_nodify(ocean)

    ocean.buses=vcat(ocean.pccs, ocean.owpps)#collects all buses

    ##########Printing
    for value in ocean.pccs
        print(value.num)
        print(" - ")
        println(value.cntr.gps)
    end
    for value in ocean.owpps
        print(value.num)
        print(" - ")
        println(value.cntr.gps)
    end
    println("GPS coordinates projected onto cartesian plane.")
    println("Axis transformed.")
    return ocean
end
function lof_nodify(ocn)
    wbnd,ebnd=lof_bndPerimeter(ocn)

end

function lof_bndPerimeter(ocn)
    #finds mid line dividing east and west regions
    ys=Array{Float64,1}()
    for xys in ocn.bndryPnts
        push!(ys,xys.cntr.xy.y)
    end
    mxY=findmax(ys)
    mnY=findmin(ys)
    midLn=reverse([[ocn.bndryPnts[mxY[2]].cntr.xy.y,ocn.bndryPnts[mnY[2]].cntr.xy.y] ones(2)]\[ocn.bndryPnts[mxY[2]].cntr.xy.x,ocn.bndryPnts[mnY[2]].cntr.xy.x])

    #sorts to east and west bounding points
    wpntsT=Array{xy,1}()
    epntsT=Array{xy,1}()
    wys=Array{Float64,1}()
    eys=Array{Float64,1}()
    push!(wpntsT,ocn.bndryPnts[mnY[2]].cntr.xy)
    push!(epntsT,ocn.bndryPnts[mnY[2]].cntr.xy)
    push!(wys, ocn.bndryPnts[mnY[2]].cntr.xy.y)
    push!(eys, ocn.bndryPnts[mnY[2]].cntr.xy.y)
    for xys in ocn.bndryPnts
        if xys.cntr.xy.x<(xys.cntr.xy.y*midLn[2]+midLn[1])
            push!(wpntsT, xys.cntr.xy)
            push!(wys, xys.cntr.xy.y)
        elseif xys.cntr.xy.x>(xys.cntr.xy.y*midLn[2]+midLn[1])
            push!(epntsT, xys.cntr.xy)
            push!(eys, xys.cntr.xy.y)
        else
        end
    end
    push!(wpntsT,ocn.bndryPnts[mxY[2]].cntr.xy)
    push!(epntsT,ocn.bndryPnts[mxY[2]].cntr.xy)
    push!(wys, ocn.bndryPnts[mxY[2]].cntr.xy.y)
    push!(eys, ocn.bndryPnts[mxY[2]].cntr.xy.y)

    #orders points ymin to ymax
    wpnts=Array{xy,1}()
    epnts=Array{xy,1}()
    for i =1:length(wys)
        index=findmin(wys)[2]
        wys[index]=Inf
        push!(wpnts,wpntsT[index])
    end
    for i =1:length(eys)
        index=findmin(eys)[2]
        eys[index]=Inf
        push!(epnts,epntsT[index])
    end

    #constructs line boundaries for limits
    wbnd=Array{line,1}()
    ebnd=Array{line,1}()
    for i = 1:length(wpnts)-1
        dummy_line=line()
        alpha_beta=reverse([[wpnts[i+1].y,wpnts[i].y] ones(2)]\[wpnts[i+1].x,wpnts[i].x])
        dummy_line.b=alpha_beta[1]
        dummy_line.m=alpha_beta[2]
        dummy_line.ymax=wpnts[i+1].y
        dummy_line.ymn=wpnts[i].y
        push!(wbnd,deepcopy(dummy_line))
    end
    for i = 1:length(epnts)-1
        dummy_line=line()
        alpha_beta=reverse([[epnts[i+1].y,epnts[i].y] ones(2)]\[epnts[i+1].x,epnts[i].x])
        dummy_line.b=alpha_beta[1]
        dummy_line.m=alpha_beta[2]
        dummy_line.ymax=epnts[i+1].y
        dummy_line.ymn=epnts[i].y
        push!(ebnd,deepcopy(dummy_line))
    end
    return wbnd,ebnd
end

function lof_getPccData(id_count)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "pcc_data")...)
    pccs=Array{bus,1}()
    for index=1:length(df[!, 1])
        dummy_bus=bus()
        dummy_bus.cntr.gps.lng=df.longitude[index]
        dummy_bus.cntr.gps.lat=df.latitude[index]
        dummy_bus.id=id_count
        push!(pccs,deepcopy(dummy_bus))
        pccs[length(pccs)].num=length(pccs)
        id_count=id_count+1
    end
    return pccs,id_count
end

function lof_getOwppData(id_count)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "owpp_data")...)
    owpps=Array{bus,1}()
    for index=1:length(df[!, 1])
        dummy_bus=bus()
        dummy_bus.cntr.gps.lng=df.longitude[index]
        dummy_bus.cntr.gps.lat=df.latitude[index]
        dummy_bus.id=id_count
        push!(dummy_bus.mvas,df.mva[index])
        push!(dummy_bus.wnds,df.name[index])
        push!(owpps,deepcopy(dummy_bus))
        owpps[length(owpps)].num=length(owpps)
        id_count=id_count+1
    end
    return owpps,id_count
end

function lof_getBndData()
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "domain_data")...)
    boundary=Array{node,1}()
    for index=1:length(df[!, 1])
        dummy_node=node()
        dummy_node.cntr.gps.lng=df.longitude[index]
        dummy_node.cntr.gps.lat=df.latitude[index]
        push!(boundary,deepcopy(dummy_node))
    end
    return boundary
end

function lof_getNoGoData(index)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "nogo_data"*string(index))...)
    nogo=Array{node,1}()
    for index=1:length(df[!, 1])
        dummy_node=node()
        dummy_node.cntr.gps.lng=df.longitude[index]
        dummy_node.cntr.gps.lat=df.latitude[index]
        push!(nogo,deepcopy(dummy_node))
    end
    return nogo
end

function lof_getSysData()
    sht = XLSX.readxlsx("layout//data.xlsx")["sys_data"]
    sys=system()
    sys.nogoNum=sht["B1"]
    sys.prec=sht["B2"]
    return sys
end
############################ GPS to cartesian transform ########################
################################################################################
#sets the gps coords that used as reference coords
function lof_bseCrd(ocean)
    base=gps()
    #for india (type layouts)
    if ocean.owpps[length(ocean.owpps)].cntr.gps.lat < ocean.pccs[length(ocean.pccs)].cntr.gps.lat
        base.lat=ocean.owpps[length(ocean.owpps)].cntr.gps.lat#base lat
        base.lng=ocean.owpps[length(ocean.owpps)].cntr.gps.lng#base long
    #for belgium (type layouts)
elseif ocean.pccs[length(ocean.pccs)].cntr.gps.lat < ocean.owpps[length(ocean.owpps)].cntr.gps.lat
        base.lat=ocean.pccs[length(ocean.pccs)].cntr.gps.lat#base lat
        base.lng=ocean.pccs[length(ocean.pccs)].cntr.gps.lng#base long
    else
        error("No proper base coordinates system established!")
    end
    return base
end

#calculates lengths based on latitude
#as lattitude changes number of km should be updated
function lof_gps2cartesian(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.cntr.xy.x=lof_deg2lgth(value.cntr.gps.lng-base.lng,lof_lg1deg(value.cntr.gps.lat,lnthLT))
        value.cntr.xy.y=lof_deg2lgth(value.cntr.gps.lat-base.lat,lnthLT)
    end
end

#rotates and slides cartesian axis
function lof_transformAxis(ocn)
    offset=lof_rotateAxis(ocn)
    lof_slideAxis(ocn,offset)
end

#finds angle to rotate and applies to owpps and pccs
#rotates axis to align n-s with y
function lof_rotateAxis(ocn)
    theta=atan((ocn.pccs[length(ocn.pccs)].cntr.xy.x-ocn.owpps[length(ocn.owpps)].cntr.xy.x)/(ocn.owpps[length(ocn.owpps)].cntr.xy.y-ocn.pccs[length(ocn.pccs)].cntr.xy.y))
    ocn.theta=theta
    offset=0.0
    offset=lof_rotateGroup_Os(ocn.owpps,theta,offset)
    offset=lof_rotateGroup_Os(ocn.pccs,theta,offset)
    ocn.offset=deepcopy(offset)
    lof_rotateGroup(ocn.bndryPnts,theta)
    for i=1:length(ocn.nogos)
        lof_rotateGroup(ocn.nogos[i].nodes,theta)
    end
    return ocn.offset
end

#loops through to apply rotations for a specified group
function lof_rotateGroup_Os(locations,theta,os)
    for value in locations
        xy=lof_rotatePnt(value.cntr.xy.x,value.cntr.xy.y,theta)
        value.cntr.xy.x=xy[1]
        value.cntr.xy.y=xy[2]
        if value.cntr.xy.x<os
            os=value.cntr.xy.x
        end
    end
    return os
end

#loops through to apply rotations for a specified group
function lof_rotateGroup(locations,theta)
    for value in locations
        xy=lof_rotatePnt(value.cntr.xy.x,value.cntr.xy.y,theta)
        value.cntr.xy.x=xy[1]
        value.cntr.xy.y=xy[2]
    end
end

#applies rotational matrix individual coordinates
function lof_rotatePnt(x,y,theta)
    co_od=[x y]
    rotated=co_od*[cos(theta) -1*sin(theta);sin(theta) cos(theta)]
    return rotated
end

#translates the entire region by specified offset
#sets unique IDs for owpps and pccs
function lof_slideAxis(ocn,os)
    for value in ocn.owpps
        value.cntr.xy.x=value.cntr.xy.x-os
    end
    for value in ocn.pccs
        value.cntr.xy.x=value.cntr.xy.x-os
    end
    for value in ocn.bndryPnts
        value.cntr.xy.x=value.cntr.xy.x-os
    end
    for i=1:length(ocn.nogos)
        for value in ocn.nogos[i].nodes
            value.cntr.xy.x=value.cntr.xy.x-os
        end
    end
end

#changes angle to an arc length
function lof_deg2lgth(d,dPl)
    return d*dPl
end

#calculates length of 1 deg of longitude at given lattitude
function lof_lg1deg(lat,lngth)
    return cos(lof_d2r(lat))*lngth
end

#Change radians to degrees
function lof_r2d(rad)
    return rad*180/pi
end

#Change degrees to radians
function lof_d2r(deg)
    return deg*pi/180
end
