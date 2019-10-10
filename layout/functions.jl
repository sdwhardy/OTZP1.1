function lof_layoutEez()
    #reading data from excel input files
    ocean=eez()
    ocean.id_count=0
    ocean.pccs,ocean.id_count=lof_getPccData(ocean.id_count)
    ocean.owpps,ocean.id_count=lof_getOwppData(ocean.id_count)
    ocean.bndryPnts=lof_getBndData()
    ocean.sys=lof_getSysData()
    ocean.finance=cstD_cfs()
    lof_nogoZones(ocean)

    #transformation from gps to cartesian
    ocean.base=lof_bseCrd(ocean)#find base coordinates
    print("Base coordinates are: ")
    println(ocean.base)
    lof_gps2cartesian4bus(ocean.owpps,ocean.base)#projects owpps onto cartesian plane
    lof_gps2cartesian4bus(ocean.pccs,ocean.base)#projects pccs onto cartesian plane
    lof_gps2cartesian4nodes(ocean.bndryPnts,ocean.base)#projects boundary points onto cartesian plane
    for i=1:ocean.sys.nogoNum
        lof_gps2cartesian4nodes(ocean.nogos[i].nodes,ocean.base)#projects boundary points onto cartesian plane
    end
    lof_transformAxis(ocean)

    #set area of owpp
    lof_setAreaOwpp(ocean)

    #set range of MV
    lof_mVrng(ocean)

    #build domain "mesh"
    lof_nodify(ocean)

    ocean.buses=vcat(ocean.pccs, ocean.owpps)#collects all buses

    ##########Printing
    for value in ocean.pccs
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    for value in ocean.owpps
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    println("GPS coordinates projected onto cartesian plane.")
    println("Axis transformed.")
    return ocean
end

function lof_nogoZones(ocn)
    for i=1:ocn.sys.nogoNum
        nogoT=nogo()
        nogoT.nodes=lof_getNoGoData(i)
        push!(ocn.nogos,deepcopy(nogoT))
    end
end

function lof_mVrng(ocn)
    for owpp in ocn.owpps
        owpp.mv_zone.radius=cstF_mVrng(owpp.zone.radius,owpp.mvas[1],owpp.wnds[1],ocn.finance,ocn.sys.prec)
        owpp.mv_zone.area=pi*(owpp.mv_zone.radius)^2
    end
end

function lof_setAreaOwpp(ocn)
    for owpp in ocn.owpps
        owpp.zone.area=owpp.mvas[1]/ocn.sys.mwPerKm
        owpp.zone.radius=sqrt(owpp.zone.area/pi)
    end
end

function lof_nodify(ocn)
    wbnd,ebnd=lof_bndPerimeter(ocn.bndryPnts)
    Allnodes=lof_addNodes(ocn,wbnd,ebnd)
    ocn.discretedom.nodes=lof_deleteNgNodes(ocn,Allnodes)


    lof_busPlaceOnNodes(ocn.owpps,ocn.discretedom.nodes)
    lof_busPlaceOnNodes(ocn.pccs,ocn.discretedom.nodes)
    #lof_nogoAreaNodes(ocn)
    lof_owppAreaNodes(ocn)
    lof_mvAreaNodes(ocn)
end

#=
function lof_nogoPerimeter(ngs)
    wbndNG=Array{line,1}()
    ebndNG=Array{line,1}()
    for i=1:length(ngs)
        wbnd_ng,ebnd_ng=lof_bndPerimeter(ngs[i].nodes)
        for ln in wbnd_ng
            push!(wbndNG,deepcopy(ln))
        end
        for ln in ebnd_ng
            push!(ebndNG,deepcopy(ln))
        end
    end
    return wbndNG, ebndNG
end
=#

function lof_nogoAreaNodes(ocn)
    for node in ocn.discretedom.nodes
        for owpp in ocn.owpps
            km2owpp=lof_pnt2pnt_dist(node.xy,owpp.node.xy)
            if km2owpp<=owpp.zone.radius
                push!(owpp.zone.pnts,node)
                if km2owpp+ocn.sys.prec>owpp.zone.radius
                    push!(owpp.zone.periPnts,node)
                end
            end
        end
    end
end


function lof_busPlaceOnNodes(buses,nodes)
    dummy_node=node()
    for bus in buses
        closestNode=Inf
        for node in nodes
            km2owpp=lof_pnt2pnt_dist(node.xy,bus.node.xy)
            if km2owpp<closestNode
                closestNode=deepcopy(km2owpp)
                dummy_node=node
            end
        end
        bus.node=deepcopy(dummy_node)
    end
end

function lof_owppAreaNodes(ocn)
    for node in ocn.discretedom.nodes
        for owpp in ocn.owpps
            km2owpp=lof_pnt2pnt_dist(node.xy,owpp.node.xy)
            if km2owpp<=owpp.zone.radius
                push!(owpp.zone.pnts,node)
                if km2owpp+ocn.sys.prec>owpp.zone.radius
                    push!(owpp.zone.periPnts,node)
                end
            end
        end
    end
end

function lof_mvAreaNodes(ocn)
    for node in ocn.discretedom.nodes
        for owpp in ocn.owpps
            km2owpp=lof_pnt2pnt_dist(node.xy,owpp.node.xy)
            if km2owpp<=owpp.mv_zone.radius
                push!(owpp.mv_zone.pnts,node)
                if km2owpp+ocn.sys.prec>owpp.mv_zone.radius
                    push!(owpp.mv_zone.periPnts,node)
                end
            end
        end
    end
end

function lof_deleteNgNodes(ocn,Allnodes)
    for ngs in ocn.nogos
        wbnd,ebnd=lof_bndPerimeter(ngs.nodes)
        ymx=wbnd[length(wbnd)].ymax
        ymn=wbnd[1].ymn
        windex=1
        eindex=1
        lnth=length(Allnodes)
        i=1
        while i<lnth
            if Allnodes[i].xy.y<=ymx && Allnodes[i].xy.y>=ymn
                while Allnodes[i].xy.y>wbnd[windex].ymax
                    windex=windex+1
                end
                while Allnodes[i].xy.y>ebnd[eindex].ymax
                    eindex=eindex+1
                end
                xmx=Allnodes[i].xy.y*ebnd[eindex].m+ebnd[eindex].b
                xmn=Allnodes[i].xy.y*wbnd[windex].m+wbnd[windex].b
                if Allnodes[i].xy.x>xmn && Allnodes[i].xy.x<xmx
                    deleteat!(Allnodes,i)
                    i=i-1
                    lnth=lnth-1
                else
                end
            else
            end
            i=i+1
        end
    end
    return deepcopy(Allnodes)
end

function lof_addNodes(ocn,wbnd,ebnd)
    dummy_nodes=Array{node,1}()
    ymx=wbnd[length(wbnd)].ymax
    ymn=wbnd[1].ymn
    prec=ocn.sys.prec
    windex=1
    eindex=1
    for y=ymn:prec:ymx
        while y>wbnd[windex].ymax
            windex=windex+1
        end
        while y>ebnd[eindex].ymax
            eindex=eindex+1
        end
        xmx=y*ebnd[eindex].m+ebnd[eindex].b
        xmn=y*wbnd[windex].m+wbnd[windex].b
        x=xmn
        while (x>=xmn && x<=xmx)
            dummy_node=node()
            dummy_node.xy.x=x
            dummy_node.xy.y=y
            push!(dummy_nodes,deepcopy(dummy_node))
            x=x+prec
        end
    end
    return dummy_nodes
end

function lof_bndPerimeter(bnd)
    #finds mid line dividing east and west regions
    ys=Array{Float64,1}()
    for xys in bnd
        push!(ys,xys.xy.y)
    end
    mxY=findmax(ys)
    mnY=findmin(ys)
    midLn=reverse([[bnd[mxY[2]].xy.y,bnd[mnY[2]].xy.y] ones(2)]\[bnd[mxY[2]].xy.x,bnd[mnY[2]].xy.x])

    #sorts to east and west bounding points
    wpntsT=Array{xy,1}()
    epntsT=Array{xy,1}()
    wys=Array{Float64,1}()
    eys=Array{Float64,1}()
    push!(wpntsT,bnd[mnY[2]].xy)
    push!(epntsT,bnd[mnY[2]].xy)
    push!(wys, bnd[mnY[2]].xy.y)
    push!(eys, bnd[mnY[2]].xy.y)
    for xys in bnd
        if xys.xy.x<(xys.xy.y*midLn[2]+midLn[1])
            push!(wpntsT, xys.xy)
            push!(wys, xys.xy.y)
        elseif xys.xy.x>(xys.xy.y*midLn[2]+midLn[1])
            push!(epntsT, xys.xy)
            push!(eys, xys.xy.y)
        else
        end
    end
    push!(wpntsT,bnd[mxY[2]].xy)
    push!(epntsT,bnd[mxY[2]].xy)
    push!(wys, bnd[mxY[2]].xy.y)
    push!(eys, bnd[mxY[2]].xy.y)

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
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
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
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.name=df.name[index]
        dummy_bus.id=id_count
        push!(dummy_bus.mvas,df.mva[index])
        push!(owpps,deepcopy(dummy_bus))
        owpps[length(owpps)].num=length(owpps)
        id_count=id_count+1
    end
    lof_getOwppWindData(owpps)
    return owpps,id_count
end

function lof_getOwppWindData(owpps)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "wind_data")...)
    for owpp in owpps
        wnd=wndF_wndPrf([getproperty(df,Symbol(owpp.name))])
        push!(owpp.wnds,deepcopy(wnd))
    end
end

function lof_getBndData()
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "domain_data")...)
    boundary=Array{node,1}()
    for index=1:length(df[!, 1])
        dummy_node=node()
        dummy_node.gps.lng=df.longitude[index]
        dummy_node.gps.lat=df.latitude[index]
        push!(boundary,deepcopy(dummy_node))
    end
    return boundary
end

function lof_getNoGoData(index)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "nogo_data"*string(index))...)
    nogo=Array{node,1}()
    for index=1:length(df[!, 1])
        dummy_node=node()
        dummy_node.gps.lng=df.longitude[index]
        dummy_node.gps.lat=df.latitude[index]
        push!(nogo,deepcopy(dummy_node))
    end
    return nogo
end

function lof_getSysData()
    sht = XLSX.readxlsx("layout//data.xlsx")["sys_data"]
    sys=system()
    sys.nogoNum=sht["B1"]
    sys.prec=sht["B2"]
    sys.mwPerKm=sht["B3"]
    return sys
end
############################ GPS to cartesian transform ########################
################################################################################
#sets the gps coords that used as reference coords
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

#calculates lengths based on latitude
#as lattitude changes number of km should be updated
function lof_gps2cartesian4nodes(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.xy.x=lof_deg2lgth(value.gps.lng-base.lng,lof_lg1deg(value.gps.lat,lnthLT))
        value.xy.y=lof_deg2lgth(value.gps.lat-base.lat,lnthLT)
    end
end

#as lattitude changes number of km should be updated
function lof_gps2cartesian4bus(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.node.xy.x=lof_deg2lgth(value.node.gps.lng-base.lng,lof_lg1deg(value.node.gps.lat,lnthLT))
        value.node.xy.y=lof_deg2lgth(value.node.gps.lat-base.lat,lnthLT)
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
    theta=atan((ocn.pccs[length(ocn.pccs)].node.xy.x-ocn.owpps[length(ocn.owpps)].node.xy.x)/(ocn.owpps[length(ocn.owpps)].node.xy.y-ocn.pccs[length(ocn.pccs)].node.xy.y))
    ocn.theta=theta
    offset=0.0
    offset=lof_rotateGroup4bus(ocn.owpps,theta,offset)
    offset=lof_rotateGroup4bus(ocn.pccs,theta,offset)
    ocn.offset=deepcopy(offset)
    lof_rotateGroup4node(ocn.bndryPnts,theta)
    for i=1:length(ocn.nogos)
        lof_rotateGroup4node(ocn.nogos[i].nodes,theta)
    end
    return ocn.offset
end

#loops through to apply rotations for a specified group
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

#loops through to apply rotations for a specified group
function lof_rotateGroup4node(locations,theta)
    for value in locations
        xy=lof_rotatePnt(value.xy.x,value.xy.y,theta)
        value.xy.x=xy[1]
        value.xy.y=xy[2]
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
        value.node.xy.x=value.node.xy.x-os
    end
    for value in ocn.pccs
        value.node.xy.x=value.node.xy.x-os
    end
    for value in ocn.bndryPnts
        value.xy.x=value.xy.x-os
    end
    for i=1:length(ocn.nogos)
        for value in ocn.nogos[i].nodes
            value.xy.x=value.xy.x-os
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


#returns the hypotenuse distance between 2 cartesian points
#minimum distance for a path is 1km
function lof_pnt2pnt_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
    if hyp < 1
        hyp=1
        #println("Arc distance is less than 1km, set to 1km.")
    end
    return hyp
end
