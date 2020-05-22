using Geodesy, Combinatorics, LinearAlgebra
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
