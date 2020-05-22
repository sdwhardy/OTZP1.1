include("functions.jl")
#This file contains the exported functions
################################################################################
########################### Finding MV max distance ############################
################################################################################
function find_max_mv_transmission(owpps,database)
    ks=get_Cost_Data()
    for owpp in owpps
        pnt0=0.001
        pnt1=deepcopy(database["cables"]["66.0"][string(owpp.mva)][length(database["cables"]["66.0"][string(owpp.mva)])].length)
        mv_cable_pnt0=deepcopy(mvac_cable(owpp.mva,pnt0,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
        mv_cable_pnt1=deepcopy(mvac_cable(owpp.mva,pnt1,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
        plat=platform()
        plat.mva=owpp.mva
        plat.wnd=owpp.wnd
        plat=cost_ac_platform(plat,ks)
        plat=adjust_base_ac_platform(plat,ks)
        np=0
        if (database["bits"]["hvac"]==true || database["bits"]["mpc_ac"]==true)
            hv_c="220.0"
        elseif (database["bits"]["hvac"]==true)
            hv_c="300.0"
        end
        hv_cable_pnt0=deepcopy(hvac_cable(owpp.mva,pnt0,owpp.wnd,database["cables"][hv_c][string(owpp.mva)],ks))
        if (mv_cable_pnt0.costs.ttl<hv_cable_pnt0.costs.ttl+plat.costs.ttl)
            hv_cable_pnt1=deepcopy(hvac_cable(owpp.mva,pnt1,owpp.wnd,database["cables"][hv_c][string(owpp.mva)],ks))
            np=pnt1
            if (mv_cable_pnt1.costs.ttl>hv_cable_pnt1.costs.ttl+plat.costs.ttl)
                while (pnt1-pnt0>1)
                    np=(pnt0+pnt1)/2
                    mv_cable_np=deepcopy(mvac_cable(owpp.mva,np,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
                    hv_cable_np=deepcopy(hvac_cable(owpp.mva,np,owpp.wnd,database["cables"][hv_c][string(owpp.mva)],ks))
                    if (mv_cable_np.costs.ttl<hv_cable_np.costs.ttl+plat.costs.ttl)
                        pnt0=deepcopy(np)
                    elseif (mv_cable_np.costs.ttl>hv_cable_np.costs.ttl+plat.costs.ttl)
                        pnt1=deepcopy(np)
                    else
                        break
                    end
                end
            end
        end
        owpp.mv_zone=deepcopy(np)
        owpp.kV=66.0
    end
    return owpps
end

################################################################################
############################## Finding lookup table ############################
################################################################################
function get_equipment_tables(km_mva_set,wnd,kv)
    cable_database=get_cable_table(km_mva_set,wnd)
    xfo_database,conv_database=get_xfoConv_table(km_mva_set,wnd)
    database=Dict{String, Dict{String,Any}}()
    push!(database,"cables"=>cable_database)
    push!(database,"transformers"=>xfo_database)
    push!(database,"converters"=>conv_database)

    #vheck if hvdc is to be included
    bits=Dict([("hvdc", true), ("mpc_ac", false), ("hvac", true)])
    push!(database,"bits"=>bits)
    hvdc_bit=check4hvdc(database,km_mva_set,wnd,kv)

    #check hvac connections
    hvac_bit=true
    if (hvdc_bit==true)
        hvac_bit=check4hvac(database,km_mva_set,wnd,kv)
    end

    #vheck if mid point compensation is to be included
    bits=Dict([("hvdc", false), ("mpc_ac", true), ("hvac", true)])
    push!(database,"bits"=>bits)
    mpc_bit=check4mpc_ac(database,km_mva_set,wnd,kv)

    #final bits sets
    bits=Dict([("hvdc", hvdc_bit), ("mpc_ac", mpc_bit), ("hvac", hvac_bit)])
    push!(database,"bits"=>bits)

    return database
end

################################################################################
######################## Optimal equipment using lookup table ##################
################################################################################
#finds the offshore transformer using the look up table
function xfo_oss(xfo0,ks,xfo_data)
    xfo=transformer()
    xfo.costs.ttl=Inf
    xfo.mva=xfo0.mva
    for xd in xfo_data
        if ((xd+10<=xfo0.mva) && (3*xd>xfo0.mva))
            xfo0.num=1
            xfo0.elec.mva=xd
            while ((xfo0.num*xfo0.elec.mva)<xfo0.mva)
                xfo0.num=xfo0.num+1
            end
            xfo0=cost_xfo_oss(xfo0,ks)
            if (xfo0.costs.ttl<xfo.costs.ttl)
                xfo=deepcopy(xfo0)
            end
        end
    end
    return xfo
end

#finds the onshore transformer using the look up table
function xfo_pcc(xfo0,ks,xfo_data)
    xfo=transformer()
    xfo.costs.ttl=Inf
    xfo.mva=xfo0.mva
    for xd in xfo_data
        if ((xd+10<=xfo0.mva) && (3*xd>xfo0.mva))
            xfo0.num=1
            xfo0.elec.mva=xd
            while ((xfo0.num*xfo0.elec.mva)<xfo0.mva)
                xfo0.num=xfo0.num+1
            end
            xfo0=cost_xfo_pcc(xfo0,ks)
            if (xfo0.costs.ttl<xfo.costs.ttl)
                xfo=deepcopy(xfo0)
            end
        end
    end
    return xfo
end


#finds the hvdc cable using the look up table
function hvdc_cable(mva,km,wnd,cable_array,ks)
    cbl=deepcopy(cable_array[1])
    cbl.wnd=wnd
    cbl.length=km
    cbl=cost_hvdc_cable(cbl,ks)
    return cbl
end

#=cable_220=deepcopy(hvac_cable(mva,190.0,wnd,database["cables"]["220.0"][string(mva)],ks))
cbl=deepcopy(database["cables"]["220.0"][string(mva)][6])=#
#finds the hvac cable using the look up table
function hvac_cable(mva,km,wnd,cable_array,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    for cbl0 in cable_array
        if (km<=cbl0.length+0.01)#+0.01 is just to be certain floating point error does not occur
            cbl=deepcopy(cbl0)
            break
        end
    end
    cbl.wnd=wnd
    cbl.length=km
    cbl.elec.mva=get_newQ_Capacity(cbl.elec.freq,km,cbl.elec.volt,cbl.elec.farrad,cbl.elec.amp)
    cbl=cost_hvac_cable(cbl,ks)
    return cbl
end

#finds the mvac cable using the look up table
function mvac_cable(mva,km,wnd,cable_array,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    for cbl0 in cable_array
        if (km<=cbl0.length+0.001)#+0.01 is just to be certain floating point error does not occur
            cbl=deepcopy(cbl0)
            break
        end
    end
    cbl.elec.mva=get_newQ_Capacity(cbl.elec.freq,km,cbl.elec.volt,cbl.elec.farrad,cbl.elec.amp)
    cbl.wnd=wnd
    cbl.length=km
    cbl=cost_mvac_cable(cbl,ks)
    return cbl
end

#Converter
#calculates the cost of an offshore HVDC converter
#each line needing a converter shares the total before the fixed cost is added at the end
function cost_hvdc_oss(conv,ks)
    #capex converter
    conv.costs.cpx=capex_hvdc(conv,ks)
    #cost of losses
    conv.costs.tlc=cost_tlc(conv,ks)
    #corrective maintenance
    conv.costs.cm=cost_cm(conv.costs.cpx,ks.opx_co)
    #eens calculation
    conv.costs.eens=cost_eens(conv,ks)
    #totals the xfo cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

#once all incoming lines with converters are combined the fixed cost is added
function adjust_base_hvdc_offshore_converter(conv,ks)
    #adjust to add base cost capex converter
    conv.costs.cpx=conv.costs.cpx+ks.conv_d
    #corrective maintenance
    conv.costs.cm=conv.costs.cm+ks.conv_d*ks.opx_co*npv_years()
    #totals the converter cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

#calculates the cost of an onshore HVDC converter
#each line needing a converter shares the total before the fixed cost is added at the end
function cost_hvdc_pcc(conv,ks)
    #capex converter
    conv.costs.cpx=capex_hvdc(conv,ks)
    #cost of losses
    conv.costs.tlc=cost_tlc(conv,ks)
    #corrective maintenance
    conv.costs.cm=cost_cm(conv.costs.cpx,ks.opx_cp)
    #eens calculation
    conv.costs.eens=cost_eens(conv,ks)
    #totals the xfo cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

#once all incoming lines with converters are combined the fixed cost is added
function adjust_base_hvdc_onshore_converter(conv,ks)
    #adjust to add base cost capex converter
    conv.costs.cpx=conv.costs.cpx+ks.conv_d
    #corrective maintenance
    conv.costs.cm=conv.costs.cm+ks.conv_d*ks.opx_cp*npv_years()
    #totals the converter cost
    conv.costs.ttl=cost_conv_sum(conv)
    return conv
end

################################### platform ###################################
#calculates incremental cost of an AC platform NOT including fixed cost
function cost_ac_platform(plat,ks)
    #capex oss
    plat.costs.cpx=capex_plat_ac(plat,ks)
    #corrective maintenance
    plat.costs.cm=cost_cm(plat.costs.cpx,ks.opx_pl)
    #totals the plat cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end

#once all incoming lines established the fixed cost is added to the of incremental sum - AC
function adjust_base_ac_platform(plat,ks)
    #adjust to add base cost capex converter
    plat.costs.cpx=plat.costs.cpx+ks.pac_f
    #corrective maintenance
    plat.costs.cm=plat.costs.cm+ks.pac_f*ks.opx_pl*npv_years()
    #totals the converter cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end

#calculates incremental cost of an DC platform NOT including fixed cost
function cost_dc_platform(plat,ks)
    #capex oss
    plat.costs.cpx=capex_plat_dc(plat,ks)
    #corrective maintenance
    plat.costs.cm=cost_cm(plat.costs.cpx,ks.opx_pl)
    #totals the plat cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end

#once all incoming lines established the fixed cost is added to the of incremental sum - DC
function adjust_base_dc_platform(plat,ks)
    #adjust to add base cost capex converter
    plat.costs.cpx=plat.costs.cpx+ks.pdc_h
    #corrective maintenance
    plat.costs.cm=plat.costs.cm+ks.pdc_h*ks.opx_pl*npv_years()
    #totals the converter cost
    plat.costs.ttl=plat.costs.cpx+plat.costs.cm
    return plat
end
