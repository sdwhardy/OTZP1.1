########################### OWPP to PCC direct connection ######################
#finds the optimal connection for OWPP to PCC
function optimal_owpp2pcc(circ,ocn)

    #MV connection
    ks=get_Cost_Data()
    km=euclidian_distance(circ.owpps[1].node.xy,circ.pcc.node.xy)
    mv_connection=Dict{String,Any}()
    push!(mv_connection,("cost"=>Inf))#xfo_for
    if (circ.owpps[1].mv_zone>km)
        #cable
        mv_cable=deepcopy(mvac_cable(circ.owpps[1].mva,km,circ.owpps[1].wnd,ocn.database["cables"]["66.0"][string(circ.owpps[1].mva)],ks))
        mv_cable.costs.grand_ttl=mv_cable.costs.ttl
        #transformer
        xfo=deepcopy(ocn.database["transformers"][string(circ.owpps[1].mva)]["onshore"])
        xfo.wnd=circ.owpps[1].wnd
        xfo.elec.hv=circ.pcc.kV
        xfo.elec.lv=mv_cable.elec.volt
        xfo=cost_xfo_pcc(xfo,ks)
        #totals and organizes mv equipment
        mv_connection_cost=mv_cable.costs.ttl+xfo.costs.ttl
        push!(mv_connection,("cable"=>mv_cable))#cable
        push!(mv_connection,("xfo_for"=>xfo))#xfo_for
        push!(mv_connection,("cost"=>mv_connection_cost))#xfo_for
        push!(mv_connection["cable"].path,circ.owpps[1].node)
        push!(mv_connection["cable"].path,circ.pcc.node)
    end
    #finds HV alternative and choses cheaper option
    hv_connection=optimal_mog2pcc(km,circ.owpps[1].mva,circ.pcc.kV,circ.owpps[1].wnd,ocn.database)
    owpp2pcc_connection_cost=findmin([mv_connection["cost"],hv_connection["cost"]])
    owpp2pcc_connection=[mv_connection,hv_connection][owpp2pcc_connection_cost[2]]
    #Build the zero length MVAC cable
    mv_cable0=deepcopy(mvac_cable(circ.owpps[1].mva,1,circ.owpps[1].wnd,ocn.database["cables"]["66.0"][string(circ.owpps[1].mva)],ks))
    mv_cable0.length=0.0
    if (mv_cable0.costs.ttl==Inf)
        mv_cable0.costs.perkm_ttl=10^9
    end
    mv_cable0.costs.grand_ttl=0.0
    push!(mv_cable0.path,circ.owpps[1].node)
    #Build the best solution circuit
    if (haskey(owpp2pcc_connection,"plat_aft_dc"))#HVDC connection
        #build the offshore MOG
        push!(circ.mog,aft_connection_DC_equipment(owpp2pcc_connection,circ.owpps[1],ocn))
        #lay the cable
        push!(mv_cable0.path,circ.mog[1].node)
        push!(circ.MVcbls,mv_cable0)
        push!(owpp2pcc_connection["cable"].path,circ.mog[1].node)
        push!(owpp2pcc_connection["cable"].path,circ.pcc.node)
        push!(circ.PCCcbls,owpp2pcc_connection["cable"])
        #add the converter onshore
        push!(circ.pcc.conv,owpp2pcc_connection["conv_for"])
        #add the transformer onshore
        push!(circ.pcc.xfmrs,owpp2pcc_connection["xfo_for"])
        ocn.owpps[circ.base.node.num].kv2pcc=owpp2pcc_connection["cable"].elec.volt
    elseif (haskey(owpp2pcc_connection,"plat_aft_ac"))#HVAC
        #Build MOG
        push!(circ.mog,aft_connection_AC_equipment(owpp2pcc_connection,circ.owpps[1],ocn))
        #lay the cable
        push!(mv_cable0.path,circ.mog[1].node)
        push!(circ.MVcbls,mv_cable0)
        push!(owpp2pcc_connection["cable"].path,circ.mog[1].node)
        push!(owpp2pcc_connection["cable"].path,circ.pcc.node)
        push!(circ.PCCcbls,owpp2pcc_connection["cable"])
        #add the transformer onshore
        push!(circ.pcc.xfmrs,owpp2pcc_connection["xfo_for"])
        ocn.owpps[circ.base.node.num].kv2pcc=owpp2pcc_connection["cable"].elec.volt
    else#MVAC
        #lay the cable
        push!(circ.MVcbls,owpp2pcc_connection["cable"])
        #add the transformer onshore
        push!(circ.pcc.xfmrs,owpp2pcc_connection["xfo_for"])
        ocn.owpps[circ.base.node.num].kv2pcc=owpp2pcc_connection["cable"].elec.volt
    end
    circ=total_circuit_cost(circ)
    return circ
end

function aft_connection_AC_equipment(connection,connect_bus,ocn)
    #create an OSS bus
    oss=bus()
    oss.kV=connection["cable"].elec.volt
    oss.mva=connection["cable"].mva
    oss.wnd=connection["cable"].wnd
    oss.node=deepcopy(connect_bus.node)
    ocn.num=ocn.num+1
    oss.node.num=deepcopy(ocn.num)
    #build the platform
    push!(oss.plat,connection["plat_aft_ac"])
    #place the transformers
    push!(oss.xfmrs,connection["xfo_aft"])
    return oss
end

function aft_connection_DC_equipment(connection,connect_bus,ocn)
    oss=bus()
    oss.kV=connection["cable"].elec.volt
    oss.mva=connection["cable"].mva
    oss.wnd=connection["cable"].wnd
    oss.node=deepcopy(connect_bus.node)
    ocn.num=ocn.num+1
    oss.node.num=deepcopy(ocn.num)
    #build the platform
    push!(oss.plat,connection["plat_aft_ac"])
    push!(oss.plat,connection["plat_aft_dc"])
    #place the transformers
    push!(oss.xfmrs,connection["xfo_aft"])
    #place the transformers
    push!(oss.conv,connection["conv_aft"])
    return oss
end


#returns the final cost of the circuit
function total_circuit_cost(circ)
    circ.cost=0.0
    #totla cable caosts
    for mv_cbl in circ.MVcbls
        circ.cost=circ.cost+mv_cbl.costs.grand_ttl
    end
    for hv_cbl in circ.HVcbls
        circ.cost=circ.cost+hv_cbl.costs.grand_ttl
    end
    for hv_cbl in circ.O2Ocbls
        circ.cost=circ.cost+hv_cbl.costs.grand_ttl
    end
    for hv_cbl in circ.PCCcbls
        circ.cost=circ.cost+hv_cbl.costs.grand_ttl
    end

    #total transformer and converter costs
    for ss in circ.mog
        for xfo in ss.xfmrs
            circ.cost=circ.cost+xfo.costs.ttl
        end
        for cv in ss.conv
            circ.cost=circ.cost+cv.costs.ttl
        end
        #platforms per ss
        for plt in ss.plat
            circ.cost=circ.cost+plt.costs.ttl
        end
    end
    #totals transformers and converters and platform costs
    for ss in circ.oss
        #transformers per ss
        for xfo in ss.xfmrs
            circ.cost=circ.cost+xfo.costs.ttl
        end
        #converters per ss
        for cv in ss.conv
            circ.cost=circ.cost+cv.costs.ttl
        end
        #platforms per ss
        for plt in ss.plat
            circ.cost=circ.cost+plt.costs.ttl
        end
    end
    #transformers per pcc
    for xfo in circ.pcc.xfmrs
        circ.cost=circ.cost+xfo.costs.ttl
    end
    #converters per pcc
    for cv in circ.pcc.conv
        circ.cost=circ.cost+cv.costs.ttl
    end
    return circ
end



######################### Main feeder connections MOG to PCC ###################
#finds the optimal connection for the main MOG to PCC connection
function optimal_mog2pcc(km,mva,kv,wnd,database)
    ks=get_Cost_Data()
    #HVAC connection - mpc handled within HVAC
    ac_connection=hvac_canditate_mog2pcc(km,mva,kv,wnd,database,ks)#returns dictionary: [platform,trans_aft,cable,trans_for]
    #HVDC connection
    dc_connection=hvdc_canditate_mog2pcc(km,mva,kv,wnd,database,ks)#returns dictionary: [platform,trans_aft,conv_aft,cable,conv_for,trans_for]

    lowest_cost=findmin([dc_connection["cost"],ac_connection["cost"]])
    optimal_equipment=[dc_connection,ac_connection][lowest_cost[2]]
    return optimal_equipment
end

#HVAC connection from an MOG to a PCC
function hvac_canditate_mog2pcc(km,mva,kv_forward,wnd,database,ks)
    #initialize dummy set
    dict_ac_connection=Dict{String,Any}()
    push!(dict_ac_connection,("cost"=>Inf))#xfo_for
    #if HVAC is selected find optimal equipments
    if (database["bits"]["hvac"]==true)
        #build candidate cables
        cable_220=deepcopy(hvac_cable(mva,km,wnd,database["cables"]["220.0"][string(mva)],ks))
        cable_400=deepcopy(hvac_cable(mva,km,wnd,database["cables"]["400.0"][string(mva)],ks))
        cable_220.costs.grand_ttl=cable_220.costs.ttl
        cable_400.costs.grand_ttl=cable_400.costs.ttl
        #if mid point compensation is checked look for cheaper cable options
        if (database["bits"]["mpc_ac"]==true)
            println("jello from the inside!")
            cable_220,cable_400=check4beter_mpc_ac_cables(mva,km,database,wnd,cable_220,cable_400,ks)
        end
        #aft candidate transformers
        xfo_aft=deepcopy(database["transformers"][string(mva)]["offshore"])
        xfo_aft.wnd=wnd
        xfo_aft=cost_xfo_oss(xfo_aft,ks)

        #platform
        plat_aft=platform()
        plat_aft.acdc="ac"
        plat_aft.mva=mva
        plat_aft.wnd=wnd
        plat_aft=cost_ac_platform(plat_aft,ks)
        plat_aft=adjust_base_ac_platform(plat_aft,ks)

        #onshore transformera
        if (kv_forward!=400.0)
            xfo_for400=deepcopy(database["transformers"][string(mva)]["onshore"])
            xfo_for400.wnd=wnd
            xfo_for400.elec.hv=220.0
            xfo_for400.elec.lv=400.0
            xfo_for400=cost_xfo_pcc(xfo_for400,ks)
        else
            xfo_for400=transformer()
            xfo_for400.costs.ttl=0
        end
        #onshore
        if (kv_forward!=220.0)
            xfo_for220=deepcopy(database["transformers"][string(mva)]["onshore"])
            xfo_for220.wnd=wnd
            xfo_for220.elec.hv=400.0
            xfo_for220.elec.lv=220.0
            xfo_for220=cost_xfo_pcc(xfo_for220,ks)
        else
            xfo_for220=transformer()
            xfo_for220.costs.ttl=0
        end
        #find best alternative PCC connection voltage
        ac_connection220_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_220.costs.grand_ttl+xfo_for220.costs.ttl
        ac_connection400_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_400.costs.grand_ttl+xfo_for400.costs.ttl
        best_voltage=findmin([ac_connection220_cost,ac_connection400_cost])
        ac_connection220=[plat_aft,xfo_aft,cable_220,xfo_for220]
        ac_connection400=[plat_aft,xfo_aft,cable_400,xfo_for400]
        ac_connection=[ac_connection220,ac_connection400][best_voltage[2]]
        cost_ac_connection=best_voltage[1]
        #store object in a dictionary
        push!(dict_ac_connection,("plat_aft_ac"=>ac_connection[1]))#plat_aft
        push!(dict_ac_connection,("xfo_aft"=>ac_connection[2]))#xfo_aft
        push!(dict_ac_connection,("cable"=>ac_connection[3]))#cable
        push!(dict_ac_connection,("xfo_for"=>ac_connection[4]))#xfo_for
        push!(dict_ac_connection,("cost"=>cost_ac_connection))#xfo_for
        dict_ac_connection["plat_aft_ac"].kv=dict_ac_connection["cable"].elec.volt
        dict_ac_connection["xfo_aft"].elec.hv=dict_ac_connection["cable"].elec.volt
    end
    return dict_ac_connection
end

#check if HVAC mid-point compensation is better option
function check4beter_mpc_ac_cables(mva,km,database,wnd,cable_220,cable_400,ks)
    cable_220mpc=deepcopy(hvac_cable(mva,km/2,wnd,database["cables"]["220.0"][string(mva)],ks))
    cable_400mpc=deepcopy(hvac_cable(mva,km/2,wnd,database["cables"]["400.0"][string(mva)],ks))
    cable_220mpc.mpc_ac=true
    cable_400mpc.mpc_ac=true
    #platform
    #compensation
    cable_220mpc.plat.mva=mva
    cable_220mpc.plat.wnd=wnd
    cable_220mpc.plat=cost_ac_platform(cable_220mpc.plat,ks)
    cable_220mpc.plat=adjust_base_ac_platform(cable_220mpc.plat,ks)
    cable_400mpc.plat=cable_220mpc.plat
    #total cost
    cable_220mpc.costs.grand_ttl=cable_220mpc.costs.ttl*2+cable_220mpc.plat.costs.ttl
    cable_400mpc.costs.grand_ttl=cable_400mpc.costs.ttl*2+cable_400mpc.plat.costs.ttl
    #compare with no compensation
    if (cable_220mpc.costs.grand_ttl<cable_220.costs.grand_ttl)
        cable_220=deepcopy(cable_220mpc)
    end
    if (cable_400mpc.costs.grand_ttl<cable_400.costs.grand_ttl)
        cable_400=deepcopy(cable_400mpc)
    end
    #return best option
    println("cable size: "*string(cable_220.size))
    println("cable num: "*string(cable_220.num))
    println("cable length: "*string(cable_220.length))
    return cable_220,cable_400
end

#hvdc mog to pcc candidate connection
function hvdc_canditate_mog2pcc(km,mva,kv_forward,wnd,database,ks)
    dict_dc_connection=Dict{String,Any}()
    push!(dict_dc_connection,("cost"=>Inf))#conv_for
    if (database["bits"]["hvdc"]==true)
        cable_dc=deepcopy(hvdc_cable(mva,km,wnd,database["cables"]["300.0"][string(mva)],ks))
        cable_dc.costs.grand_ttl=cable_dc.costs.ttl
        #aft
        #converters
        conv_aft=deepcopy(database["converters"][string(mva)]["offshore"])
        conv_aft.wnd=wnd
        conv_aft=cost_hvdc_oss(conv_aft,ks)
        conv_aft=adjust_base_hvdc_offshore_converter(conv_aft,ks)

        #transformers
        xfo_aft=transformer()
        xfo_aft=deepcopy(database["transformers"][string(mva)]["offshore"])
        xfo_aft.wnd=wnd
        xfo_aft=cost_xfo_oss(xfo_aft,ks)

        #base ac platform
        plat_aft_ac=platform()
        plat_aft_ac.acdc="ac"
        plat_aft_ac.mva=mva
        plat_aft_ac.wnd=wnd
        plat_aft_ac=cost_ac_platform(plat_aft_ac,ks)
        #hvdc portion of platform
        plat_aft_dc=platform()
        plat_aft_dc.acdc="dc"
        plat_aft_dc.mva=mva
        plat_aft_dc.wnd=wnd
        plat_aft_dc=cost_dc_platform(plat_aft_dc,ks)
        plat_aft_dc=adjust_base_dc_platform(plat_aft_dc,ks)
        #conv_aft=hvdc_offshore_station_cost(conv_aft)

        #forward
        #converter
        conv_for=deepcopy(database["converters"][string(mva)]["onshore"])
        conv_for.wnd=wnd
        conv_for=cost_hvdc_pcc(conv_for,ks)
        conv_for=adjust_base_hvdc_onshore_converter(conv_for,ks)
        #transformers
        xfo_for=deepcopy(database["transformers"][string(mva)]["onshore"])
        xfo_for.wnd=wnd
        xfo_for=cost_xfo_pcc(xfo_for,ks)
        #conv_for=hvdc_onshore_station_cost(conv_for)

        #total cost and store in a dictionary
        cost_dc_connection=plat_aft_ac.costs.ttl+plat_aft_dc.costs.ttl+xfo_aft.costs.ttl+conv_aft.costs.ttl+cable_dc.costs.grand_ttl+conv_for.costs.ttl+xfo_for.costs.ttl
        push!(dict_dc_connection,("plat_aft_ac"=>plat_aft_ac))#plat_aft
        push!(dict_dc_connection,("plat_aft_dc"=>plat_aft_dc))#plat_aft
        push!(dict_dc_connection,("xfo_aft"=>xfo_aft))#conv_aft
        push!(dict_dc_connection,("conv_aft"=>conv_aft))#conv_aft
        push!(dict_dc_connection,("cable"=>cable_dc))#cable
        push!(dict_dc_connection,("conv_for"=>conv_for))#conv_for
        push!(dict_dc_connection,("xfo_for"=>xfo_for))#conv_for
        push!(dict_dc_connection,("cost"=>cost_dc_connection))#conv_for
    end

    return dict_dc_connection
end
######################### NOTE: tested to this point OK!

########################### OWPP to OSS connection #############################
#finds the optimal connections for OWPPS an OSS
#This is circuit building
#this is a layout function
#this is only a sorting function - seperates scenarios AC-AC, AC-DC, DC-DC, DC-AC
#NOTE this si now to be tested, beware there will be many mistakes testing needs to be 100% before proceeding!!!!
function optimal_owpps2oss(mog_circ,database)
    ks=get_Cost_Data()
    owpp2oss_connections=Array{Dict{String,Dict{String,Any}},1}()
    #sort OWPP connections to be made from farthest away to closest
    lengths_owpps=Array{Tuple,1}()
    hvdc_check=false
    for owp in mog_circ.owpps
        push!(lengths_owpps,(euclidian_distance(owp.node.xy,mog_circ.mog[1].node.xy),owp))
        if (owp.kv2pcc==300.0 && database["bits"]["hvdc"]==true)
            hvdc_check=true
        end
    end
    if (hvdc_check==true && mog_circ.mog[1].kV!=300.0 && database["bits"]["hvdc"]==true)
        println("optimal_owpps2oss says: Feeder upgrade to HVDC may be needed! Circuit ID: "*string(id))
    end
    sort!(lengths_owpps,by = x-> x[1],rev=true)

    if (oss.kV==300.0)
        #HVDC Feeder
        #AC to DC options handled
        if (hvdc_check==false)
            for l_o in lengths_owpps
                push!(owpp2oss_connections,deepcopy(connection_ACowpp2ossACDC(l_o[1],l_o[2],mog_circ.mog[1],database,ks)))
            end
            mog_circ=combine_connections_ACowpp2ossDC(owpp2oss_connections,mog_circ,database)
        #DC to DC options handled here
        else
            for l_o in lengths_owpps
                push!(owpp2oss_connections,deepcopy(optimal_ACDCowpp2ossDC(l_o[1],l_o[2],mog_circ.mog[1],database,ks)))
            end
            mog_circ=combine_connections_ACDCowpp2ossDC(owpp2oss_connections,mog_circ,database)
        end
    else
        #AC to AC option handled here
        for l_o in lengths_owpps
            push!(owpp2oss_connections,deepcopy(connection_ACowpp2ossACDC(l_o[1],l_o[2],mog_circ.mog[1],database,ks)))
        end
        mog_circ=combine_connections_ACowpp2ossAC(owpp2oss_connections,mog_circ,database)
    end
    return mog_circ
end

#takes the results of each owpp to oss line and combines with the oss to pcc line to create a final circuit only for hvac/mvac connections to an HVAC MOG
function combine_connections_ACDCowpp2ossDC(owpp2oss_connections,mog_circ,database)
    #finds common equipment
    mog_circDC=combine_connections_DCowpp2ossDC(owpp2oss_connections["HVDC"],deepcopy(mog_circ),database)
    mog_circBest=combine_connections_Bestowpp2ossDC(owpp2oss_connections["Best"],deepcopy(mog_circ),database)
    if (mog_circDC.cost<mog_circBest)
        return mog_circDC
    else
        return mog_circBest
    end
end

function optimal_ACDCowpp2ossDC(km,mva,mog,database,ks)
    owpp2oss_connections=Dict{String,Dict{String,Any}}()
    owpp2oss_connectionsACDC=optimal_ACowpp2ossDC(km,mva,mog,database,ks)
    owpp2oss_connectionsDCDC=optimal_DCowpp2ossDC(km,mva,mog,database,ks)
    if (owpp2oss_connectionsACDC["cost"]<owpp2oss_connectionsDCDC["cost"])
        owpp2oss_connectionsBest=owpp2oss_connectionsACDC
    else
        owpp2oss_connectionsBest=owpp2oss_connectionsDCDC
    end
    push!(owpp2oss_connections,("HVDC"=>owpp2oss_connectionsDCDC))
    push!(owpp2oss_connections,("Best"=>owpp2oss_connectionsBest))
    return owpp2oss_connections
end

#takes the results of each owpp to oss line and combines with the oss to pcc line to create a final circuit only for hvac/mvac connections to an HVAC MOG
function combine_connections_ACowpp2ossDC(owpp2oss_connections,mog_circ,database)
    #finds common equipment
    mog_circ=combine_connections_ACowpp2ossACDC(owpp2oss_connections,mog_circ,database)
    #AC base platform
    plat_base=platform()
    plat_base.acdc="ac"
    plat_base.mva=mog_circ.PCCcbls[1].mva
    plat_base.wnd=mog_circ.PCCcbls[1].wnd
    plat_base=cost_ac_platform(plat_base,ks)
    push!(mog_circ.mog[1].plat,plat_base)
    #finds DC platform
    plat_dc=platform()
    plat_dc.acdc="dc"
    plat_dc.mva=mog_circ.PCCcbls[1].mva
    plat_dc.wnd=mog_circ.PCCcbls[1].wnd
    plat_dc=cost_dc_platform(plat_dc,ks)
    plat_dc=adjust_base_dc_platform(plat_dc,ks)
    push!(mog_circ.mog[1].plat,plat_dc)
    #Finds DC converter
    conv=deepcopy(database["converters"][string(mog_circ.PCCcbls[1].mva)]["offshore"])
    conv.wnd=mog_circ.PCCcbls[1].wnd
    conv=cost_hvdc_oss(conv,ks)
    conv=adjust_base_hvdc_offshore_converter(conv,ks)
    push!(mog_circ.mog[1].conv,conv)
    #Totals cost
    mog_circ=total_circuit_cost(mog_circ)
    return mog_circ
end

#takes the results of each owpp to oss line and combines with the oss to pcc line to create a final circuit only for hvac/mvac connections to an HVAC MOG
function combine_connections_ACowpp2ossAC(owpp2oss_connections,mog_circ,database)
    #finds common equipment
    mog_circ=combine_connections_ACowpp2ossACDC(owpp2oss_connections,mog_circ,database)
    #finds AC platform
    plat=platform()
    plat.mva=mog_circ.PCCcbls[1].mva
    plat.wnd=mog_circ.PCCcbls[1].wnd
    plat=cost_ac_platform(plat,ks)
    plat=adjust_base_ac_platform(plat,ks)
    push!(mog_circ.mog[1].plat,plat)
    mog_circ=total_circuit_cost(mog_circ)
    return mog_circ
end

function connections_DCowpp2ossDC(owpp2oss_connection,mog_circ)
    #build the offshore MOG
    push!(mog_circ.oss,aft_connection_DC_equipment(owpp2pcc_connection))
    #lay the cable
    push!(mog_circ.HVcbls,owpp2oss_connection["cable"])
    return
end

function connections_ACowpp2ossACDC(owpp2oss_connection,mog_circ)
    #create an OSS bus
    oss=bus()
    oss.kV=owpp2oss_connection["cable"].elec.volt
    oss.mva=owpp2oss_connection["cable"].mva
    oss.wnd=owpp2oss_connection["cable"].wnd
    oss.node=owpp2oss_connection["cable"].path[1]
    #build the platform
    push!(oss.plat,owpp2oss_connection["plat_aft_ac"])
    #place the transformers
    push!(oss.xfmrs,owpp2oss_connection["xfo_aft"])
    push!(mog_circ.oss,oss)
    #lay the cable
    push!(mog_circ.HVcbls,owpp2oss_connection["cable"])
    return mog_circ
end

function combine_connections_DCowpp2ossDC(owpp2oss_connections,mog_circ,database)
    for owpp2oss_connection in owpp2oss_connections
        #create an OSS bus
        mog_circ=connections_DCowpp2ossDC(owpp2oss_connection,mog_circ)
    end
    #finds AC as there is no converter at MOG it is just a collection platform for DC
    plat=platform()
    plat.mva=mog_circ.PCCcbls[1].mva
    plat.wnd=mog_circ.PCCcbls[1].wnd
    plat=cost_ac_platform(plat,ks)
    plat=adjust_base_ac_platform(plat,ks)
    push!(mog_circ.mog[1].plat,plat)
    mog_circ=total_circuit_cost(mog_circ)
    return mog_circ
end

function combine_connections_Bestowpp2ossDC(owpp2oss_connections,mog_circ,database)
    mva66=0.0
    mva220=0.0
    mva400=0.0
    wnd66=Array{wnd,1}()
    wnd220=Array{wnd,1}()
    wnd400=Array{wnd,1}()
    for owpp2oss_connection in owpp2oss_connections
        if (haskey(owpp2oss_connection, "xfo_aft"))#hvac connection
            mog_circ=connections_ACowpp2ossACDC(owpp2oss_connection,mog_circ)
            #record powers and winds per voltage level to build forward platform transformers
            if (owpp2oss_connection["cable"].elec.volt.kV==220.0)
                mva220=mva220+owpp2oss_connection["cable"].mva
                push!(wnd220,owpp2oss_connection["cable"].wnd)
            else
                mva400=mva400+owpp2oss_connection["cable"].mva
                push!(wnd400,owpp2oss_connection["cable"].wnd)
            end
        elseif (haskey(owpp2oss_connection, "conv_aft"))#hvdc connection
            mog_circ=connections_DCowpp2ossDC(owpp2oss_connection,mog_circ)
        else#mvac
            #lay the cable
            push!(mog_circ.MVcbls,owpp2oss_connection["cable"])
            #track powers and winds per 66kv to build forward platform transformers
            mva66=mva66+owpp2oss_connection["cable"].mva
            push!(wnd66,owpp2oss_connection["cable"].wnd)
        end
    end
    mog_circ.mog[1].xfmrs=Array{transformer,1}()
    if (mva66>0.0)
        xfo66=deepcopy(database["transformers"][string(mva66)]["offshore"])
        xfo66.wnd=find_netWind(wnd66)
        xfo66.elec.lv=mv_cable.elec.volt
        xfo66.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo66=cost_xfo_oss(xfo66,ks)
        push!(mog_circ.mog[1].xfmrs,xfo66)
    end
    if (mva220>0.0 && (mog_circ[1].kV != 220.0))
        xfo220=deepcopy(database["transformers"][string(mva220)]["offshore"])
        xfo220.wnd=find_netWind(wnd220)
        xfo220.elec.lv=mv_cable.elec.volt
        xfo220.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo220=cost_xfo_oss(xfo220,ks)
        push!(mog_circ.mog[1].xfmrs,xfo220)
    end
    if (mva400>0.0 && (mog_circ[1].kV != 400.0))
        xfo400=deepcopy(database["transformers"][string(mva400)]["offshore"])
        xfo400.wnd=find_netWind(wnd400)
        xfo400.elec.lv=mv_cable.elec.volt
        xfo400.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo400=cost_xfo_oss(xfo400,ks)
        push!(mog_circ.mog[1].xfmrs,xfo400)
    end
    conv_mva=mva66+mva220+mva400
    if (conv_mva>0.0)
        conv_wnd.wnd=find_netWind([find_netWind(wnd66),find_netWind(wnd220),find_netWind(wnd400)])
        #platform
        plat_base=platform()
        plat_base.acdc="ac"
        plat_base.mva=mog_circ.PCCcbls[1].mva
        plat_base.wnd=mog_circ.PCCcbls[1].wnd
        plat_base=cost_ac_platform(plat_base,ks)
        push!(mog_circ.mog[1].plat,plat_base)
        #finds DC platform
        plat=platform()
        plat.acdc="dc"
        plat.mva=conv_mva
        plat.wnd=conv_wnd
        plat=cost_dc_platform(plat,ks)
        plat=adjust_base_dc_platform(plat,ks)
        push!(mog_circ.mog[1].plat,plat)
        #Finds DC converter
        conv=deepcopy(database["converters"][string(conv_mva)]["offshore"])
        conv.wnd=conv_wnd
        conv=cost_hvdc_oss(conv,ks)
        conv=adjust_base_hvdc_offshore_converter(conv,ks)
        push!(mog_circ.mog[1].conv,conv)
        mog_circ=total_circuit_cost(mog_circ)
    else
        mog_circ.cost=Inf
    end
    return mog_circ
end

#finds common equipment of circuit for ACowpp2ossAC and ACowpp2ossDC
function combine_connections_ACowpp2ossACDC(owpp2oss_connections,mog_circ,database)
    mva66=0.0
    mva220=0.0
    mva400=0.0
    wnd66=Array{wnd,1}()
    wnd220=Array{wnd,1}()
    wnd400=Array{wnd,1}()
    for owpp2oss_connection in owpp2oss_connections
        if (owpp2oss_connection["HVAC"]["cost"]<owpp2oss_connection["MVAC"]["cost"])
            mog_circ=connections_ACowpp2ossACDC(owpp2oss_connection["HVAC"],mog_circ)
            #record powers and winds per voltage level to build forward platform transformers
            if (owpp2oss_connection["HVAC"]["cable"].elec.volt.kV==220.0)
                mva220=mva220+owpp2oss_connection["HVAC"]["cable"].mva
                push!(wnd220,owpp2oss_connection["HVAC"]["cable"].wnd)
            else
                mva400=mva400+owpp2oss_connection["HVAC"]["cable"].mva
                push!(wnd400,owpp2oss_connection["HVAC"]["cable"].wnd)
            end
        else#MV connection
            #lay the cable
            push!(mog_circ.MVcbls,owpp2oss_connection["MVAC"]["cable"])
            #track powers and winds per 66kv to build forward platform transformers
            mva66=mva66+owpp2oss_connection["MVAC"]["cable"].mva
            push!(wnd66,owpp2oss_connection["MVAC"]["cable"].wnd)
        end
    end
    mog_circ.mog[1].xfmrs=Array{transformer,1}()
    if (mva66>0.0)
        xfo66=deepcopy(database["transformers"][string(mva66)]["offshore"])
        xfo66.wnd=find_netWind(wnd66)
        xfo66.elec.lv=mv_cable.elec.volt
        xfo66.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo66=cost_xfo_oss(xfo66,ks)
        push!(mog_circ.mog[1].xfmrs,xfo66)
    end
    if (mva220>0.0 && (mog_circ[1].kV != 220.0))
        xfo220=deepcopy(database["transformers"][string(mva220)]["offshore"])
        xfo220.wnd=find_netWind(wnd220)
        xfo220.elec.lv=mv_cable.elec.volt
        xfo220.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo220=cost_xfo_oss(xfo220,ks)
        push!(mog_circ.mog[1].xfmrs,xfo220)
    end
    if (mva400>0.0 && (mog_circ[1].kV != 400.0))
        xfo400=deepcopy(database["transformers"][string(mva400)]["offshore"])
        xfo400.wnd=find_netWind(wnd400)
        xfo400.elec.lv=mv_cable.elec.volt
        xfo400.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo400=cost_xfo_oss(xfo400,ks)
        push!(mog_circ.mog[1].xfmrs,xfo400)
    end
    return mog_circ
end


#finds best alternative from 66-platform, 66-OSS-HVAC-platform, 66-OSS-HVDC-platform
#a cost function
function connection_ACowpp2ossACDC(km,owpp,mog,database,ks)
    #option 1 66-platform
    owpp2oss_connectionAC=Dict{String,Dict{String,Any}}()
    owpp2oss_connectionHVAC=Dict{String,Any}()
    owpp2oss_connectionMVAC=Dict{String,Any}()
    push!(owpp2oss_connectionHVAC,("cost"=>Inf))#cost
    push!(owpp2oss_connectionMVAC,("cost"=>Inf))#cost
    if (owpp.mv_zone>km)
        owpp2oss_connectionMVAC=connection_owpp2oss_mv2AC_feeder(owpp,km,database,ks)
        push!(owpp2oss_connectionMVAC["cable"].path,owpp.node)
        push!(owpp2oss_connectionMVAC["cable"].path,mog.node)
    else
        #option 2 66-OSS-HVAC-platform
        if (database["bits"]["hvac"]==true)
            owpp2oss_connectionHVAC=connection_owpp2oss_hv2AC_feeder(owpp,km,mog.kv,database,ks)
            push!(owpp2oss_connectionHVAC["cable"].path,owpp.node)
            push!(owpp2oss_connectionHVAC["cable"].path,mog.node)
        else
            println("Error! connection_ACowpp2ossACDC: hvac not selected and mv_zone out of range.")
        end
    end
    push!(owpp2oss_connectionAC,("HVAC"=>owpp2oss_connectionHVAC))
    push!(owpp2oss_connectionAC,("MVAC"=>owpp2oss_connectionMVAC))
    return owpp2oss_connectionAC
end

#cost function
function connection_owpp2oss_mv2AC_feeder(owpp,km,database,ks)
    #cable
    mv_cable=deepcopy(mvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))

    #totals and organizes mv equipment
    owpp2oss_connection66=Dict{String,Any}()
    push!(owpp2oss_connection66,("cable"=>mv_cable))#cable
    push!(owpp2oss_connection66,("cost"=>mv_cable.costs.ttl))#cost
    return owpp2oss_connection66
end

#cost function
function connection_owpp2oss_hv2AC_feeder(owpp,km,kv,database,ks)
    owpp2oss_connectionAC=Dict{String,Any}()
    #build candidate cables
    cable_220=deepcopy(hvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["220.0"][string(owpp.mva)],ks))
    cable_400=deepcopy(hvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["400.0"][string(owpp.mva)],ks))
    cable_220.costs.grand_ttl=cable_220.costs.ttl
    cable_400.costs.grand_ttl=cable_400.costs.ttl
    #if mid point compensation is checked look for cheaper cable options
    if (database["bits"]["mpc_ac"]==true)
        cable_220,cable_400=check4beter_mpc_ac_cables(owpp.mva,km,database,owpp.wnd,cable_220,cable_400,ks)
    else
    end
    #aft candidate transformers
    xfo_aft=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
    xfo_aft.wnd=owpp.wnd
    xfo_aft=cost_xfo_oss(xfo_aft,ks)

    #platform
    plat_aft=platform()
    plat_aft.acdc="ac"
    plat_aft.mva=owpp.mva
    plat_aft.wnd=owpp.wnd
    plat_aft=cost_ac_platform(plat_aft,ks)
    plat_aft=adjust_base_ac_platform(plat_aft,ks)

    #forward transformers
    if (kv != 220.0)
        xfo_for220=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
        xfo_for220.wnd=owpp.wnd
        xfo_for220=cost_xfo_oss(xfo_aft,ks)
    else
        xfo_for220=transformer()
        xfo_for220.costs.ttl=0.0
    end
    #forward transformers
    if (kv != 400.0)
        xfo_for400=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
        xfo_for400.wnd=owpp.wnd
        xfo_for400=cost_xfo_oss(xfo_aft,ks)
    else
        xfo_for400=transformer()
        xfo_for400.costs.ttl=0.0
    end

    #totals
    owpp2oss_connection220_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_220.costs.grand_ttl+xfo_for220.costs.ttl
    owpp2oss_connection400_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_400.costs.grand_ttl+xfo_for400.costs.ttl
    owpp2oss_connectionAC_cost=findmin([owpp2oss_connection220_cost,owpp2oss_connection400_cost])
    owpp2oss_connection220=[plat_aft,xfo_aft,cable_220,xfo_for220]
    owpp2oss_connection400=[plat_aft,xfo_aft,cable_400,xfo_for400]
    owpp2oss_connectionAC=[owpp2oss_connection220,owpp2oss_connection400][owpp2oss_connectionAC_cost[2]]

    #totals and organize hv equipment
    push!(owpp2oss_connectionAC,("plat_aft_ac"=>owpp2oss_connectionAC[1]))#aft platform
    push!(owpp2oss_connectionAC,("xfo_aft"=>owpp2oss_connectionAC[2]))#aft platform
    push!(owpp2oss_connectionAC,("cable"=>owpp2oss_connectionAC[3]))#cable
    push!(owpp2oss_connectionAC,("xfo_for"=>owpp2oss_connectionAC[4]))#xfo_for
    push!(owpp2oss_connectionAC,("cost"=>owpp2oss_connectionAC_cost[1]))#cost
    return owpp2oss_connectionAC
end

#finds best alternative from 66-platform, 66-OSS-HVAC-platform, 66-OSS-HVDC-platform
function optimal_ACowpp2ossDC(km,owpp,mog,database,ks)
    #option 1 66-platform
    owpp2oss_connectionAC=Dict{String,Any}()
    push!(owpp2oss_connectionAC,("cost"=>Inf))#cost
    if (owpp.mv_zone>km)
        owpp2oss_connectionAC=connection_owpp2oss_mv2DC_feeder(owpp,km,database,ks)
    else
        #option 2 66-OSS-HVAC-platform
        if (database["bits"]["hvac"]==true)
            owpp2oss_connectionAC=connection_owpp2oss_hv2DC_feeder(owpp,km,database,ks)
        end
    end
    push!(owpp2oss_connectionAC["cable"].path,owpp.node)
    push!(owpp2oss_connectionAC["cable"].path,mog.node)
    return owpp2oss_connectionAC
end


#finds best alternative from 66-platform, 66-OSS-HVAC-platform, 66-OSS-HVDC-platform
function optimal_DCowpp2ossDC(km,owpp,mog,database,ks)

    owpp2oss_connectionDC=Dict{String,Any}()
    push!(owpp2oss_connectionDC,("cost"=>Inf))#cost
    if (database["bits"]["hvdc"]==true)
        owpp2oss_connectionDC=connection_owpp2oss_dc2DC_feeder(owpp,km,database,ks)
    end
    #option 3 66-OSS-HVDC-platform
    push!(owpp2oss_connectionDC["cable"].path,owpp.node)
    push!(owpp2oss_connectionDC["cable"].path,mog.node)
    return owpp2oss_connectionDC
end

#option 3 of optimal_owpps2oss_DC_feeder() 66-OSS-HVDC-platform
function connection_owpp2oss_dc2DC_feeder(owpp,km,database,ks)
    owpp2oss_connectionHVDC=Dict{String,Any}()
    cable_dc=deepcopy(hvdc_cable(owpp.mva,km,owpp.wnd,database["cables"]["300.0"][string(owpp.mva)],ks))
    #platform
    plat_base=platform()
    plat_base.acdc="ac"
    plat_base.mva=owpp.mva
    plat_base.wnd=owpp.wnd
    plat_base=cost_ac_platform(plat_base,ks)
    #hvdc platform
    plat_dc=platform()
    plat_dc.acdc="dc"
    plat_dc.mva=owpp.mva
    plat_dc.wnd=owpp.wnd
    plat_dc=cost_dc_platform(plat_dc,ks)
    plat_dc=adjust_base_dc_platform(plat_dc,ks)
    #aft candidate transformers
    xfo=transformer()
    xfo=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
    xfo.wnd=owpp.wnd
    xfo.elec.hv=300.0
    xfo.elec.lv=66.0
    xfo=cost_xfo_oss(xfo,ks)
    #converter
    conv=deepcopy(database["converters"][string(owpp.mva)]["offshore"])
    conv.wnd=owpp.wnd
    conv=cost_hvdc_oss(conv,ks)

    #totals and organizes mv equipment
    owpp2oss_connectionHVDC_cost=plat_base.costs.ttl+plat_dc.costs.ttl+xfo.costs.ttl+conv.costs.ttl+cable_dc.ttl
    push!(owpp2oss_connectionHVDC,("plat_aft_ac"=>plat_base))#aft platform
    push!(owpp2oss_connectionHVDC,("plat_aft_dc"=>plat_base))#aft platform
    push!(owpp2oss_connectionHVDC,("conv_aft"=>conv))#xfo_for
    push!(owpp2oss_connectionHVDC,("xfo_aft"=>xfo))#xfo_for
    push!(owpp2oss_connectionHVDC,("cable"=>cable_dc))#cable
    push!(owpp2oss_connectionHVDC,("cost"=>owpp2oss_connectionHVDC_cost))#cost
    return owpp2oss_connectionHVDC
end
#option 2 of optimal_owpps2oss_DC_feeder() 66-OSS-HVAC-platform
function connection_owpp2oss_hv2DC_feeder(owpp,km,database,ks)
    owpp2oss_connectionAC=Dict{String,Any}()
    #build candidate cables
    cable_220=deepcopy(hvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["220.0"][string(owpp.mva)],ks))
    cable_400=deepcopy(hvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["400.0"][string(owpp.mva)],ks))
    cable_220.costs.grand_ttl=cable_220.costs.ttl
    cable_400.costs.grand_ttl=cable_400.costs.ttl
    #if mid point compensation is checked look for cheaper cable options
    if (database["bits"]["mpc_ac"]==true)
        cable_220,cable_400=check4beter_mpc_ac_cables(owpp.mva,km,database,owpp.wnd,cable_220,cable_400,ks)
    else
    end
    hv_cable_cost=argmin([cable_220.costs.grand_ttl,cable_400.costs.grand_ttl])
    hv_cable=[[cable_220,cable_400][hv_cable_cost]]
    #aft candidate transformers
    xfo_aft=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
    xfo_aft.wnd=owpp.wnd
    xfo_aft=cost_xfo_oss(xfo_aft,ks)

    #platform
    plat_aft_ac=platform()
    plat_aft_ac.acdc="ac"
    plat_aft_ac.mva=owpp.mva
    plat_aft_ac.wnd=owpp.wnd
    plat_aft_ac=cost_ac_platform(plat_aft_ac,ks)
    plat_aft_ac=adjust_base_ac_platform(plat_aft_ac,ks)

    #hvdc platform
    plat_for_dc=platform()
    plat_for_dc.acdc="dc"
    plat_for_dc.mva=owpp.mva
    plat_for_dc.wnd=owpp.wnd
    plat_for_dc=cost_dc_platform(plat_for_dc,ks)

    #transformer
    xfo_for=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
    xfo_for.wnd=owpp.wnd
    xfo_for.elec.hv=300.0
    xfo_for.elec.lv=hv_cable.elec.volt
    xfo_for=cost_xfo_oss(xfo_for,ks)

    #converter
    conv=deepcopy(database["converters"][string(owpp.mva)]["offshore"])
    conv.wnd=owpp.wnd
    conv=cost_hvdc_oss(conv,ks)

    #totals and organizes mv equipment
    owpp2oss_connectionAC_cost=plat_aft_ac.costs.ttl+xfo_aft.costs.ttl+hv_cable.costs.grand_ttl+xfo_for.costs.ttl+conv.costs.grand_ttl+plat_for_dc.costs.ttl
    push!(owpp2oss_connectionAC,("plat_aft_ac"=>plat_aft_ac))#aft platform
    push!(owpp2oss_connectionAC,("xfo_aft"=>xfo_aft))#aft platform
    push!(owpp2oss_connectionAC,("cable"=>hv_cable))#cable
    push!(owpp2oss_connectionAC,("plat_for_dc"=>plat_for_dc))#xfo_for
    push!(owpp2oss_connectionAC,("xfo_for"=>xfo_for))#xfo_for
    push!(owpp2oss_connectionAC,("conv_for"=>conv))#xfo_for
    push!(owpp2oss_connectionAC,("cost"=>owpp2oss_connectionAC_cost))#cost
    return owpp2oss_connectionAC
end

function connection_owpp2oss_mv2DC_feeder(owpp,km,database,ks)
    owpp2oss_connection66=Dict{String,Any}()
    #cable
    mv_cable=deepcopy(mvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
    mv_cable.costs.grand_ttl=mv_cable.costs.ttl

    #hvdc platform
    plat_for_dc=platform()
    plat_for_dc.acdc="dc"
    plat_for_dc.mva=owpp.mva
    plat_for_dc=cost_ac_platform(plat_for_dc,ks)

    #converter
    mv_conv=deepcopy(database["converters"][string(owpp.mva)]["offshore"])
    mv_conv.wnd=owpp.wnd
    mv_conv=cost_hvdc_oss(mv_conv,ks)

    #transformer
    xfo=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
    xfo.wnd=owpp.wnd
    xfo.elec.hv=300.0
    xfo.elec.lv=mv_cable.elec.volt
    xfo=cost_xfo_oss(xfo,ks)

    #totals and organizes mv equipment
    owpp2oss_connection66_cost=mv_cable.costs.ttl+plat_for_dc.costs.ttl+mv_conv.costs.ttl+xfo.costs.ttl
    push!(owpp2oss_connection66,("cable"=>mv_cable))#cable
    push!(owpp2oss_connection66,("plat_for_dc"=>plat_for_dc))#xfo_for
    push!(owpp2oss_connection66,("xfo_for"=>xfo_for))#xfo_for
    push!(owpp2oss_connection66,("conv_for"=>mv_conv))#xfo_for
    push!(owpp2oss_connection66,("cost"=>owpp2oss_connection66_cost))#cost
    return owpp2oss_connection66
end

#finds the optimal connection for the main MOG to PCC connection
function optimal_oss2oss(km,mva,kv_aft,mvas_aft,mvas_forward,kv_forward,wnd,database)#kvs_aft=[mva(66.0),mva(220.0),mva(400.0),mva(300.0)]
    ks=get_Cost_Data()
    #HVAC connection
    #will return AC-AC/max and min of AC-DC
    if (kv_aft != 300.0)
        connection_from_AC=hvac_canditate_oss2oss(km,mva,kv_aft,mvas_aft,mvas_forward,kv_forward,wnd,database,ks)#returns dictionary: [platform,trans_aft,cable,trans_for]
    else
        #set connection_from_AC to infinite
    end
    dc_connection,cost_dc_connection=hvdc_canditate_oss2oss(km,mva,kv_aft,mvas_aft,mvas_forward,kv_forward,wnd,database,ks)#returns dictionary: [platform,conv_aft,cable,conv_for]
    #mid point compensated HVAC connection - can be handled inside HVAC - should update
    #mpc_connection,cost_mpc_connection=mpc_canditate_oss2oss(km,mva,kvs_aft,kv_forward,wnd,database,ks)#returns dictionary: [platform,trans_aft,cable,xfo_for]

    #connection costs
    lowest_cost=findmin([cost_dc_connection,cost_ac_connection,cost_mpc_connection])
    optimal_equipment=[dc_connection,ac_connection,mpc_connection][lowest_cost[2]]
    return optimal_equipment,lowest_cost[1]
end

#HVAC connection from an OSS to an OSS
function hvac_canditate_oss2oss(km,mva,kv_aft,mvas_aft,mvas_forward,kv_forward,wnd,database,ks)
    #candidate cables
    cable_220=cable()
    cable_220.costs.ttl=Inf
    cable_400=cable()
    cable_400.costs.ttl=Inf
    #candidate aft transformers
    xfo_aft66=transformer()
    xfo_aft66.costs.ttl=0
    xfo_aft400220=transformer()
    xfo_aft400220.costs.ttl=0
    xfo_aft220400=transformer()
    xfo_aft220400.costs.ttl=0
    #candidate forward transformers
    xfo_for220=transformer()
    xfo_for220.costs.ttl=0
    xfo_for400=transformer()
    xfo_for400.costs.ttl=0
    #candidate forward transformers
    conv_for=converter()
    conv_for.costs.grand_ttl=0
    #candidate aft AC platform
    plat_aft=platform()
    plat_aft.costs.ttl=Inf
    #candidate forward upgrade if connecting to DC platform
    plat_for=platform()
    plat_for.costs.ttl=Inf
    #If AC is not included Inf price ensures it is not selected
    cost_acac_connection=Inf
    cost_acdc_connection=Inf
    acac_connection=[plat_aft,xfo_aft,cable_220,xfo_for220]
    acdc_connection=[plat_aft,xfo_aft,cable_220,conv_for]
    if (database["bits"]["hvac"]==true)#eliminate dropping from HVDC to HVAC (as length(oss0-pcc)<length(oss0-oss1-pcc)) and mva is constant. It would not make sense to downgrade
        cable_220=deepcopy(hvac_cable(mva,km,wnd,database["cables"]["220.0"][string(mva)],ks))
        cable_400=deepcopy(hvac_cable(mva,km,wnd,database["cables"]["400.0"][string(mva)],ks))
        cable_220.costs.grand_ttl=cable_220.costs.ttl
        cable_400.costs.grand_ttl=cable_400.costs.ttl
        if (database["bits"]["mpc_ac"]==true)
            cable_220mpc=deepcopy(hvac_cable(mva,km/2,wnd,database["cables"]["220.0"][string(mva)],ks))
            cable_400mpc=deepcopy(hvac_cable(mva,km/2,wnd,database["cables"]["400.0"][string(mva)],ks))
            cable_220mpc.mpc_ac=true
            cable_400mpc.mpc_ac=true
            #platforms
            #compensation
            plat_mpc=platform()
            plat_mpc.mva=mva
            plat_mpc.wnd=wnd
            plat_mpc=cost_ac_platform(plat_mpc,ks)
            plat_mpc=adjust_base_ac_platform(plat_mpc,ks)
            cable_220mpc.plat=plat_mpc
            cable_400mpc.plat=plat_mpc
            cable_220mpc.costs.grand_ttl=cable_220mpc.costs.ttl*2+plat_mpc.costs.ttl
            cable_400mpc.costs.grand_ttl=cable_400mpc.costs.ttl*2+plat_mpc.costs.ttl
            if (cable_220mpc.costs.grand_ttl<cable_220.costs.grand_ttl)
                cable_220=cable_220mpc
            end
            if (cable_400mpc.costs.grand_ttl<cable_400.costs.grand_ttl)
                cable_400=cable_400mpc
            end
        end
        #aft
        #transformers
        xfo_aft66=deepcopy(database["transformers"][string(kvs_aft["66.0"])]["offshore"])
        xfo_aft400220=deepcopy(database["transformers"][string(kvs_aft["400.0"])]["offshore"])
        xfo_aft220400=deepcopy(database["transformers"][string(kvs_aft["220.0"])]["offshore"])

        #platform
        plat_aft.acdc="ac"
        plat_aft.mva=mva
        plat_aft.wnd=wnd
        plat_aft=cost_ac_platform(plat_aft,ks)
        plat_aft=adjust_base_ac_platform(plat_aft,ks)

        #forward transformera
        if (kv!=400.0)
            xfo_for400=deepcopy(database["transformers"][string(kvs_forward["400.0"]+mva)]["offshore"])
            xfo_for400.costs.ttl=xfo_for400.costs.ttl*mva/(kvs_forward["400.0"]+mva)
        else
            xfo_for400.costs.ttl=0
        end
        #forward
        if (kv_forward!=220.0)
            xfo_for400=deepcopy(database["transformers"][string(kvs_forward["220.0"]+mva)]["offshore"])
            xfo_for400.costs.ttl=xfo_for400.costs.ttl*mva/(kvs_forward["220.0"]+mva)
        else
            xfo_for220.costs.ttl=0
        end
        if (database["bits"]["hvdc"]==true)
            #candidate forward converters
            conv_for400=deepcopy(database["converters"][string(mva)]["offshore"])
            conv_for220=deepcopy(database["converters"][string(mva)]["offshore"])
            conv_for400.xfo=xfo_for400
            conv_for220.xfo=xfo_for220
            #candidate aft DC platform
            conv_for400.plat.mva=mva
            conv_for400.plat.wnd=wnd
            conv_for400.plat=cost_dc_platform(conv_for400.plat,ks)
            conv_for400=hvdc_offshore_station_cost(conv_for400)
            conv_for220.plat.mva=mva
            conv_for220.plat.wnd=wnd
            conv_for220.plat=cost_dc_platform(conv_for220.plat,ks)
            conv_for220=hvdc_offshore_station_cost(conv_for220)
            #NOTE all the combinations 220-300, 400-300 find min of to dc circuits here ie eliminate an AC voltage for a dc connection
            #This function needs extreme testing, figure out a way to make certain!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            cost_acdc_connection=findmin([])
            acac_connection=[plat_aft,xfo_aft,cable_220,xfo_for220]
            acdc_connection=[plat_aft,xfo_aft,cable_220,conv_for]
        end

        #find best alternative PCC connection voltage
        ac_connection220_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_220.costs.ttl+xfo_for220.costs.ttl
        ac_connection400_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_400.costs.ttl+xfo_for400.costs.ttl
        best_voltage=findmin([ac_connection220_cost,ac_connection400_cost])
        ac_connection220=[plat_aft,xfo_aft,cable_220,xfo_for220]
        ac_connection400=[plat_aft,xfo_aft,cable_400,xfo_for400]
        ac_connection=[ac_connection220,ac_connection400][best_voltage[2]]
        cost_ac_connection=best_voltage[1]
    end
    dict_ac_connection=Dict{String,Any}()
    push!(dict_ac_connection,("plat_aft"=>ac_connection[1]))#plat_aft
    push!(dict_ac_connection,("xfo_aft"=>ac_connection[2]))#xfo_aft
    push!(dict_ac_connection,("cable"=>ac_connection[3]))#cable
    push!(dict_ac_connection,("xfo_for"=>ac_connection[4]))#xfo_for
    return dict_ac_connection,cost_ac_connection
end

#hvdc OSS to OSS candidate connection
function hvdc_canditate_oss2oss(km,mva,kv_aft,kv_forward,wnd,database,ks)
    #DC cables
    cable_dc=cable()
    cable_dc.costs.ttl=Inf
    #aft converter
    conv_aft=converter()
    conv_aft.costs.ttl=Inf
    #forward converter
    conv_for=converter()
    conv_for.costs.ttl=Inf
    #aft transformer
    xfo_aft=transformer()
    xfo_aft.costs.ttl=Inf
    #forward transformer
    xfo_for=transformer()
    xfo_for.costs.ttl=Inf
    #aft platform
    plat_aft=platform()
    plat_aft.costs.ttl=Inf
    #If no HVDC being considered infinite cost ensures it is not selected
    cost_dc_connection=Inf
    dc_connection=[plat_aft,xfo_aft,conv_aft,cable_dc,conv_for,xfo_for]

    if (database["bits"]["hvdc"]==true)
        cable_dc=deepcopy(hvdc_cable(mva,km,wnd,database["cables"]["300.0"][string(mva)],ks))
        #aft
        #converters
        conv_aft=deepcopy(database["converters"][string(mva)]["offshore"])
        conv_aft.wnd=wnd
        conv_aft=cost_hvdc_oss(conv_aft,ks)
        #transformers
        conv_aft.xfo=deepcopy(database["transformers"][string(mva)]["offshore"])
        conv_aft.xfo.wnd=wnd
        conv_aft.xfo=cost_xfo_oss(conv_aft.xfo,ks)
        #hvdc platform
        conv_aft.plat.mva=mva
        conv_aft.plat.wnd=wnd
        conv_aft.plat=cost_dc_platform(conv_aft.plat,ks)
        conv_aft=hvdc_offshore_station_cost(conv_aft)

        plat_aft.acdc="dc"
        plat_aft.mva=mva
        plat_aft=cost_ac_platform(plat_aft,ks)
        #forward
        #converter
        conv_for=deepcopy(database["converters"][string(mva)]["onshore"])
        conv_for.wnd=wnd
        conv_for=cost_hvdc_pcc(conv_for,ks)
        #transformers
        conv_for.xfo=deepcopy(database["transformers"][string(mva)]["onshore"])
        conv_for.xfo.wnd=wnd
        conv_for.xfo=cost_xfo_pcc(conv_for.xfo,ks)

        cost_dc_connection=plat_aft.costs.ttl+conv_aft.costs.grand_ttl+cable_dc.costs.ttl+conv_for.costs.grand_ttl
        dc_connection=[plat_aft,conv_aft,cable_dc,conv_for]
    end
    dict_dc_connection=Dict{String,Any}()
    push!(dict_dc_connection,("plat_aft"=>dc_connection[1]))#plat_aft
    push!(dict_dc_connection,("conv_aft"=>dc_connection[3]))#conv_aft
    push!(dict_dc_connection,("cable"=>dc_connection[4]))#cable
    push!(dict_dc_connection,("conv_for"=>dc_connection[5]))#conv_for
    return dict_dc_connection,cost_dc_connection
end
