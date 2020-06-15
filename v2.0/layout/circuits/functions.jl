####################### HVDC connection cost OWPP to OSS #############################
#cost function
function owpp_dc_to_oss(km,owpp,oss_location,database,ks,hv_only)
    #option 1 66-platform
    owpp2oss_connectionHVDC=Dict{String,Any}()
    push!(owpp2oss_connectionHVDC,("cost"=>Inf))#cost
    if (database["bits"]["hvdc"]==true)
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
        conv=adjust_base_hvdc_offshore_converter(conv,ks)
        #totals and organizes mv equipment
        owpp2oss_connectionHVDC_cost=plat_base.costs.ttl+plat_dc.costs.ttl+xfo.costs.ttl+conv.costs.ttl+cable_dc.costs.ttl
        push!(owpp2oss_connectionHVDC,("plat_aft_ac"=>plat_base))#aft platform
        push!(owpp2oss_connectionHVDC,("plat_aft_dc"=>plat_dc))#aft platform
        push!(owpp2oss_connectionHVDC,("conv_aft"=>conv))#xfo_for
        push!(owpp2oss_connectionHVDC,("xfo_aft"=>xfo))#xfo_for
        push!(owpp2oss_connectionHVDC,("cable"=>cable_dc))#cable
        push!(owpp2oss_connectionHVDC,("cost"=>owpp2oss_connectionHVDC_cost))#cost
    end
    return owpp2oss_connectionHVDC
end

function oss_dc_to_oss(km,oss,wnd_power,database,ks)
    #option 1 66-platform
    owpp2oss_connectionHVDC=Dict{String,Any}()
    push!(owpp2oss_connectionHVDC,("cost"=>Inf))#cost
    if (database["bits"]["hvdc"]==true)
        cable_dc=deepcopy(hvdc_cable(oss.mva,km,oss.wnd,database["cables"]["300.0"][string(oss.mva)],ks))
        #platform
        plat_base=platform()
        plat_base.acdc="ac"
        plat_base.mva=oss.mva
        plat_base.wnd=oss.wnd
        plat_base=cost_ac_platform(plat_base,ks)
        #hvdc platform
        plat_dc=platform()
        plat_dc.acdc="dc"
        plat_dc.mva=oss.mva
        plat_dc.wnd=oss.wnd
        plat_dc=cost_dc_platform(plat_dc,ks)
        plat_dc=adjust_base_dc_platform(plat_dc,ks)
        #aft candidate transformers
        xfos=set_oss_transformers(300.0,wnd_power,database,ks)
        xfo_cost=0
        for xfo in xfos
            xfo_cost=xfo_cost+xfo.costs.ttl
        end
        #converter
        conv=set_oss_converters(300.0,wnd_power,database,ks)
        if (length(conv)==0)
            conv=[converter()]
            conv[1].costs.ttl=0.0
        end
        #totals and organizes mv equipment
        owpp2oss_connectionHVDC_cost=plat_base.costs.ttl+plat_dc.costs.ttl+xfo_cost+conv[1].costs.ttl+cable_dc.costs.ttl
        push!(owpp2oss_connectionHVDC,("plat_aft_ac"=>plat_base))#aft platform
        push!(owpp2oss_connectionHVDC,("plat_aft_dc"=>plat_dc))#aft platform
        push!(owpp2oss_connectionHVDC,("conv_aft"=>conv[1]))#xfo_for
        push!(owpp2oss_connectionHVDC,("xfo_aft"=>xfos))#xfo_for
        push!(owpp2oss_connectionHVDC,("cable"=>cable_dc))#cable
        push!(owpp2oss_connectionHVDC,("cost"=>owpp2oss_connectionHVDC_cost))#cost
    end
    return owpp2oss_connectionHVDC
end
########################## AC Connection OWPP to OSS ###########################
#finds best AC (HV or MV) connection from MVAC to HVAC OSS
function owpp_ac_to_oss(km,owpp,oss_location,database,ks,hv_only)
    #option 1 66-platform
    owpp2oss_connectionHVACs=Dict{String,Dict{String,Any}}()
    owpp2oss_connectionMVAC=Dict{String,Any}()
    owpp2oss_connectionHVAC=Dict{String,Any}()
    push!(owpp2oss_connectionMVAC,("cost"=>Inf))#cost
    push!(owpp2oss_connectionHVAC,("cost"=>Inf))#cost
    push!(owpp2oss_connectionHVACs,("220.0"=>owpp2oss_connectionHVAC))
    push!(owpp2oss_connectionHVACs,("400.0"=>owpp2oss_connectionHVAC))
    if (owpp.mv_zone>km && hv_only==false)
        owpp2oss_connectionMVAC=owpp_mvac_to_oss(owpp,km,database,ks)
        push!(owpp2oss_connectionMVAC["cable"].path,owpp.node)
        #push!(owpp2oss_connectionMVAC["cable"].path,oss_location)
    else
        #option 2 66-OSS-HVAC-platform
        if (database["bits"]["hvac"]==true)
            owpp2oss_connectionHVACs=owpp_hvac_to_oss(owpp,km,database,ks)
        else
            println("Error! connection_ACowpp2ossACDC: hvac not selected and mv_zone out of range.")
        end
    end
    push!(owpp2oss_connectionHVACs,("66.0"=>owpp2oss_connectionMVAC))
    return owpp2oss_connectionHVACs
end

####################### MVAC connection cost OWPP to OSS #############################
#cost function
function owpp_mvac_to_oss(owpp,km,database,ks)
    #cable
    mv_cable=deepcopy(mvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["66.0"][string(owpp.mva)],ks))
    #totals and organizes mv equipment
    owpp2oss_connection66=Dict{String,Any}()
    push!(owpp2oss_connection66,("cable"=>mv_cable))#cable
    push!(owpp2oss_connection66,("cost"=>mv_cable.costs.ttl))#cost
    return owpp2oss_connection66
end
####################### HVAC connection cost OWPP to OSS #############################
#cost function
function owpp_hvac_to_oss(owpp,km,database,ks)
    dict_ac_connection220=Dict{String,Any}()
    dict_ac_connection400=Dict{String,Any}()
    push!(dict_ac_connection220,("cost"=>Inf))#xfo_for
    push!(dict_ac_connection400,("cost"=>Inf))#xfo_for
    if (database["bits"]["hvac"]==true)
        #build candidate cables
        cable_220=deepcopy(hvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["220.0"][string(owpp.mva)],ks))
        cable_400=deepcopy(hvac_cable(owpp.mva,km,owpp.wnd,database["cables"]["400.0"][string(owpp.mva)],ks))
        cable_220.costs.grand_ttl=cable_220.costs.ttl
        cable_400.costs.grand_ttl=cable_400.costs.ttl
        #if mid point compensation is checked look for cheaper cable options
        cable_220,cable_400=check4beter_mpc_ac_cables(owpp.mva,km,database,owpp.wnd,cable_220,cable_400,ks)
        #aft candidate transformers
        xfo_aft=deepcopy(database["transformers"][string(owpp.mva)]["offshore"])
        xfo_aft.wnd=owpp.wnd
        xfo_aft=cost_xfo_oss(xfo_aft,ks)

        #platform
        plat_aft=platform()
        plat_aft.acdc="ac"
        plat_aft.mva=owpp.mva
        plat_aft.wnd=owpp.wnd
        plat_aft=cost_ac_platform(plat_aft,ks)#OK
        plat_aft=adjust_base_ac_platform(plat_aft,ks)

        #totals
        owpp2oss_connection220_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_220.costs.grand_ttl
        owpp2oss_connection400_cost=plat_aft.costs.ttl+xfo_aft.costs.ttl+cable_400.costs.grand_ttl

        #totals and organize hv equipment
        push!(dict_ac_connection220,("plat_aft_ac"=>plat_aft))#aft platform
        push!(dict_ac_connection220,("xfo_aft"=>xfo_aft))#aft platform
        push!(dict_ac_connection220,("cable"=>cable_220))#cable
        push!(dict_ac_connection220,("cost"=>owpp2oss_connection220_cost[1]))#cost

        push!(dict_ac_connection400,("plat_aft_ac"=>plat_aft))#aft platform
        push!(dict_ac_connection400,("xfo_aft"=>xfo_aft))#aft platform
        push!(dict_ac_connection400,("cable"=>cable_400))#cable
        push!(dict_ac_connection400,("cost"=>owpp2oss_connection400_cost[1]))#cost
    end
    dict_ac_connection=Dict{String,Dict{String,Any}}()
    push!(dict_ac_connection,("220.0"=>dict_ac_connection220))
    push!(dict_ac_connection,("400.0"=>dict_ac_connection400))
    return dict_ac_connection
end

function oss_hvac_to_oss(oss,wnd_power,km,database,ks)
    dict_ac_connection220=Dict{String,Any}()
    dict_ac_connection400=Dict{String,Any}()
    push!(dict_ac_connection220,("cost"=>Inf))#xfo_for
    push!(dict_ac_connection400,("cost"=>Inf))#xfo_for
    if (database["bits"]["hvac"]==true)
        #build candidate cables
        cable_220=deepcopy(hvac_cable(oss.mva,km,oss.wnd,database["cables"]["220.0"][string(oss.mva)],ks))
        cable_400=deepcopy(hvac_cable(oss.mva,km,oss.wnd,database["cables"]["400.0"][string(oss.mva)],ks))
        cable_220.costs.grand_ttl=cable_220.costs.ttl
        cable_400.costs.grand_ttl=cable_400.costs.ttl
        #if mid point compensation is checked look for cheaper cable options
        cable_220,cable_400=check4beter_mpc_ac_cables(oss.mva,km,database,oss.wnd,cable_220,cable_400,ks)
        #aft candidate transformers
        xfo_aft220s=set_oss_transformers(220.0,wnd_power,database,ks)
        xfo_aft400s=set_oss_transformers(400.0,wnd_power,database,ks)
        xfo220_cost=0
        xfo400_cost=0
        for xfo220 in xfo_aft220s
            xfo220_cost=xfo220_cost+xfo220.costs.ttl
        end
        for xfo400 in xfo_aft400s
            xfo400_cost=xfo400_cost+xfo400.costs.ttl
        end

        #platform
        plat_aft=platform()
        plat_aft.acdc="ac"
        plat_aft.mva=oss.mva
        plat_aft.wnd=oss.wnd
        plat_aft=cost_ac_platform(plat_aft,ks)#would need to be changed
        plat_aft=adjust_base_ac_platform(plat_aft,ks)

        #totals
        owpp2oss_connection220_cost=plat_aft.costs.ttl+xfo220_cost+cable_220.costs.grand_ttl
        owpp2oss_connection400_cost=plat_aft.costs.ttl+xfo400_cost+cable_400.costs.grand_ttl

        #totals and organize hv equipment
        push!(dict_ac_connection220,("plat_aft_ac"=>plat_aft))#aft platform
        push!(dict_ac_connection220,("xfo_aft"=>xfo_aft220s))#aft platform
        push!(dict_ac_connection220,("cable"=>cable_220))#cable
        push!(dict_ac_connection220,("cost"=>owpp2oss_connection220_cost[1]))#cost

        push!(dict_ac_connection400,("plat_aft_ac"=>plat_aft))#aft platform
        push!(dict_ac_connection400,("xfo_aft"=>xfo_aft400s))#aft platform
        push!(dict_ac_connection400,("cable"=>cable_400))#cable
        push!(dict_ac_connection400,("cost"=>owpp2oss_connection400_cost[1]))#cost
    end
    dict_ac_connection=Dict{String,Dict{String,Any}}()
    push!(dict_ac_connection,("220.0"=>dict_ac_connection220))
    push!(dict_ac_connection,("400.0"=>dict_ac_connection400))
    return dict_ac_connection
end
############# checks for Mid -point compensation ####################
#check if HVAC mid-point compensation is better option
function check4beter_mpc_ac_cables(mva,km,database,wnd,cable_220,cable_400,ks)
    if (database["bits"]["mpc_ac"]==true)
        cable_220mpc=deepcopy(hvac_cable(mva,km/2,wnd,database["cables"]["220.0"][string(mva)],ks))
        cable_400mpc=deepcopy(hvac_cable(mva,km/2,wnd,database["cables"]["400.0"][string(mva)],ks))
        cable_220mpc.mpc_ac=true
        cable_400mpc.mpc_ac=true
        #platform
        #compensation
        #NOTE CHANGE PLAT COST 5
        #cable_220mpc.plat.mva=mva
        cable_220mpc.plat.mva=mva/2
        cable_220mpc.plat.wnd=wnd
        cable_220mpc.plat=cost_ac_platform(cable_220mpc.plat,ks)#would need to be changed
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
    end
    #return best option
    return cable_220,cable_400
end

function owpp_acdc_to_oss(km,owpp,oss_location,database,ks,hv_only)
    owpp2oss_connections=owpp_ac_to_oss(km,owpp,oss_location,database,ks,hv_only)
    owpp2oss_connections_dc=owpp_dc_to_oss(km,owpp,oss_location,database,ks,hv_only)
    push!(owpp2oss_connections,("300.0"=>owpp2oss_connections_dc))
    return owpp2oss_connections
end

function oss_acdc_to_oss(km,oss,wnd_power,database,ks)
    owpp2oss_connections=oss_hvac_to_oss(oss,wnd_power,km,database,ks)
    owpp2oss_connections_dc=oss_dc_to_oss(km,oss,wnd_power,database,ks)
    push!(owpp2oss_connections,("300.0"=>owpp2oss_connections_dc))
    return owpp2oss_connections
end

#creates the a basic OSS node
function create_oss_node(mva,wnd,connect_bus,ocn)
    oss=bus()
    oss.mva=mva
    oss.wnd=wnd
    oss.node=deepcopy(connect_bus)
    ocn.num=ocn.num+1
    oss.node.num=deepcopy(ocn.num)
    return oss
end


function owpps2oss(hv_connections,oss_location,circ,ocn,hv_only)
    #find base pcc connection to start
    #record the transmission voltage and copy variables
    database=ocn.database
    ks=get_Cost_Data()
    owpp2oss_connections=Array{Dict{String,Dict{String,Any}},1}()
    #create an OSS bus
    push!(circ.mog,deepcopy(create_oss_node(circ.mva,circ.wnd,oss_location,ocn)))
    #sort OWPP connections to be made from farthest away to closest
    lengths_owpps=Array{Tuple,1}()
    for owp in circ.owpps
        push!(lengths_owpps,(euclidian_distance(owp.node.xy,oss_location.xy),owp))
    end
    sort!(lengths_owpps,by = x-> x[1],rev=true)
    #find the possible connections
    for l_o in lengths_owpps
        push!(owpp2oss_connections,deepcopy(owpp_acdc_to_oss(l_o[1],l_o[2],oss_location,database,ks,hv_only)))
    end
        #Make 220kV and 400kV variations
    circ=find_optimal_circuit(lengths_owpps,owpp2oss_connections,hv_connections,ocn,circ,ks)
    return circ
end

function find_optimal_circuit(lengths_owpps,aft_cons,for_cons,ocn,topo,ks)
    in_nodes=build_input_nodes(aft_cons)
    out_nodes=build_output_nodes(for_cons)
    candidates=build_candidates(in_nodes,out_nodes,ocn.database,ks)
    sol=solve_tnep(in_nodes,candidates,out_nodes,ks)
    topo=build_solution_circuit_owpps2mog(lengths_owpps,topo,sol,aft_cons,for_cons,ocn,ks)
    return topo
end

function find_optimal_circuitO2O(lengths_owpps,aft_cons,for_cons,ocn,topo,ks)
    in_nodes=build_input_nodes(aft_cons)
    out_nodes=build_output_nodes(for_cons)
    candidates=build_candidates(in_nodes,out_nodes,ocn.database,ks)
    sol=solve_tnep(in_nodes,candidates,out_nodes,ks)
    topo=build_solution_circuit_oss2mog(lengths_owpps,topo,sol,aft_cons,for_cons,ocn,ks)
    return topo
end

function build_solution_circuit_oss2mog(lengths_owpps,topo,sol,aft_cons,for_cons,ocn,ks)
    #Attach output circuit
    wnd_power=create_wind_and_power_dict()
    #Attach inputs
    inputs=sol["input"]
    for kv in ["66.0","220.0","400.0","300.0"]
        for input in inputs[kv]
            topo,wnd_power=best_oss2oss(aft_cons[input[2]][kv],topo,lengths_owpps[input[2]][2],ocn.database,ks,wnd_power,ocn)
        end
    end

    #build platform, transformers, converters and lay export cable
    output=sol["output"]
    if haskey(output,"220.0")
        topo.mog[1].kV=220.0
    elseif haskey(output,"400.0")
        topo.mog[1].kV=400.0
    elseif haskey(output,"300.0")
        topo.mog[1].kV=300.0
        #find converters
        topo.mog[1].conv=set_oss_convertersO2O(topo.mog[1].kV,wnd_power,ocn.database,ks)
    else
        println("Mondai ga arimasu!")
    end
    #find transformers
    topo.mog[1].xfmrs=set_oss_transformers(topo.mog[1].kV,wnd_power,ocn.database,ks)
    #set platform structure
    topo=set_mog_platformO2O(topo,ks)
    return topo
end

function build_solution_circuit_owpps2mog(lengths_owpps,topo,sol,aft_cons,for_cons,ocn,ks)
    #Attach output circuit
    wnd_power=create_wind_and_power_dict()
    #Attach inputs
    inputs=sol["input"]
    for kv in ["66.0","220.0","400.0","300.0"]
        for input in inputs[kv]
            topo,wnd_power=best_owpp2oss(aft_cons[input[2]][kv],topo,lengths_owpps[input[2]][2],ocn.database,ks,wnd_power,ocn)
        end
    end

    #build platform, transformers, converters and lay export cable
    output=sol["output"]
    if haskey(output,"220.0")
        topo.mog[1].kV=220.0
        push!(topo.PCCcbls,for_cons["220.0"]["cable"])
        if (for_cons["220.0"]["xfo_for"].num>0)
            push!(topo.pcc.xfmrs,for_cons["220.0"]["xfo_for"])
        end
    elseif haskey(output,"400.0")
        topo.mog[1].kV=400.0
        push!(topo.PCCcbls,for_cons["400.0"]["cable"])
        if (for_cons["400.0"]["xfo_for"].num>0)
            push!(topo.pcc.xfmrs,for_cons["400.0"]["xfo_for"])
        end
    elseif haskey(output,"300.0")
        topo.mog[1].kV=300.0
        push!(topo.PCCcbls,for_cons["300.0"]["cable"])
        push!(topo.pcc.xfmrs,for_cons["300.0"]["xfo_for"])
        push!(topo.pcc.conv,for_cons["300.0"]["conv_for"])
        #find converters
        topo.mog[1].conv=set_oss_converters(topo.mog[1].kV,wnd_power,ocn.database,ks)
    else
        println("Mondai ga arimasu!")
    end
    #find transformers
    topo.mog[1].xfmrs=set_oss_transformers(topo.mog[1].kV,wnd_power,ocn.database,ks)
    #set platform structure
    topo=set_mog_platform(topo,ks)
    #connect PCC cable
    push!(topo.PCCcbls[1].path,topo.mog[1].node)
    push!(topo.PCCcbls[1].path,topo.pcc.node)
    topo=total_circuit_cost(topo)
    return topo
end

function set_mog_platform(topo,ks)
    #find platform
    #platform
    plat_base=platform()
    plat_base.acdc="ac"
    #NOTE CHANGE PLAT COST 3
    #plat_base.mva=topo.PCCcbls[1].mva
    #for
    xfo_mva=0
    for xfo in topo.mog[1].xfmrs
        xfo_mva=xfo_mva+xfo.mva
    end
    no_xfo_mva=topo.PCCcbls[1].mva-xfo_mva
    plat_base.mva=xfo_mva+no_xfo_mva/2
    ########################
    plat_base.wnd=topo.PCCcbls[1].wnd
    plat_base=cost_ac_platform(plat_base,ks)#need to be changed
    if (length(topo.mog[1].conv)>0)
        #hvdc platform
        plat_dc=platform()
        plat_dc.acdc="dc"
        plat_dc.mva=topo.mog[1].conv[1].mva
        plat_dc.wnd=topo.mog[1].conv[1].wnd
        plat_dc=cost_dc_platform(plat_dc,ks)
        plat_dc=adjust_base_dc_platform(plat_dc,ks)
        topo.mog[1].plat=[plat_base,plat_dc]
    else
        plat_base=adjust_base_ac_platform(plat_base,ks)
        topo.mog[1].plat=[plat_base]
    end
    return topo
end

function set_mog_platformO2O(topo,ks)
    #find platform
    #platform
    plat_base=platform()
    plat_base.acdc="ac"
    #NOTE CHANGE PLAT COST 2
    #plat_base.mva=topo.mva
    #for
    xfo_mva=0
    for xfo in topo.mog[1].xfmrs
        xfo_mva=xfo_mva+xfo.mva
    end
    no_xfo_mva=topo.mva-xfo_mva
    plat_base.mva=xfo_mva+no_xfo_mva/2
    ########################
    plat_base.wnd=topo.wnd
    plat_base=cost_ac_platform(plat_base,ks)#would need to be changed
    if (length(topo.mog[1].conv)>0)
        #hvdc platform
        plat_dc=platform()
        plat_dc.acdc="dc"
        plat_dc.mva=topo.mog[1].conv[1].mva
        plat_dc.wnd=topo.mog[1].conv[1].wnd
        plat_dc=cost_dc_platform(plat_dc,ks)
        topo.mog[1].plat=[plat_base,plat_dc]
    else
        topo.mog[1].plat=[plat_base]
    end
    return topo
end

function solve_tnep(in_nodes,candidates,out_nodes,ks)
    m = Model(Cbc.Optimizer)
    set_optimizer_attribute(m, "logLevel", 0)
    ac_oss=ks.pac_f+ks.pac_f*ks.opx_pl*npv_years()
    dc_oss=ks.pdc_h+ks.pdc_h*ks.opx_pl*npv_years()+ks.conv_d+ks.conv_d*ks.opx_co*npv_years()

    #input nodes
    in_tup66=Array{Tuple,1}()
    in_powers66=Array{Float64,1}()
    in_costs66=Array{Float64,1}()
    for in_node in in_nodes["66.0"]
        push!(in_tup66,(in_node["mva"],in_node["line"]))
        push!(in_powers66,in_node["mva"])
        push!(in_costs66,in_node["cost"])
    end
    @variable(m, x66[1:length(in_powers66)],Bin)
    input_lines66=Dict{String,VariableRef}()
    for i=1:1:length(in_tup66)
        push!(input_lines66,(string(in_tup66[i][2])=>x66[i]))
    end

    in_tup220=Array{Tuple,1}()
    in_powers220=Array{Float64,1}()
    in_costs220=Array{Float64,1}()
    for in_node in in_nodes["220.0"]
        push!(in_tup220,(in_node["mva"],in_node["line"]))
        push!(in_powers220,in_node["mva"])
        push!(in_costs220,in_node["cost"])
    end
    @variable(m, x220[1:length(in_powers220)],Bin)
    input_lines220=Dict{String,VariableRef}()
    for i=1:1:length(in_tup220)
        push!(input_lines220,(string(in_tup220[i][2])=>x220[i]))
    end

    in_tup400=Array{Tuple,1}()
    in_powers400=Array{Float64,1}()
    in_costs400=Array{Float64,1}()
    for in_node in in_nodes["400.0"]
        push!(in_tup400,(in_node["mva"],in_node["line"]))
        push!(in_powers400,in_node["mva"])
        push!(in_costs400,in_node["cost"])
    end
    @variable(m, x400[1:length(in_powers400)],Bin)
    input_lines400=Dict{String,VariableRef}()
    for i=1:1:length(in_tup400)
        push!(input_lines400,(string(in_tup400[i][2])=>x400[i]))
    end

    in_tup300=Array{Tuple,1}()
    in_powers300=Array{Float64,1}()
    in_costs300=Array{Float64,1}()
    for in_node in in_nodes["300.0"]
        push!(in_tup300,(in_node["mva"],in_node["line"]))
        push!(in_powers300,in_node["mva"])
        push!(in_costs300,in_node["cost"])
    end
    @variable(m, x300[1:length(in_powers300)],Bin)
    input_lines300=Dict{String,VariableRef}()
    for i=1:1:length(in_tup300)
        push!(input_lines300,(string(in_tup300[i][2])=>x300[i]))
    end


    number_of_lines=findmax([length(x66),length(x220),length(x400),length(x300)])[1]
    input_lines=Dict{String,Vector{VariableRef}}()
    for i=1:1:number_of_lines
        in_line=Vector{VariableRef}()
        if haskey(input_lines66,string(i))
            push!(in_line,input_lines66[string(i)])
        end
        if haskey(input_lines220,string(i))
            push!(in_line,input_lines220[string(i)])
        end
        if haskey(input_lines400,string(i))
            push!(in_line,input_lines400[string(i)])
        end
        if haskey(input_lines300,string(i))
            push!(in_line,input_lines300[string(i)])
        end
        push!(input_lines,(string(i)=>in_line))
    end
    for i=1:1:number_of_lines
        @constraint(m, sum(input_lines[string(i)])==1)
    end



    #output nodes
    if (length(out_nodes["220.0"])>0)
        out_powers220=[out_nodes["220.0"]["mva"]]
        out_costs220=[out_nodes["220.0"]["cost"]]
    else
        out_powers220=[]
        out_costs220=[]
    end
    if (length(out_nodes["400.0"])>0)
        out_powers400=[out_nodes["400.0"]["mva"]]
        out_costs400=[out_nodes["400.0"]["cost"]]
    else
        out_powers400=[]
        out_costs400=[]
    end
    if (length(out_nodes["300.0"])>0)
        out_powers300=[out_nodes["300.0"]["mva"]]
        out_costs300=[out_nodes["300.0"]["cost"]]
    else
        out_powers300=[]
        out_costs300=[]
    end
    @variable(m, z220[1:1],Bin)
    @variable(m, z400[1:1],Bin)
    @variable(m, z300[1:1],Bin)
    y=Vector{VariableRef}()
    candidate_costs=Vector{Float32}()

    @variable(m, dc,Bin)
    dc_bins=Vector{VariableRef}()
    dc_powers=Vector{Float32}()
    for kv in ["66.0220.0","66.0400.0","66.0300.0","220.0220.0","220.0400.0","220.0300.0","400.0220.0","400.0400.0","400.0300.0","300.0220.0","300.0400.0","300.0300.0"]
        for candidate in candidates[kv]
            if (kv=="66.0220.0")
                push!(x66,@variable(m,binary=true))
                push!(z220,x66[length(x66)])
                push!(in_powers66,-1*candidate["mva"])
                push!(out_powers220,-1*candidate["mva"])
                push!(y,x66[length(x66)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="66.0400.0")
                push!(x66,@variable(m,binary=true))
                push!(z400,x66[length(x66)])
                push!(in_powers66,-1*candidate["mva"])
                push!(out_powers400,-1*candidate["mva"])
                push!(y,x66[length(x66)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="66.0300.0")
                push!(x66,@variable(m,binary=true))
                push!(z300,x66[length(x66)])
                push!(in_powers66,-1*candidate["mva"])
                push!(out_powers300,-1*candidate["mva"])
                push!(y,x66[length(x66)])
                push!(candidate_costs,candidate["cost"])
                push!(dc_bins,x66[length(x66)])
                push!(dc_powers,candidate["mva"])
            elseif (kv=="220.0220.0")
                push!(x220,@variable(m,binary=true))
                push!(z220,x220[length(x220)])
                push!(in_powers220,-1*candidate["mva"])
                push!(out_powers220,-1*candidate["mva"])
                push!(y,x220[length(x220)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="220.0400.0")
                push!(x220,@variable(m,binary=true))
                push!(z400,x220[length(x220)])
                push!(in_powers220,-1*candidate["mva"])
                push!(out_powers400,-1*candidate["mva"])
                push!(y,x220[length(x220)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="220.0300.0")
                push!(x220,@variable(m,binary=true))
                push!(z300,x220[length(x220)])
                push!(in_powers220,-1*candidate["mva"])
                push!(out_powers300,-1*candidate["mva"])
                push!(y,x220[length(x220)])
                push!(candidate_costs,candidate["cost"])
                push!(dc_bins,x220[length(x220)])
                push!(dc_powers,candidate["mva"])
            elseif (kv=="400.0220.0")
                push!(x400,@variable(m,binary=true))
                push!(z220,x400[length(x400)])
                push!(in_powers400,-1*candidate["mva"])
                push!(out_powers220,-1*candidate["mva"])
                push!(y,x400[length(x400)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="400.0400.0")
                push!(x400,@variable(m,binary=true))
                push!(z400,x400[length(x400)])
                push!(in_powers400,-1*candidate["mva"])
                push!(out_powers400,-1*candidate["mva"])
                push!(y,x400[length(x400)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="400.0300.0")
                push!(x400,@variable(m,binary=true))
                push!(z300,x400[length(x400)])
                push!(in_powers400,-1*candidate["mva"])
                push!(out_powers300,-1*candidate["mva"])
                push!(y,x400[length(x400)])
                push!(candidate_costs,candidate["cost"])
                push!(dc_bins,x400[length(x400)])
                push!(dc_powers,candidate["mva"])
            elseif (kv=="300.0220.0")
                push!(x300,@variable(m,binary=true))
                push!(z220,x300[length(x300)])
                push!(in_powers300,-1*candidate["mva"])
                push!(out_powers220,-1*candidate["mva"])
                push!(y,x300[length(x300)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="300.0400.0")
                push!(x300,@variable(m,binary=true))
                push!(z400,x300[length(x300)])
                push!(in_powers300,-1*candidate["mva"])
                push!(out_powers400,-1*candidate["mva"])
                push!(y,x300[length(x300)])
                push!(candidate_costs,candidate["cost"])
            elseif (kv=="300.0300.0")
                push!(x300,@variable(m,binary=true))
                push!(z300,x300[length(x300)])
                push!(in_powers300,-1*candidate["mva"])
                push!(out_powers300,-1*candidate["mva"])
                push!(y,x300[length(x300)])
                push!(candidate_costs,candidate["cost"])
            else
                println("No match?!?")
            end
        end
    end
    @constraint(m, sum(in_powers66[i]*x66[i] for i=1:length(in_powers66))==0.0)
    @constraint(m, sum(in_powers220[i]*x220[i] for i=1:length(in_powers220))==0.0)
    @constraint(m, sum(in_powers400[i]*x400[i] for i=1:length(in_powers400))==0.0)
    @constraint(m, sum(in_powers300[i]*x300[i] for i=1:length(in_powers300))==0.0)

    @constraint(m, sum(out_powers220[i]*z220[i] for i=1:length(out_powers220))==0.0)
    @constraint(m, sum(out_powers400[i]*z400[i] for i=1:length(out_powers400))==0.0)
    @constraint(m, sum(out_powers300[i]*z300[i] for i=1:length(out_powers300))==0.0)

    @constraint(m, z220[1]+z400[1]+z300[1]==1)

    dc_ttl=sum(dc_powers)
    @constraint(m, sum(dc_powers[i]*dc_bins[i] for i=1:length(dc_powers))-dc_ttl*dc<=0.0)

    @objective(m, Min,dc*dc_oss-dc*ac_oss+ac_oss+sum(candidate_costs[i]*y[i] for i=1:length(candidate_costs))+sum(in_costs66[i]*x66[i] for i=1:length(in_costs66))+sum(in_costs220[i]*x220[i] for i=1:length(in_costs220))+sum(in_costs400[i]*x400[i] for i=1:length(in_costs400))+sum(in_costs300[i]*x300[i] for i=1:length(in_costs300))+sum(out_costs220[i]*z220[i] for i=1:length(out_costs220))+sum(out_costs400[i]*z400[i] for i=1:length(out_costs400))+sum(out_costs300[i]*z300[i] for i=1:length(out_costs300)))
    optimize!(m)
    _in=Dict{String,Vector{Tuple}}()
    in_66=Array{Tuple,1}()
    for (i,b) in enumerate(x66[1:length(in_tup66)])
        if (JuMP.value(b)>0.0)
            push!(in_66,(in_tup66[i][1],in_tup66[i][2]))
        end
    end
    push!(_in,("66.0"=>in_66))

    in_220=Array{Tuple,1}()
    for (i,b) in enumerate(x220[1:length(in_tup220)])
        if (JuMP.value(b)>0.0)
            push!(in_220,(in_tup220[i][1],in_tup220[i][2]))
        end
    end
    push!(_in,("220.0"=>in_220))
    in_costs220
    in_400=Array{Tuple,1}()
    for (i,b) in enumerate(x400[1:length(in_tup400)])
        if (JuMP.value(b)>0.0)
            push!(in_400,(in_tup400[i][1],in_tup400[i][2]))
        end
    end
    push!(_in,("400.0"=>in_400))

    in_300=Array{Tuple,1}()
    for (i,b) in enumerate(x300[1:length(in_tup300)])
        if (JuMP.value(b)>0.0)
            push!(in_300,(in_tup300[i][1],in_tup300[i][2]))
        end
    end
    push!(_in,("300.0"=>in_300))

    _out=Dict{String,Float64}()
    out_220=Float64
    if (JuMP.value(z220[1])>0.0)
        out_220=out_powers220[1]
        push!(_out,("220.0"=>out_220))
    end

    out_400=Float64
    if (JuMP.value(z400[1])>0.0)
        out_400=out_powers400[1]
        push!(_out,("400.0"=>out_400))
    end


    out_300=Float64
    if (JuMP.value(z300[1])>0.0)
        out_300=out_powers300[1]
        push!(_out,("300.0"=>out_300))
    end


    _sol=Dict{String,Dict{String,Any}}()
    push!(_sol,"input"=>_in)
    push!(_sol,"output"=>_out)
    return _sol
end

function build_candidates(in_nodes,out_nodes,database,ks)
    ps_dict=Dict{String,Array{Any,1}}()
    for kv in ["66.0","220.0","400.0","300.0"]
        ps=Array{Any,1}()
        power_combos=ps=Array{Any,1}()
        for in_node in in_nodes[kv]
            push!(ps,deepcopy(in_node["mva"]))
        end
        power_sets=collect(combinations(ps))
        for power_set in power_sets
            push!(power_combos,sum(power_set))
        end
        unique!(power_combos)
        push!(ps_dict,(kv=>deepcopy(power_combos)))
    end
    power_cost=Dict{String,Any}()
    candidates=Dict{String,Array{Dict{String,Any}}}()
    push!(candidates,("66.0220.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("66.0400.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("66.0300.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("220.0220.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("220.0400.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("220.0300.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("400.0220.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("400.0400.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("400.0300.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("300.0220.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("300.0400.0"=>Array{Dict{String,Any},1}()))
    push!(candidates,("300.0300.0"=>Array{Dict{String,Any},1}()))
    for kv_in in ["66.0","220.0","400.0","300.0"]
        for p in ps_dict[kv_in]
            push!(power_cost,("mva"=>deepcopy(p)))
            #candidate equipment
            #transformer
            xfo=deepcopy(database["transformers"][string(p)]["offshore"])
            #converter
            conv=deepcopy(database["converters"][string(p)]["offshore"])
            conv.wnd=xfo.wnd
            conv=cost_hvdc_oss(conv,ks)
            #platforms
            plat_ac=platform()
            plat_ac.mva=xfo.mva
            plat_ac.wnd=xfo.wnd
            plat_dc=deepcopy(plat_ac)
            plat_ac=cost_ac_platform(plat_ac,ks)
            plat_dc=cost_dc_platform(plat_dc,ks)
            for kv_out in ["220.0","400.0","300.0"]
                if (kv_out=="220.0" || kv_out=="400.0")#AC feeder
                    if (kv_in==kv_out)#No transformation
                        #push!(power_cost,("cost"=>0.0))
                        #NOTE CHANGE PLAT COST 1
                        push!(power_cost,("cost"=>plat_ac.costs.ttl/2))
                        #push!(power_cost,("cost"=>plat_ac.costs.ttl))
                        push!(candidates[kv_in*kv_out],deepcopy(power_cost))
                    elseif (kv_in=="300.0")#DC input
                    else#transformer required
                        push!(power_cost,("cost"=>plat_ac.costs.ttl+xfo.costs.ttl))
                        push!(candidates[kv_in*kv_out],deepcopy(power_cost))
                    end
                else
                    if (kv_in=="300.0")
                        #push!(power_cost,("cost"=>0.0))
                        #NOTE CHANGE PLAT COST 4
                        push!(power_cost,("cost"=>plat_ac.costs.ttl/2))
                        #push!(power_cost,("cost"=>plat_ac.costs.ttl))
                        push!(candidates[kv_in*kv_out],deepcopy(power_cost))
                    else
                        push!(power_cost,("cost"=>plat_ac.costs.ttl+plat_dc.costs.ttl+xfo.costs.ttl+conv.costs.ttl))
                        push!(candidates[kv_in*kv_out],deepcopy(power_cost))
                    end
                end
            end
        end
    end
    return candidates
end

function build_output_nodes(for_cons)
    out_220=Dict{String,Any}()
    out_400=Dict{String,Any}()
    out_300=Dict{String,Any}()
    output_nodes=Dict{String,Dict{String,Any}}()
    if (haskey(for_cons["220.0"],"cable"))
        if (for_cons["220.0"]["cost"] != Inf)
            push!(out_220,("cost"=>deepcopy(for_cons["220.0"]["cost"])))
            push!(out_220,("mva"=>deepcopy(for_cons["220.0"]["cable"].mva)))
        end
    end
    if (haskey(for_cons["400.0"],"cable"))
        if (for_cons["400.0"]["cost"] != Inf)
            push!(out_400,("cost"=>deepcopy(for_cons["400.0"]["cost"])))
            push!(out_400,("mva"=>deepcopy(for_cons["400.0"]["cable"].mva)))
        end
    end
    if (haskey(for_cons["300.0"],"cable"))
        if (for_cons["300.0"]["cost"] != Inf)
            push!(out_300,("cost"=>deepcopy(for_cons["300.0"]["cost"])))
            push!(out_300,("mva"=>deepcopy(for_cons["300.0"]["cable"].mva)))
        end
    end
    push!(output_nodes,("220.0"=>out_220))
    push!(output_nodes,("400.0"=>out_400))
    push!(output_nodes,("300.0"=>out_300))
    return output_nodes
end

function build_input_nodes(aft_cons)
    #Build input nodal structures
    in_66=Dict{String,Any}()
    in_220=Dict{String,Any}()
    in_400=Dict{String,Any}()
    in_300=Dict{String,Any}()
    in_66s=Dict{String,Any}()
    in_220s=Dict{String,Any}()
    in_400s=Dict{String,Any}()
    in_66s=Array{Dict{String,Any},1}()
    in_220s=Array{Dict{String,Any},1}()
    in_400s=Array{Dict{String,Any},1}()
    in_300s=Array{Dict{String,Any},1}()
    input_nodes=Dict{String,Array{Dict{String,Any}}}()
    for (i,ac) in enumerate(aft_cons)
        if (haskey(ac,"66.0"))
            if (haskey(ac["66.0"],"cable"))
                if (ac["66.0"]["cost"] != Inf)
                    push!(in_66,("cost"=>deepcopy(ac["66.0"]["cost"])))
                    push!(in_66,("mva"=>deepcopy(ac["66.0"]["cable"].mva)))
                    push!(in_66,("line"=>deepcopy(i)))
                    push!(in_66s,deepcopy(in_66))
                end
            end
        end
        if (haskey(ac,"220.0"))
            if (haskey(ac["220.0"],"cable"))
                if (ac["220.0"]["cost"] != Inf)
                    push!(in_220,("cost"=>deepcopy(ac["220.0"]["cost"])))
                    push!(in_220,("mva"=>deepcopy(ac["220.0"]["cable"].mva)))
                    push!(in_220,("line"=>deepcopy(i)))
                    push!(in_220s,deepcopy(in_220))
                end
            end
        end
        if (haskey(ac,"400.0"))
            if (haskey(ac["400.0"],"cable"))
                if (ac["400.0"]["cost"] != Inf)
                    push!(in_400,("cost"=>deepcopy(ac["400.0"]["cost"])))
                    push!(in_400,("mva"=>deepcopy(ac["400.0"]["cable"].mva)))
                    push!(in_400,("line"=>deepcopy(i)))
                    push!(in_400s,deepcopy(in_400))
                end
            end
        end
        if (haskey(ac,"300.0"))
            if (haskey(ac["300.0"],"cable"))
                if (ac["300.0"]["cost"] != Inf)
                    push!(in_300,("cost"=>deepcopy(ac["300.0"]["cost"])))
                    push!(in_300,("mva"=>deepcopy(ac["300.0"]["cable"].mva)))
                    push!(in_300,("line"=>deepcopy(i)))
                    push!(in_300s,deepcopy(in_300))
                end
            end
        end
    end
    push!(input_nodes,("66.0"=>in_66s))
    push!(input_nodes,("220.0"=>in_220s))
    push!(input_nodes,("400.0"=>in_400s))
    push!(input_nodes,("300.0"=>in_300s))
    return input_nodes
end

########################### OWPP to PCC direct connection ######################
#finds the possible connections for MOG to PCC
function mog2pcc_possibilities_wMOG(oss_location,circ,ocn)
    ks=get_Cost_Data()
    km=euclidian_distance(oss_location.xy,circ.pcc.node.xy)
    hv_connections=optimal_mog2pcc_wMOG(km,circ.mva,circ.pcc.kV,circ.wnd,ocn.database)
    return hv_connections
end

function mog2pcc_possibilities_noMOG(oss_location,circ,ocn)
    ks=get_Cost_Data()
    km=euclidian_distance(oss_location.xy,circ.pcc.node.xy)
    hv_connections=optimal_mog2pcc_noMOG(km,circ.mva,circ.pcc.kV,circ.wnd,ocn.database)
    return hv_connections
end

function owpp2pcc(circ,ocn)
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
    hv_connections=optimal_mog2pcc_wMOG(km,circ.owpps[1].mva,circ.pcc.kV,circ.owpps[1].wnd,ocn.database)
    cheapest=argmin([hv_connections["300.0"]["cost"],hv_connections["220.0"]["cost"],hv_connections["400.0"]["cost"]])
    hv_connection=[hv_connections["300.0"],hv_connections["220.0"],hv_connections["400.0"]][cheapest]
    owpp2pcc_connection_cost=findmin([mv_connection["cost"],hv_connection["cost"]])
    owpp2pcc_connection=[mv_connection,hv_connection][owpp2pcc_connection_cost[2]]
    #Build the zero length MVAC cable
    mv_cable0=add_0km_66kVcable_owpp2oss(circ.owpps[1].mva,circ.owpps[1].wnd,circ.owpps[1].node,ocn.database,ks)
    #Build the best solution circuit
    if (haskey(owpp2pcc_connection,"plat_aft_dc"))#HVDC connection
        #build the offshore MOG
        push!(circ.mog,aft_connection_DC_equipment(owpp2pcc_connection,circ.owpps[1].node,ocn))
        #lay the cable
        push!(mv_cable0.path,circ.mog[1].node)
        push!(circ.MVcbls,mv_cable0)
        push!(owpp2pcc_connection["cable"].path,circ.pcc.node)
        push!(circ.PCCcbls,owpp2pcc_connection["cable"])
        #add the converter onshore
        push!(circ.pcc.conv,owpp2pcc_connection["conv_for"])
        #add the transformer onshore
        push!(circ.pcc.xfmrs,owpp2pcc_connection["xfo_for"])
        ocn.owpps[circ.base.node.num].kv2pcc=owpp2pcc_connection["cable"].elec.volt
    elseif (haskey(owpp2pcc_connection,"plat_aft_ac"))#HVAC
        #Build MOG
        push!(circ.mog,aft_connection_AC_equipment(owpp2pcc_connection,circ.owpps[1].node,ocn))
        #lay the cable
        push!(mv_cable0.path,circ.mog[1].node)
        push!(circ.MVcbls,mv_cable0)
        push!(owpp2pcc_connection["cable"].path,circ.pcc.node)
        push!(circ.PCCcbls,owpp2pcc_connection["cable"])
        #add the transformer onshore
        if (owpp2pcc_connection["xfo_for"].costs.ttl>0.0)
            push!(circ.pcc.xfmrs,owpp2pcc_connection["xfo_for"])
        end
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

#Builds a 0km size cable from the OWPP to the OSS
function add_0km_66kVcable_owpp2oss(mva,wnd,source_node,database,ks)
    mv_cable0=deepcopy(mvac_cable(mva,0,wnd,database["cables"]["66.0"][string(mva)],ks))
    push!(mv_cable0.path,source_node)
    return mv_cable0
end

function aft_connection_AC_equipment(connection,connect_bus,ocn)
    #create an OSS bus
    oss=bus()
    oss.kV=connection["cable"].elec.volt
    oss.mva=connection["cable"].mva
    oss.wnd=connection["cable"].wnd
    oss.node=deepcopy(connect_bus)
    ocn.num=ocn.num+1
    oss.node.num=deepcopy(ocn.num)
    #build the platform
    push!(oss.plat,connection["plat_aft_ac"])
    #place the transformers
    push!(oss.xfmrs,connection["xfo_aft"])
    #attach the cable
    push!(connection["cable"].path,oss.node)
    return oss
end


function aft_connection_AC_equipmentO2O(connection,oss,ocn)
    #create an OSS bus
    oss.kV=connection["cable"].elec.volt
    #build the platform
    oss.plat=[]
    push!(oss.plat,connection["plat_aft_ac"])
    #place the transformers
    oss.xfmrs=[]
    oss.xfmrs=connection["xfo_aft"]
    #attach the cable
    push!(connection["cable"].path,oss.node)
    return oss
end

function aft_connection_DC_equipment(connection,connect_bus,ocn)
    oss=bus()
    oss.kV=connection["cable"].elec.volt
    oss.mva=connection["cable"].mva
    oss.wnd=connection["cable"].wnd
    oss.node=deepcopy(connect_bus)
    ocn.num=ocn.num+1
    oss.node.num=deepcopy(ocn.num)
    #build the platform
    push!(oss.plat,connection["plat_aft_ac"])
    push!(oss.plat,connection["plat_aft_dc"])
    #place the transformers
    push!(oss.xfmrs,connection["xfo_aft"])
    #place the converters
    push!(oss.conv,connection["conv_aft"])
    #attach the cable
    push!(connection["cable"].path,oss.node)
    return oss
end

function aft_connection_DC_equipmentO2O(connection,oss,ocn)
    oss.kV=connection["cable"].elec.volt
    #build the platform
    oss.plat=[]
    push!(oss.plat,connection["plat_aft_ac"])
    push!(oss.plat,connection["plat_aft_dc"])
    #place the transformers
    oss.xfmrs=[]
    oss.xfmrs=connection["xfo_aft"]
    #place the converters
    oss.conv=[]
    push!(oss.conv,connection["conv_aft"])
    #attach the cable
    push!(connection["cable"].path,oss.node)
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
#check
function optimal_mog2pcc_noMOG(km,mva,kv,wnd,database)
    ks=get_Cost_Data()
    #HVAC connection - mpc handled within HVAC
    hvac_connections=hvac_canditate_mog2pcc_noMOG(km,mva,kv,wnd,database,ks)#returns dictionary: [platform,trans_aft,cable,trans_for]
    #HVDC connection
    dc_connection=hvdc_canditate_mog2pcc_noMOG(km,mva,kv,wnd,database,ks)#returns dictionary: [platform,trans_aft,conv_aft,cable,conv_for,trans_for]
    push!(hvac_connections,("300.0"=>dc_connection))

    return hvac_connections
end

#check
function optimal_mog2pcc_wMOG(km,mva,kv,wnd,database)
    ks=get_Cost_Data()
    #HVAC connection - mpc handled within HVAC
    hvac_connections=hvac_canditate_mog2pcc_wMOG(km,mva,kv,wnd,database,ks)#returns dictionary: [platform,trans_aft,cable,trans_for]
    #HVDC connection
    dc_connection=hvdc_canditate_mog2pcc_wMOG(km,mva,kv,wnd,database,ks)#returns dictionary: [platform,trans_aft,conv_aft,cable,conv_for,trans_for]
    push!(hvac_connections,("300.0"=>dc_connection))

    return hvac_connections
end

function hvac_canditate_mog2pcc_wMOG(km,mva,kv_forward,wnd,database,ks)
    #if HVAC is selected find optimal equipments
    #buil forward most equipment
    dict_ac_connection=hvac_canditate_mog2pcc_noMOG(km,mva,kv_forward,wnd,database,ks)
    if (database["bits"]["hvac"]==true)
        #add the aft platform and transformer
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


        #store object in a dictionary
        push!(dict_ac_connection["220.0"],("plat_aft_ac"=>plat_aft))#plat_aft
        push!(dict_ac_connection["220.0"],("xfo_aft"=>xfo_aft))#xfo_aft
        push!(dict_ac_connection["400.0"],("plat_aft_ac"=>plat_aft))#plat_aft
        push!(dict_ac_connection["400.0"],("xfo_aft"=>xfo_aft))#xfo_aft

        dict_ac_connection["220.0"]["cost"]=dict_ac_connection["220.0"]["cost"]+plat_aft.costs.ttl+xfo_aft.costs.ttl
        dict_ac_connection["400.0"]["cost"]=dict_ac_connection["400.0"]["cost"]+plat_aft.costs.ttl+xfo_aft.costs.ttl
    end
    return dict_ac_connection
end

#HVAC connection from an MOG to a PCC
function hvac_canditate_mog2pcc_noMOG(km,mva,kv_forward,wnd,database,ks)
    #initialize dummy set
    dict_ac_connection220=Dict{String,Any}()
    dict_ac_connection400=Dict{String,Any}()
    push!(dict_ac_connection220,("cost"=>Inf))#xfo_for
    push!(dict_ac_connection400,("cost"=>Inf))#xfo_for
    #if HVAC is selected find optimal equipments
    if (database["bits"]["hvac"]==true)
        #build candidate cables
        cable_220=deepcopy(hvac_cable(mva,km,wnd,database["cables"]["220.0"][string(mva)],ks))
        cable_400=deepcopy(hvac_cable(mva,km,wnd,database["cables"]["400.0"][string(mva)],ks))
        cable_220.costs.grand_ttl=cable_220.costs.ttl
        cable_400.costs.grand_ttl=cable_400.costs.ttl
        #if mid point compensation is checked look for cheaper cable options
        if (database["bits"]["mpc_ac"]==true)
            cable_220,cable_400=check4beter_mpc_ac_cables(mva,km,database,wnd,cable_220,cable_400,ks)
        end

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
        ac_connection220_cost=cable_220.costs.grand_ttl+xfo_for220.costs.ttl
        ac_connection400_cost=cable_400.costs.grand_ttl+xfo_for400.costs.ttl
        #store object in a dictionary
        push!(dict_ac_connection220,("cable"=>cable_220))#cable
        push!(dict_ac_connection220,("xfo_for"=>xfo_for220))#xfo_for
        push!(dict_ac_connection220,("cost"=>ac_connection220_cost))#xfo_for

        #store object in a dictionary
        push!(dict_ac_connection400,("cable"=>cable_400))#cable
        push!(dict_ac_connection400,("xfo_for"=>xfo_for400))#xfo_for
        push!(dict_ac_connection400,("cost"=>ac_connection400_cost))#xfo_for

        dict_ac_connection=Dict{String,Dict{String,Any}}()
        push!(dict_ac_connection,("220.0"=>dict_ac_connection220))
        push!(dict_ac_connection,("400.0"=>dict_ac_connection400))
    end
    return dict_ac_connection
end

#hvdc mog to pcc candidate connection
function hvdc_canditate_mog2pcc_noMOG(km,mva,kv_forward,wnd,database,ks)
    dict_dc_connection=Dict{String,Any}()
    push!(dict_dc_connection,("cost"=>Inf))#conv_for
    if (database["bits"]["hvdc"]==true)
        cable_dc=deepcopy(hvdc_cable(mva,km,wnd,database["cables"]["300.0"][string(mva)],ks))
        cable_dc.costs.grand_ttl=cable_dc.costs.ttl

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
        cost_dc_connection=cable_dc.costs.grand_ttl+conv_for.costs.ttl+xfo_for.costs.ttl
        push!(dict_dc_connection,("cable"=>cable_dc))#cable
        push!(dict_dc_connection,("conv_for"=>conv_for))#conv_for
        push!(dict_dc_connection,("xfo_for"=>xfo_for))#conv_for
        push!(dict_dc_connection,("cost"=>cost_dc_connection))#conv_for
    end

    return dict_dc_connection
end

#hvdc mog to pcc candidate connection
function hvdc_canditate_mog2pcc_wMOG(km,mva,kv_forward,wnd,database,ks)
    dict_dc_connection=hvdc_canditate_mog2pcc_noMOG(km,mva,kv_forward,wnd,database,ks)
    if (database["bits"]["hvdc"]==true)
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

        #total cost and store in a dictionary
        push!(dict_dc_connection,("plat_aft_ac"=>plat_aft_ac))#plat_aft
        push!(dict_dc_connection,("plat_aft_dc"=>plat_aft_dc))#plat_aft
        push!(dict_dc_connection,("xfo_aft"=>xfo_aft))#conv_aft
        push!(dict_dc_connection,("conv_aft"=>conv_aft))#conv_aft
        dict_dc_connection["cost"]=dict_dc_connection["cost"]+plat_aft_ac.costs.ttl+plat_aft_dc.costs.ttl+xfo_aft.costs.ttl+conv_aft.costs.ttl
    end
    return dict_dc_connection
end

function best_oss2oss(cheapest_connection_2mog,circ,oss,database,ks,wnd_power,ocn)
    if (cheapest_connection_2mog["cable"].elec.volt==300.0)
        #DC transfer
        #build OSS
        push!(circ.oss,aft_connection_DC_equipmentO2O(cheapest_connection_2mog,oss,ocn))
        #lay the cable
        push!(cheapest_connection_2mog["cable"].path,circ.mog[1].node)
        push!(circ.O2Ocbls,cheapest_connection_2mog["cable"])
        #track the wind
        wnd_power=track_wind_and_power(wnd_power,cheapest_connection_2mog["cable"])
    elseif (cheapest_connection_2mog["cable"].elec.volt!=66.0)
        #Builtd the OSS
        push!(circ.oss,aft_connection_AC_equipmentO2O(cheapest_connection_2mog,oss,ocn))
        #lay the cable
        push!(cheapest_connection_2mog["cable"].path,circ.mog[1].node)
        push!(circ.O2Ocbls,cheapest_connection_2mog["cable"])
        #track the wind
        wnd_power=track_wind_and_power(wnd_power,cheapest_connection_2mog["cable"])
    else
    end
    return circ,wnd_power
end

function best_owpp2oss(cheapest_connection_2oss,circ,owp,database,ks,wnd_power,ocn)
    if (cheapest_connection_2oss["cable"].elec.volt==300.0)
        #DC transfer
        #build OSS
        push!(circ.oss,aft_connection_DC_equipment(cheapest_connection_2oss,owp.node,ocn))
        #lay the cable
        mv_cable0=add_0km_66kVcable_owpp2oss(owp.mva,owp.wnd,owp.node,database,ks)
        push!(mv_cable0.path,deepcopy(circ.oss[length(circ.oss)].node))
        push!(circ.MVcbls,mv_cable0)
        push!(cheapest_connection_2oss["cable"].path,circ.mog[1].node)
        push!(circ.HVcbls,cheapest_connection_2oss["cable"])
        #track the wind
        wnd_power=track_wind_and_power(wnd_power,cheapest_connection_2oss["cable"])
    elseif (cheapest_connection_2oss["cable"].elec.volt!=66.0)
        #Builtd the OSS
        push!(circ.oss,aft_connection_AC_equipment(cheapest_connection_2oss,owp.node,ocn))
        #lay the cable
        mv_cable0=add_0km_66kVcable_owpp2oss(owp.mva,owp.wnd,owp.node,database,ks)
        push!(mv_cable0.path,circ.oss[length(circ.oss)].node)
        push!(circ.MVcbls,mv_cable0)
        push!(cheapest_connection_2oss["cable"].path,circ.mog[1].node)
        push!(circ.HVcbls,cheapest_connection_2oss["cable"])
        #track the wind
        wnd_power=track_wind_and_power(wnd_power,cheapest_connection_2oss["cable"])
    else
        push!(cheapest_connection_2oss["cable"].path,circ.mog[1].node)
        push!(circ.MVcbls,cheapest_connection_2oss["cable"])
        wnd_power=track_wind_and_power(wnd_power,circ.MVcbls[length(circ.MVcbls)])
    end
    return circ,wnd_power
end

function set_oss_transformers(hv,wnd_power,database,ks)
    xfo_array=Array{transformer,1}()
    for kv in ["66.0","220.0","400.0"]
        if (wnd_power[kv]["mva"]>0.0 && hv != tryparse(Float32,kv))
            xfo=deepcopy(database["transformers"][string(wnd_power[kv]["mva"])]["offshore"])
            xfo.wnd=find_netWind(wnd_power[kv]["wind"])
            xfo.elec.lv=tryparse(Float32,kv)
            xfo.elec.hv=hv
            xfo=cost_xfo_oss(xfo,ks)
            push!(xfo_array,xfo)
        end
    end
    return xfo_array
end

function set_oss_converters(hv,wnd_power,database,ks)
    mva=wnd_power["66.0"]["mva"]+wnd_power["220.0"]["mva"]+wnd_power["400.0"]["mva"]
    if (mva>0.0)
        wnd=find_netWind([find_netWind(wnd_power["66.0"]["wind"]),find_netWind(wnd_power["220.0"]["wind"]),find_netWind(wnd_power["400.0"]["wind"])])
        conv=deepcopy(database["converters"][string(mva)]["offshore"])
        conv.wnd=wnd
        conv=cost_hvdc_oss(conv,ks)
        conv=adjust_base_hvdc_offshore_converter(conv,ks)
        return [conv]
    else
        return []
    end
end

function set_oss_convertersO2O(hv,wnd_power,database,ks)
    mva=wnd_power["66.0"]["mva"]+wnd_power["220.0"]["mva"]+wnd_power["400.0"]["mva"]
    if (mva>0.0)
        wnd=find_netWind([find_netWind(wnd_power["66.0"]["wind"]),find_netWind(wnd_power["220.0"]["wind"]),find_netWind(wnd_power["400.0"]["wind"])])
        conv=deepcopy(database["converters"][string(mva)]["offshore"])
        conv.wnd=wnd
        conv=cost_hvdc_oss(conv,ks)
        return [conv]
    else
        return []
    end
end

function create_wind_and_power_dict()
    volt_wnd_power=Dict{String,Dict{String,Any}}()
    for kv in ["66.0","220.0","400.0","300.0"]
        mva_wnd_dict=Dict{String,Any}()
        push!(mva_wnd_dict,("mva"=>0.0))#cost
        push!(mva_wnd_dict,("wind"=>Array{wind,1}()))#cost
        push!(volt_wnd_power,(kv=>deepcopy(mva_wnd_dict)))
    end
    return volt_wnd_power
end

function track_wind_and_power(volt_wnd_power,cbl)
    if (cbl.elec.volt==300.0)
        volt_wnd_power["300.0"]["mva"]=volt_wnd_power["300.0"]["mva"]+cbl.mva
        push!(volt_wnd_power["300.0"]["wind"],cbl.wnd)
    elseif (cbl.elec.volt==400.0)
        volt_wnd_power["400.0"]["mva"]=volt_wnd_power["400.0"]["mva"]+cbl.mva
        push!(volt_wnd_power["400.0"]["wind"],cbl.wnd)
    elseif (cbl.elec.volt==220.0)
        volt_wnd_power["220.0"]["mva"]=volt_wnd_power["220.0"]["mva"]+cbl.mva
        push!(volt_wnd_power["220.0"]["wind"],cbl.wnd)
    else
        volt_wnd_power["66.0"]["mva"]=volt_wnd_power["66.0"]["mva"]+cbl.mva
        push!(volt_wnd_power["66.0"]["wind"],cbl.wnd)
    end
    return volt_wnd_power
end

############################### Below is deprecated ############################

#=
#takes the best MOG to PCC connection and loads it into a circuit
function mog2pcc(oss_location,hv_connection,circ,ocn)
    if (haskey(hv_connection,"plat_aft_dc"))#HVDC connection
        #build the offshore MOG
        push!(circ.mog,aft_connection_DC_equipment(hv_connection,oss_location,ocn))
        #lay the cable
        push!(hv_connection["cable"].path,circ.mog[1].node)
        push!(hv_connection["cable"].path,circ.pcc.node)
        push!(circ.PCCcbls,hv_connection["cable"])
        #add the converter onshore
        push!(circ.pcc.conv,hv_connection["conv_for"])
        #add the transformer onshore
        push!(circ.pcc.xfmrs,hv_connection["xfo_for"])
    else#HVAC
        #Build MOG
        push!(circ.mog,aft_connection_AC_equipment(hv_connection,oss_location,ocn))
        #lay the cable
        push!(hv_connection["cable"].path,circ.pcc.node)
        push!(circ.PCCcbls,hv_connection["cable"])
        #add the transformer onshore
        if (hv_connection["xfo_for"].costs.ttl>0.0)
            push!(circ.pcc.xfmrs,hv_connection["xfo_for"])
        end
    end
    circ=total_circuit_cost(circ)
    return circ
end
=#
#=
function best_mog2pcc(connection,connect_bus,ocn)
    #create an OSS bus
    oss=bus()
    oss.kV=connection["cable"].elec.volt
    oss.mva=connection["cable"].mva
    oss.wnd=connection["cable"].wnd
    oss.node=deepcopy(connect_bus)
    ocn.num=ocn.num+1
    oss.node.num=deepcopy(ocn.num)
    #build the platform
    push!(oss.plat,connection["plat_aft_ac"])
    #place the transformers
    push!(oss.xfmrs,connection["xfo_aft"])
    #attach the cable
    push!(connection["cable"].path,oss.node)
    return oss
end
=#

########################### OWPP to OSS connection #############################
#finds the optimal connections for OWPPS an OSS
#This is circuit building
#this is a layout function
#this is only a sorting function - seperates scenarios AC-AC, AC-DC, DC-DC, DC-AC
#NOTE this si now to be tested, beware there will be many mistakes testing needs to be 100% before proceeding!!!!
#=function owpps2oss_setup(hv_connections,oss_location,circ,ocn)
    cheapest=argmin([hv_connections["300.0"]["cost"],hv_connections["220.0"]["cost"],hv_connections["400.0"]["cost"]])
    cheapest_connection_2pcc=[hv_connections["300.0"],hv_connections["220.0"],hv_connections["400.0"]][cheapest]
    cheapest_connection_voltage=cheapest_connection_2pcc["cable"].elec.volt
    #build the MOG
    push!(circ.mog,aft_connection_AC_equipment(cheapest_connection_2pcc,oss_location,ocn))

    #sort OWPP connections to be made from farthest away to closest
    lengths_owpps=Array{Tuple,1}()
    hvdc_check=false
    for owp in circ.owpps
        push!(lengths_owpps,(euclidian_distance(owp.node.xy,oss_location.xy),owp))
        if (owp.kv2pcc==300.0 && ocn.database["bits"]["hvdc"]==true)
            hvdc_check=true
        end
    end
    if (hvdc_check==true && cheapest_connection_2pcc["cable"].elec.volt!=300.0 && ocn.database["bits"]["hvdc"]==true)
        println("optimal_owpps2oss says: Feeder upgrade to HVDC may be needed! Circuit ID: "*string(circ.id))
        cheapest_connection_2pcc=hv_connections["300.0"]
    end
    sort!(lengths_owpps,by = x-> x[1],rev=true)
    return lengths_owpps,hvdc_check,cheapest_connection_2pcc,circ
end

function combine_owpps2oss_ac2ac(circ,hv_connections,owpp2oss_connections,wnd_power,database,ks,hv_only)
    circ220=deepcopy(circ)
    circ400=deepcopy(circ)
    hv_connections["400.0"]["cable"].path=[circ.mog[1].node,circ.pcc.node]
    hv_connections["220.0"]["cable"].path=[circ.mog[1].node,circ.pcc.node]
    circ400.mog[1].kV=400.0
    circ400.pcc.xfmrs=[]
    circ220.mog[1].kV=220.0
    circ220.pcc.xfmrs=[]
    if (hv_connections["400.0"]["xfo_for"].costs.ttl>0.0)
        push!(circ400.pcc.xfmrs,hv_connections["400.0"]["xfo_for"])
    end
    if (hv_connections["220.0"]["xfo_for"].costs.ttl>0.0)
        push!(circ220.pcc.xfmrs,hv_connections["220.0"]["xfo_for"])
    end
    circ400.PCCcbls=[hv_connections["400.0"]["cable"]]
    circ220.PCCcbls=[hv_connections["220.0"]["cable"]]
    #find transformers for each variation
    circ220.mog[1].xfmrs=set_oss_transformers(circ220.mog[1].kV,wnd_power,database,ks)
    circ220=total_circuit_cost(circ220)
    circ400.mog[1].xfmrs=set_oss_transformers(circ400.mog[1].kV,wnd_power,database,ks)
    circ400=total_circuit_cost(circ400)
    if (hv_only==true)
        circ220=check_4_cable_swaps(circ220,owpp2oss_connections,deepcopy(wnd_power),database,ks)
        circ400=check_4_cable_swaps(circ400,owpp2oss_connections,deepcopy(wnd_power),database,ks)
    end
    cheapest=findmin([circ220.cost,circ400.cost])
    circ=deepcopy([circ220,circ400][cheapest[2]])
    return circ
end



function combine_owpps2oss_dc2dc(circ,dc_connection,owpp2oss_connections,wnd_power,database,ks,hv_only)
    #build common PCC componebts
    circ.mog[1].kV=300.0
    circ.pcc.xfmrs=[]
    circ.pcc.conv=[]
    push!(circ.pcc.xfmrs,dc_connection["xfo_for"])
    push!(circ.pcc.conv,dc_connection["conv_for"])

    #transfer the cable
    dc_connection["cable"].path=[circ.mog[1].node,circ.pcc.node]
    circ.PCCcbls=[dc_connection["cable"]]

    #split into 2 copies
    circ300=deepcopy(circ)
    circ300.cost=Inf
    circOpt=deepcopy(circ)

    #find optimal upstream connections
    #find transformers
    circOpt.mog[1].xfmrs=set_oss_transformers(circOpt.mog[1].kV,wnd_power,database,ks)
    #find converters
    circOpt.mog[1].conv=set_oss_converters(circOpt.mog[1].kV,wnd_power,database,ks)
    #find platform
    #platform
    plat_base=platform()
    plat_base.acdc="ac"
    plat_base.mva=circOpt.PCCcbls[1].mva
    plat_base.wnd=circOpt.PCCcbls[1].wnd
    plat_base=cost_ac_platform(plat_base,ks)
    if (length(circOpt.mog[1].conv)>0)
        #hvdc platform
        plat_dc=platform()
        plat_dc.acdc="dc"
        plat_dc.mva=circOpt.mog[1].conv[1].mva
        plat_dc.wnd=circOpt.mog[1].conv[1].wnd
        plat_dc=cost_dc_platform(plat_dc,ks)
        plat_dc=adjust_base_dc_platform(plat_dc,ks)
        circOpt.mog[1].plat=[plat_base,plat_dc]
    else
        plat_base=adjust_base_ac_platform(plat_base,ks)
    end

    #Check if meshed DC grid is best
    if ((length(circOpt.MVcbls)==length(circOpt.oss)) && (length(circOpt.mog[1].conv)>0))
        for (index, oss) in enumerate(circ300.oss)
            #remove old platform
            circ300.oss[index].plat=[]
            #Build HVDC platform upstream
            #hvac base
            plat_base=platform()
            plat_base.acdc="ac"
            plat_base.mva=circ300.HVcbls[index].mva
            plat_base.wnd=circ300.HVcbls[index].wnd
            plat_base=cost_ac_platform(plat_base,ks)
            #hvdc platform
            plat_dc=platform()
            plat_dc.acdc="dc"
            plat_dc.mva=circ300.HVcbls[index].mva
            plat_dc.wnd=circ300.HVcbls[index].wnd
            plat_dc=cost_dc_platform(plat_dc,ks)
            plat_dc=adjust_base_dc_platform(plat_dc,ks)
            circ300.oss[index].plat=deepcopy([plat_base,plat_dc])

            #add converter to platform
            conv=deepcopy(database["converters"][string(circ300.HVcbls[index].mva)]["offshore"])
            conv.wnd=circ300.HVcbls[index].wnd
            conv=cost_hvdc_oss(conv,ks)
            conv=adjust_base_hvdc_offshore_converter(conv,ks)
            push!(circ300.oss[index].conv,conv)

            #replace cables
            owpp2oss_connections[index]["300.0"]["cable"].path=circ300.HVcbls[index].path
            circ300.HVcbls[index]=deepcopy(owpp2oss_connections[index]["300.0"]["cable"])
        end
        #platform
        plat_base=platform()
        plat_base.acdc="ac"
        plat_base.mva=circOpt.PCCcbls[1].mva
        plat_base.wnd=circOpt.PCCcbls[1].wnd
        plat_base=cost_ac_platform(plat_base,ks)
        plat_base=adjust_base_ac_platform(plat_base,ks)
        circ300.mog[1].plat=[plat_base]
        circ300=total_circuit_cost(circ300)
    end
    circOpt=total_circuit_cost(circOpt)
    cheapest=argmin([circ300.cost,circOpt.cost])
    circ=deepcopy([circ300,circOpt][cheapest])
    return circ
end

function check_4_cable_swaps(circ,connections,wnd_power,database,ks)
    circ_copy=deepcopy(circ)
    wnd_power_copy=deepcopy(wnd_power)
    for (index,hv_cbl) in enumerate(circ_copy.HVcbls)
        if (hv_cbl.elec.volt!=circ.PCCcbls[1].elec.volt)
            wnd_power_copy=modify_wind_and_power_dict(wnd_power_copy,hv_cbl,connections[index][string(circ.PCCcbls[1].elec.volt)]["cable"])
            connections[index][string(circ.PCCcbls[1].elec.volt)]["cable"].path=circ_copy.HVcbls[index].path
            circ_copy.HVcbls[index]=connections[index][string(circ.PCCcbls[1].elec.volt)]["cable"]
            circ_copy.mog[1].xfmrs=set_oss_transformers(circ_copy.mog[1].kV,wnd_power_copy,database,ks)
            circ_copy=total_circuit_cost(circ_copy)
            if (circ_copy.cost<circ.cost)
                circ=deepcopy(circ_copy)
                wnd_power=deepcopy(wnd_power_copy)
            else
                circ_copy=deepcopy(circ)
                wnd_power_copy=deepcopy(wnd_power)
            end
        end
    end
    return circ
end



function modify_wind_and_power_dict(volt_wnd_power,cbl_old,cbl_new)
    if (cbl_old.elec.volt==300.0)
        volt_wnd_power["300.0"]["mva"]=volt_wnd_power["300.0"]["mva"]-cbl_old.mva
        volt_wnd_power["300.0"]["wind"]=delete_wind_entry(cbl_old,volt_wnd_power["300.0"]["wind"])
    elseif (cbl_old.elec.volt==400.0)
        volt_wnd_power["400.0"]["mva"]=volt_wnd_power["400.0"]["mva"]-cbl_old.mva
        volt_wnd_power["400.0"]["wind"]=delete_wind_entry(cbl_old,volt_wnd_power["400.0"]["wind"])
    elseif (cbl_old.elec.volt==220.0)
        volt_wnd_power["220.0"]["mva"]=volt_wnd_power["220.0"]["mva"]-cbl_old.mva
        volt_wnd_power["220.0"]["wind"]=delete_wind_entry(cbl_old,volt_wnd_power["220.0"]["wind"])
    else
        volt_wnd_power["66.0"]["mva"]=volt_wnd_power["66.0"]["mva"]-cbl_old.mva
        volt_wnd_power["66.0"]["wind"]=delete_wind_entry(cbl_old,volt_wnd_power["66.0"]["wind"])
    end
    volt_wnd_power=track_wind_and_power(volt_wnd_power,cbl_new)
    return volt_wnd_power
end

#delete a wind entry of the array
function delete_wind_entry(cbl_old,wnd_array)
    for (index,wnd) in enumerate(wnd_array)
        if (wnd.delta+1e-5>cbl_old.wnd.delta && wnd.delta-1e-5<cbl_old.wnd.delta && wnd.lf+1e-5>cbl_old.wnd.lf && wnd.lf-1e-5<cbl_old.wnd.lf)
            wnd_array=deleteat!(wnd_array,index)
            break
        end
    end
    return wnd_array
end=#

#=function connection_HVACDCowpp2ossDC(km,mva,mog,database,ks)
    owpp2oss_connections=Dict{String,Dict{String,Any}}()
    owpp2oss_connectionsACDC=optimal_HVACowpp2ossDC(km,mva,mog,database,ks)
    owpp2oss_connectionsDCDC=optimal_DCowpp2ossDC(km,mva,mog,database,ks)
    if (owpp2oss_connectionsACDC["cost"]<owpp2oss_connectionsDCDC["cost"])
        owpp2oss_connectionsBest=owpp2oss_connectionsACDC
    else
        owpp2oss_connectionsBest=owpp2oss_connectionsDCDC
    end
    push!(owpp2oss_connections,("HVDC"=>owpp2oss_connectionsDCDC))
    push!(owpp2oss_connections,("Best"=>owpp2oss_connectionsBest))
    return owpp2oss_connections
end=#
#=
function connection_HVACowpp2ossACDC(km,owpp,mog,database,ks)
    #option 1 66-platform
    owpp2oss_connectionAC=Dict{String,Dict{String,Any}}()
    owpp2oss_connectionHVAC=Dict{String,Any}()
    owpp2oss_connectionMVAC=Dict{String,Any}()
    push!(owpp2oss_connectionHVAC,("cost"=>Inf))#cost
    push!(owpp2oss_connectionMVAC,("cost"=>Inf))#cost

    #option 2 66-OSS-HVAC-platform
    if (database["bits"]["hvac"]==true)
        owpp2oss_connectionHVAC=owpp_hvac_to_oss(owpp,km,mog.kV,database,ks)
        push!(owpp2oss_connectionHVAC["cable"].path,owpp.node)
        push!(owpp2oss_connectionHVAC["cable"].path,mog.node)
    else
        println("Error! connection_HVACowpp2ossACDC: hvac not selected.")
    end
    push!(owpp2oss_connectionAC,("HVAC"=>owpp2oss_connectionHVAC))
    push!(owpp2oss_connectionAC,("MVAC"=>owpp2oss_connectionMVAC))
    return owpp2oss_connectionAC
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
function combine_connections_ACDCowpp2ossDC(owpp2oss_connections,mog_circ,database,ks)
    #finds common equipment
    mog_circDC=combine_connections_DCowpp2ossDC(owpp2oss_connections["HVDC"],deepcopy(mog_circ),database,ks)
    mog_circBest=combine_connections_Bestowpp2ossDC(owpp2oss_connections["Best"],deepcopy(mog_circ),database)
    if (mog_circDC.cost<mog_circBest)
        return mog_circDC
    else
        return mog_circBest
    end
end


#finds common equipment of circuit for ACowpp2ossAC and ACowpp2ossDC
function combine_connections_ACowpp2ossACDC(owpp2oss_connections,mog_circ,database)
    mva66=0.0
    mva220=0.0
    mva400=0.0
    wnd66=Array{wind,1}()
    wnd220=Array{wind,1}()
    wnd400=Array{wind,1}()
    for owpp2oss_connection in owpp2oss_connections
        if (owpp2oss_connection["HVAC"]["cost"]<owpp2oss_connection["MVAC"]["cost"])
            mog_circ=connections_ACowpp2ossACDC(owpp2oss_connection["HVAC"],mog_circ)
            #record powers and winds per voltage level to build forward platform transformers
            if (owpp2oss_connection["HVAC"]["cable"].elec.volt==220.0)
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
        xfo66.elec.lv=mog_circ.MVcbls[1].elec.volt
        xfo66.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo66=cost_xfo_oss(xfo66,ks)
        push!(mog_circ.mog[1].xfmrs,xfo66)
    end
    if (mva220>0.0 && (mog_circ.mog[1].kV != 220.0))
        xfo220=deepcopy(database["transformers"][string(mva220)]["offshore"])
        xfo220.wnd=find_netWind(wnd220)
        xfo220.elec.lv=mog_circ.MVcbls[1].elec.volt
        xfo220.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo220=cost_xfo_oss(xfo220,ks)
        push!(mog_circ.mog[1].xfmrs,xfo220)
    end
    if (mva400>0.0 && (mog_circ.mog[1].kV != 400.0))
        xfo400=deepcopy(database["transformers"][string(mva400)]["offshore"])
        xfo400.wnd=find_netWind(wnd400)
        xfo400.elec.lv=mog_circ.MVcbls[1].elec.volt
        xfo400.elec.hv=mog_circ.PCCcbls[1].elec.volt
        xfo400=cost_xfo_oss(xfo400,ks)
        push!(mog_circ.mog[1].xfmrs,xfo400)
    end
    return mog_circ
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

function combine_connections_DCowpp2ossDC(owpp2oss_connections,mog_circ,database,ks)
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

function connections_DCowpp2ossDC(owpp2oss_connection,mog_circ)
    #build the offshore MOG
    push!(mog_circ.oss,aft_connection_DC_equipment(owpp2pcc_connection))
    #lay the cable
    push!(mog_circ.HVcbls,owpp2oss_connection["cable"])
    return
end

function combine_connections_Bestowpp2ossDC(owpp2oss_connections,mog_circ,database,ks)
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

#takes the results of each owpp to oss line and combines with the oss to pcc line to create a final circuit only for hvac/mvac connections to an HVAC MOG
function combine_connections_ACowpp2ossAC(owpp2oss_connections,mog_circ,database,ks)
    #finds common equipment
    mog_circ=combine_connections_ACowpp2ossACDC(owpp2oss_connections,mog_circ,database,ks)
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
=#
#=
#NOTE Nothing from below should be involved
############################################################################################
#finds best alternative from 66-platform, 66-OSS-HVAC-platform, 66-OSS-HVDC-platform
function optimal_HVACowpp2ossDC(km,owpp,mog,database,ks)
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




function owpps2oss(hv_connections,oss_location,circ,ocn,hv_only)
    #find base pcc connection to start
    #record the transmission voltage and copy variables
    database=ocn.database
    ks=get_Cost_Data()
    wnd_power=create_wind_and_power_dict()
    owpp2oss_connections=Array{Dict{String,Dict{String,Any}},1}()
    #record the transmission voltage and copy variables
    lengths_owpps,hvdc_check,cheapest_connection_2pcc,circ=owpps2oss_setup(hv_connections,oss_location,circ,ocn)

    if (cheapest_connection_2pcc["cable"].elec.volt==300.0)#HVDC Feeder
        #AC/DC to DC options handled
        #l_o=lengths_owpps[2]
        for l_o in lengths_owpps
            push!(owpp2oss_connections,deepcopy(owpp_acdc_to_oss(l_o[1],l_o[2],oss_location,database,ks,hv_only)))
            cheapest=argmin([owpp2oss_connections[length(owpp2oss_connections)]["66.0"]["cost"],owpp2oss_connections[length(owpp2oss_connections)]["220.0"]["cost"],owpp2oss_connections[length(owpp2oss_connections)]["400.0"]["cost"],owpp2oss_connections[length(owpp2oss_connections)]["300.0"]["cost"]])
            cheapest_connection_2oss=[owpp2oss_connections[length(owpp2oss_connections)]["66.0"],owpp2oss_connections[length(owpp2oss_connections)]["220.0"],owpp2oss_connections[length(owpp2oss_connections)]["400.0"],owpp2oss_connections[length(owpp2oss_connections)]["300.0"]][cheapest]
            #save the set of best connections and track the wind and power
            circ,wnd_power=best_owpp2oss(cheapest_connection_2oss,circ,l_o[2],database,ks,wnd_power,ocn)
        end
        circ=combine_owpps2oss_dc2dc(circ,cheapest_connection_2pcc,owpp2oss_connections,wnd_power,database,ks,hv_only)
    #end
    else#AC to AC option handled here
        #l_o=lengths_owpps[8]
        for l_o in lengths_owpps
            push!(owpp2oss_connections,deepcopy(owpp_ac_to_oss(l_o[1],l_o[2],oss_location,database,ks,hv_only)))
            cheapest=argmin([owpp2oss_connections[length(owpp2oss_connections)]["66.0"]["cost"],owpp2oss_connections[length(owpp2oss_connections)]["220.0"]["cost"],owpp2oss_connections[length(owpp2oss_connections)]["400.0"]["cost"]])
            cheapest_connection_2oss=[owpp2oss_connections[length(owpp2oss_connections)]["66.0"],owpp2oss_connections[length(owpp2oss_connections)]["220.0"],owpp2oss_connections[length(owpp2oss_connections)]["400.0"]][cheapest]
            #save the set of best connections and track the wind and power
            circ,wnd_power=best_owpp2oss(cheapest_connection_2oss,circ,l_o[2],database,ks,wnd_power,ocn)
        end
        #Make 220kV and 400kV variations
        circ=combine_owpps2oss_ac2ac(circ,hv_connections,owpp2oss_connections,wnd_power,database,ks,hv_only)
    end
    return circ
end
=#
