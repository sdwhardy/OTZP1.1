#using JuliaInterpreter
#using Debugger, ProxSDP
using DataFrames, XLSX, CSV, JLD2, FileIO
using StatsPlots, Plots
using Polynomials, TypedPolynomials, SpecialFunctions
using LinearAlgebra, Polyhedra, SetProg, MathOptInterface
using JuMP, Ipopt, Mosek, MosekTools

include("wind/struct.jl")
include("cost/struct.jl")#
include("eqp/struct.jl")#
include("optimization/struct.jl")#
include("layout/struct.jl")#
#include("topology/struct.jl")

include("cost/data.jl")#
include("eqp/data.jl")#

include("eens/functions.jl")
include("wind/functions.jl")
include("cost/functions.jl")#
include("layout/functions.jl")#
include("layout/functions0size.jl")#
include("eqp/functions.jl")#
include("post_process/functions.jl")#costs
include("astar/functions.jl")
include("optimization/functions.jl")
include("optimization/functions0size.jl")
include("optimization/elipse_functions.jl")
include("topology/functions.jl")




#=
    #with path finding
    @time ocean=lof_layoutEez_basis()
    @time ocean.owpps=lof_order2Pcc(ocean,ocean.pccs[2])
    @time ocean=lof_layoutEez_expand(ocean,ocean.pccs[2])
    @time ocean.circuits=opt_hvOSSplacement(ocean,ocean.pccs[2])
    @time ocean.circuits=opt_mvOSSplacement(ocean,ocean.owpps,ocean.pccs[2])
    ppf_testing(ocean)
=#
function main()
    #without path finding, zero size
    @time ocean=lof_layoutEez_basis()
    @time ocean.owpps=lof_order2Pcc(ocean,ocean.pccs[2])
    @time ocean,oe_length=lof_layoutEez_expand_testing(ocean,ocean.pccs[2])
    number_of_owpps=length(ocean.owpps)
    #ocean=load("tempFiles/data/solutions/clean_ocean8.jld2")["ocean"]
    ppf_saveSystem(ocean,"clean_ocean"*string(number_of_owpps))

    @time ocean.circuits=opt_hvOSSplacement_zero_size(ocean,ocean.pccs[2])
    ppf_saveSystem(ocean,"hv_ocean"*string(number_of_owpps))
    @time opt_compoundOSS(ocean)
    ppf_saveSystem(ocean,"compound_hv_ocean"*string(number_of_owpps))
    @time best_full_syss,ocean.circuits=opt_rollUp(ocean)
    ppf_saveSystem(ocean,"compound_hv_rolledUp_ocean"*string(number_of_owpps))
    ppf_saveCircuit(best_full_syss,"best_full_hv_systems"*string(number_of_owpps))

    ocean_hv_nodes=load("tempFiles/data/solutions/hv_ocean"*string(number_of_owpps)*".jld2")["ocean"].discretedom.nodes
    ocean=load("tempFiles/data/solutions/clean_ocean"*string(number_of_owpps)*".jld2")["ocean"]
    @time ocean.circuits=opt_mvOSSplacement_zero_size(ocean,ocean.owpps,ocean.pccs[2],ocean_hv_nodes)
    ocean_hv=load("tempFiles/data/solutions/hv_ocean"*string(number_of_owpps)*".jld2")["ocean"]
    ocean=opt_updateMVocean(ocean,ocean_hv,oe_length)
    ocean_hv=eez()
    ppf_saveSystem(ocean,"mv_ocean"*string(number_of_owpps))

    @time opt_compoundOSS(ocean)
    ppf_saveSystem(ocean,"compound_mv_ocean"*string(number_of_owpps))
    @time best_full_syss,ocean.circuits=opt_rollUp(ocean)
    ppf_saveSystem(ocean,"compound_mv_rolledUp_ocean"*string(number_of_owpps))
    ppf_saveCircuit(best_full_syss,"best_full_mv_systems"*string(number_of_owpps))
    ocean_hv=load("tempFiles/data/solutions/compound_hv_rolledUp_ocean"*string(number_of_owpps)*".jld2")["ocean"]

    ocean.circuits=optSysCompare(ocean_hv.circuits,ocean.circuits)
    best_full_syss,ocean.circuits=opt_rollUp(ocean)
    ppf_saveSystem(ocean,"compound_hvmv_rolledUp_ocean"*string(number_of_owpps))
    ppf_saveCircuit(best_full_syss,"best_full_hvmv_systems"*string(number_of_owpps))

    return ocean,best_full_syss
end
@time ocean,best_full_syss=main()
ppf_equipment_OSS_MOG(ocean,best_full_syss[1])



#################
ocean_hv=load("tempFiles/data/solutions/hv_ocean8.jld2")["ocean"]
@time opt_compoundOSS(ocean_hv)
#################

ocean_chv=load("tempFiles/data/solutions/compound_hv_ocean8.jld2")["ocean"]

ocean_mv=load("tempFiles/data/solutions/mv_ocean8.jld2")["ocean"]
ocean_cmv=load("tempFiles/data/solutions/compound_mv_ocean8.jld2")["ocean"]
bc_cmv=load("tempFiles/data/circuits/circ_best_full_hvmv_systems8.jld2")["circuits"]

ppf_equipment_OSS_MOG(ocn,oss_system)
ppf_equipment_OSS_MOG(ocean_cmv,bc_cmv[2])
ppf_testing(ocean_cmv)
function cbl_count(ocnhv)
    cbles=0
    oss=0
    for (i,circ) in enumerate(ocnhv)
        cbles=cbles+length(circ.pcc_cbls)+length(circ.oss2oss_cbls)+length(circ.owp_HVcbls)+length(circ.owp_MVcbls)
        oss=oss+length(circ.osss_mog)+length(circ.osss_owp)
        println(string(i)*" - "*string(oss))
        oss=0
    end
    println(cbles)
    println(oss)
end
cbl_count(bc_cmv)

ppf_equipment_OSS_MOG(oceanvv,best_full_syss[1])
ppf_equipment_OSS_MOG(oceanhv,best_full_syss[1])
    #combining/refinnig solution

    #ppf_saveSystem(ocean,"HV_comp_ocean")
    @time best_full_sysshv,oceanhv.circuits=opt_rollUp(oceanhv)
    #ocean=lof_layoutEez_expand_testing(ocean,ocean.pccs[2])
    #ppf_layout_testing(oceanMV)

    ppf_equipment_OSS_MOG(oceanhv,oceanhv.circuits[255])
    for (i,crc) in enumerate(best_full_syss)
        println(string(i)*" - "*string(length(crc.osss_mog)+length(crc.osss_owp)))
    end
    ppf_equipment(ocean,oss_system)
    oceanHV=load("tempFiles/data/solutions/HV_comp_ocean.jld2")["ocean"]
    ocn_l=load("tempFiles/data/solutions/lowCB.jld2")["ocean"]
    ocn_h=load("tempFiles/data/solutions/highCB.jld2")["ocean"]

oceanhv.circuits=optSysCompare(oceanmv.circuits,oceanhv.circuits)

    ppf_printOcnXY_cables(ocn_l,ocn_l.circuits[28])
    printLines(ocean)
    gr()
    gui()

    ppf_equipment(ocean,ocean.circuits[3])
    ppf_equipment_OSS_MOG(ocean,ocean.circuits[135])


    ppf_saveCircuit(best_full_syss,"best_full_syss_ocean")
    ppf_saveSystem(ocean,"best_full_syss_ocean")
#ocean.discretedom.nodes[85]
    #start=ocean.owpps[3].node
    #goal=ocean.pccs[1].node
    #path=as_Astar(start,goal,ocean.discretedom.nodes)
    #ppf_printOcnGPS(ocean)
    #ppf_printOcnXY(ocean,path,start)

    #ppf_printOcnXY(ocean)
    #return ocean
end
main()
plotly()
plot([1,2,3,4,5],[1,2,3,4,5])

ppf_testing(ocean)

main()
ocean.circuits[15]
sea.circuits[14]

l=5
km=l
S=1500
kv=66
cstF_HvCblallKvo2p(l,S,wnd,ocean.finance)
cb=cstF_MvCbl(l,S,kv,wnd,cstD_cfs())
kv=132
cbo=cstF_HvCblo2o(l,S,kv,wnd,cstD_cfs())
kv=220
cb1=cstF_HvCblo2p(l,S,kv,ocean.owpps[1].wnd,ocean.finance)
kv=400
cb2=cstF_HvCblo2o(l,S,kv,ocean.owpps[1].wnd,ocean.finance)
cstF_xfo_oss(S,ocean.owpps[1].wnd,ocean.finance)
S=100
cstF_xfo_oss(S,ocean.owpps[1].wnd,ocean.finance)
cstF_xfo_pcc(S,ocean.circuits[1].osss_mog[1].wnd,ocean.finance)
cstF_xfo_pcc(power_sum,wind_sum,ocean.finance)


println(cstF_MvCbl3366(5,250,ocean.owpps[8].wnd,cstD_cfs()).costs.ttl)
cstF_HvCblallKvo2p(12,750,ocean.owpps[8].wnd,cstD_cfs(),ocean.pccs[2])
chv,xhv=cstF_HvCblallKvo2o(12,750,ocean.owpps[6].wnd,cstD_cfs(),220)
cstF_xfo_oss(750,ocean.owpps[6].wnd,cstD_cfs(),66,220)
