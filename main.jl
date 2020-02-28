
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

ocean_hv=load("tempFiles/data/solutions/hv_ocean8.jld2")["ocean"]
@time opt_compoundOSS(ocean_hv)
@time ocean_hvc = opt_readjust_circuits(ocean_hvc)#255 = 37.924910 seconds

ocean_hvc=load("tempFiles/data/solutions/compound_hv_ocean8.jld2")["ocean"]
ppf_testing(ocean_hvc)
ppf_equipment_OSS_MOG(ocean_hvc,ocean_hvc.circuits[28])
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
    ocean.pccs[2].node.xy.y=15
    number_of_owpps=length(ocean.owpps)
    #ocean=load("tempFiles/data/solutions/clean_ocean8.jld2")["ocean"]
    ppf_saveSystem(ocean,"clean_ocean"*string(number_of_owpps))

    @time ocean.circuits=opt_hvOSSplacement_zero_size(ocean,ocean.pccs[2])
    ppf_saveSystem(ocean,"hv_ocean"*string(number_of_owpps))
    @time opt_compoundOSS_zeroSize(ocean)
    ppf_saveSystem(ocean,"compound_hv_ocean"*string(number_of_owpps))
    @time ocean = opt_readjust_circuits(ocean)

    @time best_full_syss,ocean.circuits=opt_rollUp(ocean)
    ppf_saveSystem(ocean,"compound_hv_rolledUp_ocean"*string(number_of_owpps))
    ppf_saveCircuit(best_full_syss,"best_full_hv_systems"*string(number_of_owpps))

    ocean_hv_nodes=load("tempFiles/data/solutions/hv_ocean"*string(number_of_owpps)*".jld2")["ocean"].discretedom.nodes
    ocean=load("tempFiles/data/solutions/clean_ocean"*string(number_of_owpps)*".jld2")["ocean"]
    @time ocean.circuits=opt_mvOSSplacement_zero_size(ocean,ocean.owpps,ocean.pccs[2],ocean.discretedom.nodes)
    ocean_hv=load("tempFiles/data/solutions/hv_ocean"*string(number_of_owpps)*".jld2")["ocean"]
    ocean=opt_updateMVocean(ocean,ocean_hv,oe_length)
    ocean_hv=eez()
    @time ocean = opt_readjust_circuits(ocean)
    ppf_saveSystem(ocean,"mv_ocean"*string(number_of_owpps))

    @time opt_compoundOSS_zeroSize(ocean)
    @time ocean = opt_readjust_circuits(ocean)
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
ppf_equipment_OSS_MOG(ocean_hv,ocean_hv.circuits[234])


#################
ocean=load("tempFiles/data/solutions/clean_ocean8.jld2")["ocean"]
ocean_hv=load("tempFiles/data/solutions/hv_ocean8.jld2")["ocean"]
@time opt_compoundOSS(ocean_hv)
#################
ocean_cmv=load("tempFiles/data/solutions/compound_mv_rolledUp_readjusted_ocean8.jld2")["ocean"]
@time ocean_cmv = opt_readjust_circuits(ocean_cmv)
best_full_syss_cmv,ocean_cmv.circuits=opt_rollUp(ocean_cmv)
ppf_saveSystem(ocean_cmv,"compound_mv_rolledUp_readjusted_ocean8")
ppf_saveCircuit(best_full_syss_cmv,"best_full_syss_reajusted_mv8")

ocean_chv=load("tempFiles/data/solutions/compound_hv_rolledUp_readjusted_ocean8.jld2")["ocean"]
@time ocean_chv = opt_readjust_circuits(ocean_chv)
best_full_syss_chv,ocean_chv.circuits=opt_rollUp(ocean_chv)
ppf_saveSystem(ocean_chv,"compound_hv_rolledUp_readjusted_ocean8")
ppf_saveCircuit(best_full_syss_chv,"best_full_syss_reajusted_hv8")

ocean_cmv.circuits=optSysCompare(ocean_chv.circuits,ocean_cmv.circuits)
best_full_syss_mvhv,ocean_cmv.circuits=opt_rollUp(ocean_cmv)
ppf_saveSystem(ocean_cmv,"compound_mvhv_rolledUp_readjusted_ocean8")
ppf_saveCircuit(best_full_syss_mvhv,"best_full_syss_reajusted_mvhv8")

##################################################################################################
bsf_hv=load("tempFiles/data/circuits/circ_best_full_syss_reajusted_hv8.jld2")["circuits"]
bsf_mv=load("tempFiles/data/circuits/circ_best_full_syss_reajusted_mv8.jld2")["circuits"]
bsf_mvhv=load("tempFiles/data/circuits/circ_best_full_syss_reajusted_mvhv8.jld2")["circuits"]

ocean_chv=load("tempFiles/data/solutions/compound_hv_ocean8.jld2")["ocean"]
ocean_cmv=load("tempFiles/data/solutions/compound_mv_ocean8.jld2")["ocean"]
@time ocean_cmv = opt_readjust_circuits(ocean_cmv)
@time ocean_chv = opt_readjust_circuits(ocean_chv)
push!(bsf_hv,ocean_chv.circuits[255])
push!(bsf_mv,ocean_cmv.circuits[255])

ocean_chv=load("tempFiles/data/solutions/compound_hv_rolledUp_readjusted_ocean8.jld2")["ocean"]
ocean_cmv=load("tempFiles/data/solutions/compound_mv_rolledUp_readjusted_ocean8.jld2")["ocean"]
push!(bsf_hv,ocean_chv.circuits[255])
push!(bsf_mv,ocean_cmv.circuits[255])
bsf_mvhv=combineAndrank(bsf_mv,bsf_hv,bsf_mvhv)
ppf_equipment_OSS_MOG(ocean_cmv,bsf_mvhv[1])
for (i,c) in enumerate(bsf_mvhv)
    println(string(i)*" - "*string(c.cost))
end
##################################################################################################
ppf_equipment_OSS_MOG(ocean_chv,bsf_mvhv[1])
println(ocean_chv.circuits[1].cost+ocean_chv.circuits[2].cost+ocean_chv.circuits[4].cost+ocean_chv.circuits[8].cost+ocean_chv.circuits[16].cost+ocean_chv.circuits[32].cost+ocean_chv.circuits[64].cost+ocean_chv.circuits[128].cost)
ppf_equipment_OSS_MOG(ocean_chv,best_full_syss_chv[1])
ppf_equipment_OSS_MOG(ocean_cmv,best_full_syss_cmv[1])
combineAndrank(best_full_syss_cmv,best_full_syss_chv,best_full_syss_mvhv)
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
cost_bd(bsf_mvhv)
function cost_bd(ocnhv)
    c_loss=0
    c_cm=0
    c_eens=0
    c_cpx=0
    c_q=0

    x_cpx=0
    x_eens=0
    x_cm=0
    x_loss=0


    println("c_cpx c_loss c_eens c_cm c_q x_cpx x_loss x_eens x_cm")
    for (i,circ) in enumerate(ocnhv)
        for cbl in circ.owp_MVcbls
            c_cpx=c_cpx+cbl.costs.cbc
            c_q=c_q+cbl.costs.qc
            c_loss=c_loss+cbl.costs.rlc
            c_cm=c_cm+cbl.costs.cm
            c_eens=c_eens+cbl.costs.eens
        end
        for cbl in circ.owp_HVcbls

            c_cpx=c_cpx+cbl.costs.cbc
            c_q=c_q+cbl.costs.qc
            c_loss=c_loss+cbl.costs.rlc
            c_cm=c_cm+cbl.costs.cm
            c_eens=c_eens+cbl.costs.eens
            #println(c_cpx+c_q+c_loss+c_cm+c_eens)
        end
        for cbl in circ.oss2oss_cbls

            c_cpx=c_cpx+cbl.costs.cbc
            c_q=c_q+cbl.costs.qc
            c_loss=c_loss+cbl.costs.rlc
            c_cm=c_cm+cbl.costs.cm
            c_eens=c_eens+cbl.costs.eens
            #println(c_cpx+c_q+c_loss+c_cm+c_eens)
        end
        for cbl in circ.pcc_cbls
            c_cpx=c_cpx+cbl.costs.cbc
            c_q=c_q+cbl.costs.qc
            c_loss=c_loss+cbl.costs.rlc
            c_cm=c_cm+cbl.costs.cm
            c_eens=c_eens+cbl.costs.eens
            #println(c_cpx+c_q+c_loss+c_cm+c_eens)
        end
        for mog in circ.osss_mog
            x_cpx=x_cpx+mog.base_cost
            for x in mog.xfmrs
                x_cpx=x_cpx+x.costs.cpx
                x_loss=x_loss+x.costs.tlc
                x_cm=x_cm+x.costs.cm
                x_eens=x_eens+x.costs.eens
            end
        end
        for mog in circ.osss_owp
            x_cpx=x_cpx+mog.base_cost
            for x in mog.xfmrs
                x_cpx=x_cpx+x.costs.cpx
                x_loss=x_loss+x.costs.tlc
                x_cm=x_cm+x.costs.cm
                x_eens=x_eens+x.costs.eens
            end
        end
        println(string(i)*" "*string(c_cpx)*" "*string(c_loss)*" "*string(c_eens)*" "*string(c_cm)*" "*string(c_q)*" "*string(x_cpx)*" "*string(x_loss)*" "*string(x_eens)*" "*string(x_cm))
        c_loss=0
        x_loss=0
        c_cm=0
        x_cm=0
        c_eens=0
        x_eens=0
        c_cpx=0
        x_cpx=0
        c_q=0
    end

end

circ=bsf_mvhv[1]
cbl_count(bsf_mvhv)
lof_pnt2pnt_dist(ocean_cmv.owpps[1].node.xy,ocean_cmv.owpps[2].node.xy)_
ppf_equipment_OSS_MOG(ocean_cmv,bsf_mv[8])
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

l=26.26666136
km=l
S=2*250
kv=220
wa=wind[]
wo=wind()
wo.pu=zeros(Float32,8759)
wo.ce=zeros(Float32,8759)
wo.delta=0
wo.lf=0
push!(wa,ocean.owpps[2].wnd)
push!(wa,ocean.owpps[1].wnd)

wnd=opt_Wsum(deepcopy(wo),wa)
println(cstF_HvCblo2p(l,S,kv,wnd,cstD_cfs()).costs.ttl)

println(cstF_MvCbl(l,S,kv,ocean.owpps[1].wnd,cstD_cfs()).costs.ttl)
cstF_HvCblallKvo2p(l,S,wnd,ocean.finance)
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
