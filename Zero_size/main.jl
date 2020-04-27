using DataFrames, XLSX, CSV, JLD2, FileIO
using StatsPlots, Plots
using Polynomials, TypedPolynomials, SpecialFunctions
using LinearAlgebra, Polyhedra, SetProg, MathOptInterface
using JuMP, Ipopt, Mosek, MosekTools

include("struct.jl")#
include("economics/functions.jl")#
include("layout/functions.jl")#
include("optimization/functions.jl")#
include("post_process/functions.jl")#



function main()
    ocean=lof_layoutEez_basis()
    ocean.owpps=lof_order2Pcc(ocean,ocean.pccs[2])
    ocean=lof_layoutEez_expand(ocean,ocean.pccs[2])
    ocean.hv_circuits=opt_hvOSSplacement(ocean,ocean.pccs[2])
    ocean.mv_circuits=opt_mvOSSplacement(ocean,ocean.pccs[2])

    #ocean.mv_circuits=opt_updateMVocean(ocean)
    #ocean_clean=deepcopy(ocean)
    #ocean=ocean_clean
    ocean.hv_circuits,ocean.hvc_pct= opt_readjust_circuits(ocean,ocean.hv_circuits,ocean.hvc_pct)
    ocean.mv_circuits,ocean.mvc_pct = opt_readjust_circuits(ocean,ocean.mv_circuits,ocean.mvc_pct)
    #ocean.mv_circuits=ocean.hv_circuits
    #ocean.mv_circuits,ocean.hv_circuits=check_ids(ocean.mv_circuits,ocean.hv_circuits)
    #ppf_saveSystem(ocean,"clean_ocean7")
    #ocean.hv_circuits=ocean.mv_circuits
    #ocean_old=deepcopy(ocean)
    #ocean=ocean_old
    @time circuits=opt_compoundOSS(ocean)
    #@time circuits,circuitsHV,circuitsMV=opt_compoundOSSMVHV(ocean)
#c[224]=deepcopy(circuits)
#ch=deepcopy(circuitsHV)
#cm=deepcopy(circuitsMV)
    circuitsHV,ocean.hvc_pct= opt_readjust_circuits(ocean,circuitsHV,ocean.hvc_pct)
    circuitsMV,ocean.hvc_pct= opt_readjust_circuits(ocean,circuitsMV,ocean.hvc_pct)
    #@time circuits,ocean.hvc_pct=opt_readjust_circuits(ocean,circuits,ocean.hvc_pct)
    #circuits_old0=deepcopy(circuits)
    #circuits=circuits_old0
    #ppf_saveSystem(circuits,"bsf_mvhv2")127,247
    #ppf_saveSystem(ocean,"ocean")
    #best_fullHV_syss,circuits=opt_rollUp(ocean,circuits)
    #circuits_old2=deepcopy(circuits)
    #best_full_syss,circuits=opt_rollUp_partial(circuits, ceil(Int,length(circuits)/2))
    #best_full_syssHV,circuitsHV=opt_rollUp_partial(circuitsHV, 127)
    #best_full_syssMV[1],circuitsMV=opt_rollUp_partial(circuitsMV, 127)
    #best_fullMV_syss,ocean.mv_circuits=opt_rollUp(ocean,ocean.mv_circuits)
    #mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
    #best_fullMVHV_syss,mvhv_circuits=opt_rollUp(ocean,mvhv_circuits)
    #bsf_mvhv=combineAndrank(best_fullMV_syss,best_fullHV_syss,best_fullMVHV_syss)
    #ppf_saveCircuit(bsf_mvhv,"bsf_mvhv")
    @time FR,circuits=opt_rollUp(circuits)
    return ocean, circuits,FR
end
#ppf_saveCircuit(circuits,"circuits24_4")
#circuits=deepcopy(bsf_mvhv2[1:254])
#ocean11=ocean
@time ocean, bsf_mvhv2,fr=main()
FR=opt_rollUp(bsf_mvhv)
ppf_saveCircuit(fr,"fastHVonly_1inQ")

ppf_saveSystem(ocean,"clean_ocean7")
bsf_mvhv2,ocean.hvc_pct= opt_readjust_circuits(ocean,[bsf_mvhv2],ocean.hvc_pct)
best_fullHV_syss,ocean.hv_circuits=opt_rollUp_firstInLine(ocean,ocean.hv_circuits)
best_fullMV_syss,ocean.mv_circuits=opt_rollUp_firstInLine(ocean,ocean.mv_circuits)
mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
best_fullMVHV_syss,mvhv_circuits=opt_rollUp_firstInLine(ocean,mvhv_circuits)
bsf_mvhv=combineAndrank_id(best_full_syssMV,best_full_syssHV,best_full_syss)
best_full_syss=keep_unic(best_full_syss)

for i=1:length(ocean.mv_circuits)
    if (ocean.mv_circuits[i][1].id != ocean11.mv_circuits[i][1].id)
        println("Full: "*string(ocean.mv_circuits[i][1].id)*" partial: "*string(ocean11.mv_circuits[i][1].id))
    end
end
ocean=load("Zero_size/tempfiles/data/solutions/clean_ocean.jld2")
bsf_mvhv=load("Zero_size/tempFiles/data/circuits/circ_bsf_mvhv2.jld2")["circuits"]
#25.499868+0.12514786+41.13849

circuits_old=deepcopy(circuits)
ppf_testing(bsf_mvhv)
ppf_testing(ocean.mv_circuits[14])
ppf_printIt(ocean,fr[1])

ppf_cbl_count(ocean.mv_circuits)
ppf_printCost(bsf_mvhv)
partial_sets_15=deepcopy(Q_2bd)
ppf_equipment_OSS_MOG(ocean,bsf_mvhv[1])
ppf_equipment_OSS_MOG(ocean,fr[1])
#615.003-mv, 615.3868-hv,615.1236 -mhv
mvc=ocean.mv_circuits
ocean
ocean=load("Zero_size/tempFiles/data/solutions/clean_ocean.jld2")["ocean"]
ocean.mv_circuits[252][14]=opt_circOpt(ocean.mv_circuits[252][14],ocean)
longest_cable=lof_pnt2pnt_dist(ocean.owpps[length(ocean.owpps)].node.xy,ocean.pccs[length(ocean.pccs)].node.xy)
new_coords=opt_reAdjust_oss(ocean.mv_circuits[252][14],ocean.owpps[1].mv_zone,ocean.sys.mvCl,10e-6)
ocean.mv_circuits[252][14]=opt_reAdjust_cbls(ocean.mv_circuits[252][14],new_coords,ocean,longest_cable)

ocean.mv_circuits[3]
bsf_mvhv[9]


for (i,qcs) in enumerate(ocean.mv_circuits[252])
    println(i)
    for (i,qc) in enumerate(ocean.mv_circuits[252])

        println(string(i)*" - "*string(length(qc.osss_owp)+length(qc.osss_mog)))
    end
    println()
end
