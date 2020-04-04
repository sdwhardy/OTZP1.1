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
    #ocean.mv_circuits=opt_mvOSSplacement(ocean,ocean.pccs[2])
    #ocean.mv_circuits=opt_updateMVocean(ocean)
    #ocean.mv_circuits,ocean.hv_circuits=check_ids(ocean.mv_circuits,ocean.hv_circuits)
    ocean.hv_circuits = opt_readjust_circuits(ocean,ocean.hv_circuits)
    #ocean.mv_circuits = opt_readjust_circuits(ocean,ocean.mv_circuits)
    #@time ocean=opt_compoundOSS(ocean)
    #ppf_saveSystem(ocean,"clean_ocean")
    best_fullHV_syss,ocean.hv_circuits=opt_rollUp(ocean,ocean.hv_circuits)
    #best_fullMV_syss,ocean.mv_circuits=opt_rollUp(ocean,ocean.mv_circuits)
    #mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
    #best_fullMVHV_syss,mvhv_circuits=opt_rollUp(ocean,mvhv_circuits)
    #bsf_mvhv=combineAndrank(best_fullMV_syss,best_fullHV_syss,best_fullMVHV_syss)
    #ppf_saveCircuit(bsf_mvhv,"bsf_mvhv")
    return ocean, best_fullHV_syss
end
#ocean11=ocean
@time ocean, bsf_mvhv2=main()

best_fullHV_syss,ocean.hv_circuits=opt_rollUp_firstInLine(ocean,ocean.hv_circuits)
best_fullMV_syss,ocean.mv_circuits=opt_rollUp_firstInLine(ocean,ocean.mv_circuits)
mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
best_fullMVHV_syss,mvhv_circuits=opt_rollUp_firstInLine(ocean,mvhv_circuits)
bsf_mvhv=combineAndrank_id(best_fullMV_syss,best_fullHV_syss,best_fullMVHV_syss)
bsf_mvhv=keep_unic(bsf_mvhv)

for i=1:length(ocean.mv_circuits)
    if (ocean.mv_circuits[i][1].id != ocean11.mv_circuits[i][1].id)
        println("Full: "*string(ocean.mv_circuits[i][1].id)*" partial: "*string(ocean11.mv_circuits[i][1].id))
    end
end

bsf_mvhv=load("Zero_size/tempFiles/data/circuits/circ_bsf_mvhv.jld2")["circuits"]
#25.499868+0.12514786+41.13849


ppf_testing(bsf_mvhv)
ppf_testing(ocean.mv_circuits)
ppf_printIt(ocean,bsf_mvhv)
ppf_cbl_count(bsf_mvhv)
ppf_printCost(bsf_mvhv)
ppf_equipment_OSS_MOG(ocean,bsf_mvhv[3575])
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
