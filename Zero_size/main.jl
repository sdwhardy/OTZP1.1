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
    ocean.mv_circuits=opt_updateMVocean(ocean)
    ocean.hv_circuits = opt_readjust_circuits(ocean,ocean.hv_circuits)
    ocean.mv_circuits = opt_readjust_circuits(ocean,ocean.mv_circuits)
    ppf_saveSystem(ocean,"baseMVHV_ocean_lm")
    #ppf_saveSystem(ocean,"zero_ocean001")
    #ocean=load("Zero_size/tempfiles/data/solutions/zero_ocean001.jld2")["ocean"]
    @time opt_compoundOSS(ocean)

    ppf_saveSystem(ocean,"compMVHV_ocean_lm")
    best_fullHV_syss,ocean.hv_circuits=opt_rollUp(ocean,ocean.hv_circuits)
    best_fullMV_syss,ocean.mv_circuits=opt_rollUp(ocean,ocean.mv_circuits)
    ppf_saveSystem(ocean,"roldUpMVHV_ocean_lm")
    ppf_saveCircuit(best_fullMV_syss,"sol_MV_lm")
    ppf_saveCircuit(best_fullHV_syss,"sol_HV_lm")
    mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
    best_fullMVHV_syss,ocean.mv_circuits=opt_rollUp(ocean,mvhv_circuits)
    bsf_mvhv=combineAndrank(best_fullMV_syss,best_fullHV_syss,best_fullMVHV_syss)
    ppf_saveCircuit(bsf_mvhv,"sol_MVHV_lm")
    return ocean, bsf_mvhv
end
function main_testing()
    ocean=load("Zero_size/tempfiles/data/solutions/opt_ocean001.jld2")["ocean"]
    @time opt_compoundOSS(ocean)

    #ppf_saveSystem(ocean,"zero_oceanXLarge1")
    best_fullHV_syss,ocean.hv_circuits=opt_rollUp(ocean,ocean.hv_circuits)
    best_fullMV_syss,ocean.mv_circuits=opt_rollUp(ocean,ocean.mv_circuits)
    mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
    best_fullMVHV_syss,ocean.mv_circuits=opt_rollUp(ocean,mvhv_circuits)
    bsf_mvhv=combineAndrank(best_fullMV_syss,best_fullHV_syss,best_fullMVHV_syss)
    return ocean, bsf_mvhv
end
@time ocean, bsf_mvhv=main_testing()
@time ocean, bsf_mvhv=main()
ppf_saveSystem(ocean,"zero_oceanMed")
ppf_saveCircuit(bsf_mvhv,"Med")
ppf_testing(bsf_mvhv)
ppf_testing(ocean.mv_circuits)
ppf_equipment_OSS_MOG(ocean,bsfm[1])
gui()
ppf_cbl_count(bsf)
function cst_total(circ)
    cst=0
    for (i,c) in enumerate(circ)
        if (i==1 || i==2 || i==4 || i==8 || i==16 || i==32 || i==64 || i==128)
            cst=cst+c.cost
        end
        println(cst)
    end

end
cst_total(ocean.hv_circuits)

ocean=load("Zero_size/tempfiles/data/solutions/roldUpMVHV_ocean.jld2")["ocean"]
bsfh=load("Zero_size/tempfiles/data/circuits/circ_sol_HV.jld2")["circuits"]
bsfm=load("Zero_size/tempfiles/data/circuits/circ_sol_MV_no.jld2")["circuits"]
bsf=load("Zero_size/tempfiles/data/circuits/circ_XLarge.jld2")["circuits"]
ppf_cbl_count(bsfm)
for (i,c) in enumerate(bsfm)
    println(string(i)*" - "*string(c.cost))
end
