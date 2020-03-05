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
#    ppf_saveSystem(ocean,"zero_ocean001")
#    ocean=load("Zero_size/tempfiles/data/solutions/zero_ocean001.jld2")["ocean"]
    @time opt_compoundOSS(ocean)


    best_fullHV_syss,ocean.hv_circuits=opt_rollUp(ocean,ocean.hv_circuits)
    best_fullMV_syss,ocean.mv_circuits=opt_rollUp(ocean,ocean.mv_circuits)
    mvhv_circuits=optSysCompare(ocean.hv_circuits,ocean.mv_circuits)
    best_fullMVHV_syss,ocean.mv_circuits=opt_rollUp(ocean,mvhv_circuits)
    bsf_mvhv=combineAndrank(best_fullMV_syss,best_fullHV_syss,best_fullMVHV_syss)
    return ocean, bsf_mvhv
end
@time ocean, bsf_mvhv=main()
ppf_testing(ocean.hv_circuits)
ppf_testing(ocean.mv_circuits)
ppf_equipment_OSS_MOG(ocean,bsf_mvhv[1])
ppf_cbl_count(bsf_mvhv)

for c in chc_circ
    println(c.cost)
end
