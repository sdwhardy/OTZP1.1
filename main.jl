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
include("eqp/functions.jl")#
include("post_process/functions.jl")#costs
include("astar/functions.jl")
include("optimization/functions.jl")
include("optimization/elipse_functions.jl")
include("topology/functions.jl")



function main()

    @time ocean=lof_layoutEez_basis()
    @time ocean.owpps=lof_order2Pcc(ocean,ocean.pccs[2])
    @time ocean=lof_layoutEez_expand(ocean,ocean.pccs[2])
    @time ocean.circuits=opt_hvOSSplacement(ocean,ocean.pccs[2])
    ppf_testing(ocean)
    ppf_equipment(ocean,ocean.circuits[127])
    #ppf_saveSystem(ocean,side*"HV")
    @time ocean.circuits=opt_mvOSSplacement(ocean,ocean.owpps,ocean.pccs[2])
    ppf_testing(ocean)
    #ppf_saveSystem(ocean,side*"MV")
    @time opt_compoundOSS(ocean)
    ppf_testing(ocean)
    #ppf_saveSystem(ocean,"Combo_8owpps_boo")
    @time best_full_syss,ocean.circuits=opt_rollUp(ocean)
    #ocean=lof_layoutEez_expand_testing(ocean,ocean.pccs[2])
    #ppf_layout_testing(ocean)

    ocean=load("tempFiles/data/solutions/Combo_8owpps.jld2")["ocean"]
    ocn_l=load("tempFiles/data/solutions/lowCB.jld2")["ocean"]
    ocn_h=load("tempFiles/data/solutions/highCB.jld2")["ocean"]

    ppf_printOcnXY_cables(ocn_l,ocn_l.circuits[51])
    printLines(ocean)
    gr()
    gui()

    ppf_equipment(ocean,ocean.circuits[20])
    ppf_equipment(ocean,best_full_syss[1])
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
