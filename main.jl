#using JuliaInterpreter
#using Debugger
using DataFrames,XLSX,CSV
using StatsPlots, SpecialFunctions
using Polynomials
using JuMP, Ipopt


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
include("topology/functions.jl")



function main()
    @time ocean=lof_layoutEez()
    @time ocean.circuits=opt_mvOSSplacement(ocean,ocean.owpps,ocean.pccs[2])
    #Working here returning the wrong equipment I think (drunk desu) choosing hv cables when mv are better
    test[1]=opt_hvOSSplacement(ocean,ocean.pccs[2])
    ppf_printOcnXY(ocean,ocean.circuits[1].pths)

    #start=ocean.owpps[3].node
    #goal=ocean.pccs[1].node
    #path=as_Astar(start,goal,ocean.discretedom.nodes)
    #ppf_printOcnGPS(ocean)
    #ppf_printOcnXY(ocean,path,start)

    #ppf_printOcnXY(ocean)
end
function testing(ocn,tst)
    mv=0
    hv=0
    for i=1:length(ocn.circuits)
        if (ocn.circuits[i].cost<tst[i].cost)
            mv=mv+1
            println(string(i)*" MV: "*string(ocn.circuits[i].binary)*" - "*string(ocn.circuits[i].cost))
        else
            #println("HV: "*string(ocn.circuits[i].binary)*" - "*string(tst[i].cost))
            hv=hv+1
        end
    end
    println(mv)
    println(hv)
end
testing(ocean,test)
main()


l=5
km=l
S=250
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
cstF_xfo_pcc(S,ocean.owpps[1].wnd,ocean.finance)
