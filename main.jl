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
    @time ocean.circuits=opt_hvOSSplacement(ocean,ocean.pccs[2])
    @time ocean.circuits=opt_compoundOSS(ocean,ocean.pccs[2])
    #ppf_printOcnXY(ocean,ocean.circuits[31].pths)

    #start=ocean.owpps[3].node
    #goal=ocean.pccs[1].node
    #path=as_Astar(start,goal,ocean.discretedom.nodes)
    #ppf_printOcnGPS(ocean)
    #ppf_printOcnXY(ocean,path,start)

    #ppf_printOcnXY(ocean)
end
ocn=main()

function testing(ocn)
    mv=0
    hv=0
    total=0
    for i=1:length(ocn.circuits)
        println(string(ocn.circuits[i].decimal)*") Cst:"*string(ocn.circuits[i].cost)*", mvC:"*string(length(ocn.circuits[i].owp_MVcbls))*", hvC:"*string(length(ocn.circuits[i].owp_HVcbls))*", oss:"*string(length(ocn.circuits[i].osss_owp))*", mog:"*string(length(ocn.circuits[i].osss_mog))*", o2oC:"*string(length(ocn.circuits[i].oss2oss_cbls)))
        total=total+ocn.circuits[i].cost
        #if (ocn.circuits[i].cost>=tst[i].cost)
        #    mv=mv+1
            #println(string(i)*" MV: "*string(tst[i].cost)*" - "*string(bst_sys[i].cost))
        #else
            #println("HV: "*string(ocn.circuits[i].binary)*" - "*string(tst[i].cost))
        #    hv=hv+1
        #end
    end
    println(total)
end
sea=eez()
sea.circuits=deepcopy(ocn)
testing(sea)
ppf_printOcnXY_cables(ocean,sea.circuits[31])
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
cstF_xfo_pcc(S,ocean.circuits[1].osss_mog[1].wnd,ocean.finance)
cstF_xfo_pcc(power_sum,wind_sum,ocean.finance)
