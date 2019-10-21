using JuliaInterpreter
using Debugger
using DataFrames,XLSX,CSV
using StatsPlots, SpecialFunctions
using Polynomials


include("wind/struct.jl")
include("cost/struct.jl")#
include("eqp/struct.jl")#
include("layout/struct.jl")#

include("cost/data.jl")#
include("eqp/data.jl")#

include("eens/functions.jl")
include("wind/functions.jl")
include("cost/functions.jl")#
include("layout/functions.jl")#
include("eqp/functions.jl")#
include("post_process/functions.jl")#costs
include("astar/functions.jl")
function main()
    ocean=lof_layoutEez()
    start=ocean.owpps[3].node
    goal=ocean.pccs[1].node
    path=as_Astar(start,goal,ocean)
    #ppf_printOcnGPS(ocean)
    ppf_printOcnXY(ocean,path,start)

    ppf_printOcnXY(ocean)
end

main()


l=50
S=500
kv=66

cb=cstF_MvCbl(l,S,kv,ocean.owpps[1].wnds[1],ocean.finance)
kv=132
cbo=cstF_HvCblo2o(l,S,kv,ocean.owpps[1].wnds[length(ocean.owpps[1].wnds)],ocean.finance)
kv=220
cb1=cstF_HvCblo2p(l,S,kv,ocean.owpps[1].wnds[1],ocean.finance)
kv=400
cb2=cstF_HvCblo2o(l,S,kv,ocean.owpps[1].wnds[1],ocean.finance)
cstF_xfo_oss(0,ocean.owpps[1].wnds[1],ocean.finance)
cstF_xfo_pcc(500,ocean.owpps[1].wnds[1],ocean.finance)
