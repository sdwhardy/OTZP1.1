using DataFrames,XLSX
using StatsPlots, SpecialFunctions

include("layout/struct.jl")#
include("layout/functions.jl")#
include("post_process/functions.jl")#costs

function main()
    ocean=lof_layoutEez()
    ppf_printOcnGPS(ocean)
    ppf_printOcnXY(ocean)
    println(ocean.pccs[length(ocean.pccs)])
    println(ocean.owpps[length(ocean.owpps)])
    print(ocean.id_count)
end

main()
