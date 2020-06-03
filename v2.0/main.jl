include("economics/main.jl")
include("layout/main.jl")
include("post_process/main.jl")

function main(pccs,owpps)
    pcc,owpps,offset=utm_gps2xy(pccs[2],owpps,32,true)
    owpps[1].node.xy.y=25
    owpps[1].node.xy.x=155
    owpps[1].mva=800
    owpps[2].node.xy.y=185
    owpps[2].node.xy.x=12
    owpps[2].mva=500
    owpps[3].node.xy.y=150
    owpps[3].node.xy.x=50
    owpps[3].mva=500
    owpps[4].node.xy.y=75
    owpps[4].node.xy.x=175
    owpps[4].mva=500
    owpps[5].node.xy.y=20
    owpps[5].node.xy.x=200
    owpps[5].mva=800
    owpps[6].node.xy.y=160
    owpps[6].node.xy.x=50
    owpps[6].mva=500
    owpps[7].node.xy.y=230
    owpps[7].node.xy.x=0
    owpps[7].mva=500
    owpps[8].node.xy.y=225
    owpps[8].node.xy.x=50
    owpps[8].mva=500

    #=owpps[1].node.xy.x=abs(owpps[1].node.xy.x)


    owpps[2].node.xy.x=abs(owpps[2].node.xy.x)


    owpps[3].node.xy.x=abs(owpps[3].node.xy.x)


    owpps[4].node.xy.x=abs(owpps[4].node.xy.x)

    owpps[5].node.xy.x=abs(owpps[5].node.xy.x)

    owpps[6].node.xy.x=abs(owpps[6].node.xy.x)

    owpps[7].node.xy.x=abs(owpps[7].node.xy.x)

    owpps[8].node.xy.x=abs(owpps[8].node.xy.x)=#

    @time ocean=ocean_layout(pcc,owpps,offset)#takes xy coordinates to allow easy testing directly in cartesian
    @time B=make_set_B(ocean.owpps,pcc)
    ocean=make_set_Tb(ocean.owpps,ocean.pcc,B,ocean)
    return ocean
end
pccs=get_PccData()
owpps=get_OwppData()
@time ocn=main(pccs,owpps)
@time ocn.hv_circuits=finalize_circuit_layout(ocn.hv_circuits,ocn)
@time ocn.mv_circuits=finalize_circuit_layout(ocn.mv_circuits,ocn)
new_coords=finalize_mog_location(ocn.mv_circuits[25][1],1e-2)

#testing
#owpps[7].node.xy.x=abs(owpps[7].node.xy.x) - for testing purposes
for circ in ocn.hv_circuits
    println(string(circ[1].decimal)*" - "*string(circ[1].cost))
end
save("v2.0/tempfiles/ocean/b4Opt_tnep_ocean.jld2", "ocean",ocn)
old_meth=load("v2.0/tempfiles/ocean/b4Opt_tnep_ocean.jld2")["ocean"]
i=26
plot_circuit(ocn.hv_circuits[i][1])

plot_circuit(old_meth.hv_circuits[i][1])

ocn400=deepcopy(ocn)
owpp2oss_connectionMVAC=owpp_mvac_to_oss(ocn.hv_circuits[255][1].owpps[1],15.69,ocn.database,ks)

for (index,o220) in enumerate(ocn.hv_circuits)
    if (ocn.hv_circuits[index][1].cost<=ocn400.hv_circuits[index][1].cost)
        println(index)
    end
end

plot_ocean(ocn)
ocn.owpps[1].mva=500
ocn.owpps[1].node.xy.x=0
ocn.owpps[1].node.xy.y=250

for owpp in owpps
    owpp.node.xy.y=abs(owpp.node.xy.y)
end
owpps[7].node.xy.y=-1*abs(owpps[7].node.xy.y)
function printit(B)
    A=range(1,stop=255)
    for a in A
        isthere=false
        for b in B
            if (b[2]==a)
                isthere=true
            end
        end
        if (isthere==false)
            println(a)
        end
    end
end
printit(B)
ocean.owpps[3].node.xy.x=abs(ocean.owpps[3].node.xy.x)
plot_PCCnOWPP(ocn)
