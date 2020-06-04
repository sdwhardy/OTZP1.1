include("economics/main.jl")
include("layout/main.jl")
include("post_process/main.jl")
include("greedy/main.jl")
#limit to 2000MW single nodes
function main(pccs,owpps)
    pcc,owpps,offset=utm_gps2xy(pccs[2],owpps,32,true)#base belgian case do nothing
    #owpps=test_case001(owpps)#hvdc 1
    #owpps=test_case002(owpps)#hvdc 2
    #owpps=test_case003(owpps)#mid-point
    #owpps=test_case004(owpps)#all positive belgian case
    @time ocean=ocean_layout(pcc,owpps,offset)#takes xy coordinates to allow easy testing directly in cartesian
    return ocean
end
pccs=get_PccData()
owpps=get_OwppData()
@time ocn=main(pccs,owpps)
save("v2.0/tempfiles/ocean/greedy_tnep_ocean.jld2", "ocean",ocn)
#testing
#owpps[7].node.xy.x=abs(owpps[7].node.xy.x) - for testing purposes

plot_ocean(ocn)
old_meth=load("v2.0/tempfiles/ocean/b4Opt_tnep_ocean.jld2")["ocean"]


for circ in ocn.mv_circuits
    println(string(circ[1].decimal)*" - "*string(circ[1].cost))
end
println()
for circ in old_meth.hv_circuits
    println(string(circ[1].decimal)*" - "*string(circ[1].cost))
end

i=231
plot_circuit(ocn.mv_circuits[i][1])

plot_circuit(old_meth.mv_circuits[i][1])

ocn400=deepcopy(ocn)
owpp2oss_connectionMVAC=owpp_mvac_to_oss(ocn.hv_circuits[255][1].owpps[1],15.69,ocn.database,ks)

for (index,o220) in enumerate(ocn.hv_circuits)
    if (length(o220[1].owpps)>4)
        println(index)
    end
end


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


function test_case001(owpps)#hvdc
    owpps[1].node.xy.y=275
    owpps[1].node.xy.x=-30
    owpps[1].mva=250
    owpps[2].node.xy.y=185
    owpps[2].node.xy.x=30
    owpps[2].mva=250
    owpps[3].node.xy.y=250
    owpps[3].node.xy.x=50
    owpps[3].mva=2000
    owpps[4].node.xy.y=175
    owpps[4].node.xy.x=-50
    owpps[4].mva=250
    owpps[5].node.xy.y=200
    owpps[5].node.xy.x=25
    owpps[5].mva=250
    owpps[6].node.xy.y=225
    owpps[6].node.xy.x=-25
    owpps[6].mva=250
    owpps[7].node.xy.y=200
    owpps[7].node.xy.x=-50
    owpps[7].mva=250
    owpps[8].node.xy.y=235
    owpps[8].node.xy.x=50
    owpps[8].mva=250
    return owpps
end

function test_case002(owpps)#hvdc
    owpps[1].node.xy.y=275
    owpps[1].node.xy.x=30
    owpps[1].mva=2000
    owpps[2].node.xy.y=185
    owpps[2].node.xy.x=30
    owpps[2].mva=2000
    owpps[3].node.xy.y=250
    owpps[3].node.xy.x=50
    owpps[3].mva=2000
    owpps[4].node.xy.y=175
    owpps[4].node.xy.x=50
    owpps[4].mva=1000
    owpps[5].node.xy.y=200
    owpps[5].node.xy.x=25
    owpps[5].mva=1000
    owpps[6].node.xy.y=225
    owpps[6].node.xy.x=25
    owpps[6].mva=500
    owpps[7].node.xy.y=200
    owpps[7].node.xy.x=50
    owpps[7].mva=250
    owpps[8].node.xy.y=235
    owpps[8].node.xy.x=50
    owpps[8].mva=600
    return owpps
end

function test_case003(owpps)#mid point compensation
    owpps[1].node.xy.y=185
    owpps[1].node.xy.x=60
    owpps[1].mva=600
    owpps[2].node.xy.y=165
    owpps[2].node.xy.x=40
    owpps[2].mva=400
    owpps[3].node.xy.y=155
    owpps[3].node.xy.x=55
    owpps[3].mva=1250
    owpps[4].node.xy.y=175
    owpps[4].node.xy.x=45
    owpps[4].mva=800
    owpps[5].node.xy.y=200
    owpps[5].node.xy.x=25
    owpps[5].mva=300
    owpps[6].node.xy.y=125
    owpps[6].node.xy.x=35
    owpps[6].mva=500
    owpps[7].node.xy.y=135
    owpps[7].node.xy.x=50
    owpps[7].mva=1000
    owpps[8].node.xy.y=145
    owpps[8].node.xy.x=30
    owpps[8].mva=750
    return owpps
end


function test_case004(owpps)
    owpps[1].node.xy.x=abs(owpps[1].node.xy.x)
    owpps[2].node.xy.x=abs(owpps[2].node.xy.x)
    owpps[3].node.xy.x=abs(owpps[3].node.xy.x)
    owpps[4].node.xy.x=abs(owpps[4].node.xy.x)
    owpps[5].node.xy.x=abs(owpps[5].node.xy.x)
    owpps[6].node.xy.x=abs(owpps[6].node.xy.x)
    owpps[7].node.xy.x=abs(owpps[7].node.xy.x)
    owpps[8].node.xy.x=abs(owpps[8].node.xy.x)
    return owpps
end
