include("economics/main.jl")
include("layout/main.jl")
include("post_process/main.jl")
include("greedy/main.jl")
#limit to 2000MW single nodes
function main(pccs,owpps)
    pccs=get_PccData()
    owpps=get_OwppData()
    pcc,owpps,offset=utm_gps2xy(pccs[2],owpps,32,true)#base belgian case do nothing
    #owpps=test_case001(owpps)#hvdc 1
    #owpps=test_case002(owpps)#hvdc 2
    #owpps=test_case003(owpps)#mid-point
    #owpps=test_case004(owpps)#all positive belgian case
    #owpps=test_case005(owpps)
    @time ocean=ocean_layout(pcc,owpps,offset)#takes xy coordinates to allow easy testing directly in cartesian
    @time TH=greedy_search(ocean)
    @time FR=FR_from_TH(deepcopy(TH),ones(Int8,1,length(ocean.owpps)),ocean)
    return FR,TH
end
pccs=get_PccData()
owpps=get_OwppData()
@time fr,th=main(pccs,owpps)
(ones(Int8,1,length(owpps))==test)
TH1=deepcopy(TH)
th=TH_reduction(TH1,length(owpps))
fr=fr_from_th(th,ocean)

ocean=deepcopy(TH1)
save("v2.0/tempfiles/ocean/TH_220.jld2", "TH",TH)
save("v2.0/tempfiles/ocean/FR_set_forward_cable.jld2", "fr",fr)
#testing
#owpps[7].node.xy.x=abs(owpps[7].node.xy.x) - for testing purposes
for cbls in ocean.database["cables"]["220.0"]["300.0"]
    #for cbl in cbls
        println(cbls.costs.ttl/(cbls.mva*cbls.length)*1000)
    #end
end
plot_ocean(ocean)
ocean=load("v2.0/tempfiles/ocean/ocean16.jld2")["ocean"]
TH[10]
cable_220=deepcopy(hvac_cable(500.0,24.928,FR[2].PCCcbls[3].wnd,ocean.database["cables"]["220.0"][string(500.0)],ks))
for (i,circ) in enumerate(fr)
    if (length(circ.O2Ocbls)>0)
        println(i)
    end
    #println(string(i)*" "*string(circ.cost))
end
println()
for circ in old_meth.hv_circuits
    println(string(circ[1].decimal)*" - "*string(circ[1].cost))
end

i=231
239.87378+5.2676396+53.50167
fr=best_tHs_from_TH_parts(r)
FR220=deepcopy(fr[24])
plotly()
p=plot()
plot_circuit(p,FR[1])
gui()
for (i,circ) in enumerate(r)
    if (circ[1].cost<ocean.mv_circuits[i][1].cost)
        println(i)
    end
end
for (i,fr) in enumerate(fr)
    if (length(fr.O2Ocbls)>0)
        println(i)
    end
end
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
ocean.database["cables"]
ks=get_Cost_Data()
mva=1000.0
for mva in [250.0,500.0,750.0,1000.0,1500.0,2000.0,4000.0]
    plat_base=platform()
    plat_base.acdc="ac"
    plat_base.mva=mva
    plat_base.wnd=ocean.owpps[3].wnd
    plat_base=cost_ac_platform(plat_base,ks)
    plat_base=adjust_base_ac_platform(plat_base,ks)
    println(plat_base.costs.ttl)
end

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
    owpps[1].mva=4000.0
    owpps[2].mva=500.0
    owpps[3].mva=500.0
    owpps[4].mva=500.0
    owpps[5].mva=500.0
    owpps[6].mva=500.0
    owpps[7].mva=500.0
    owpps[8].mva=500.0
    return owpps
end

function test_case005(owpps)
    owpps=owpps[1:4]
#=    owpps[5].node.xy.y=owpps[1].node.xy.y
    owpps[6].node.xy.y=owpps[2].node.xy.y
    owpps[7].node.xy.y=owpps[3].node.xy.y
    owpps[8].node.xy.y=owpps[4].node.xy.y
    owpps[5].node.xy.x=abs(owpps[1].node.xy.x-5)
    owpps[6].node.xy.x=abs(owpps[2].node.xy.x-5)
    owpps[7].node.xy.x=abs(owpps[3].node.xy.x-5)
    owpps[8].node.xy.x=abs(owpps[4].node.xy.x-5)=#
    for owp in owpps[1:4]
        owp_copy=deepcopy(owp)
        owp_copy.node.xy.x=abs(owp.node.xy.x)-15
        push!(owpps,owp_copy)
    end
    #=for owp in owpps[1:4]
        owp_copy=deepcopy(owp)
        owp_copy.node.xy.x=abs(owp.node.xy.x)-5
        push!(owpps,owp_copy)
    end=#
    #=owpps[1].mva=500
    owpps[2].mva=500
    owpps[3].mva=500
    owpps[4].mva=500
    owpps[5].mva=500
    owpps[6].mva=500
    owpps[7].mva=500
    owpps[8].mva=500=#
    return owpps
end
