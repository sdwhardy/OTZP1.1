include("economics/main.jl")
include("layout/main.jl")
include("post_process/main.jl")

function main(pccs,owpps)
    pcc,owpps,offset=utm_gps2xy(pccs[2],owpps,32,true)
    owpps[1].node.xy.y=250
    owpps[1].mva=500
    ocean=ocean_layout(pcc,owpps,offset)#takes xy coordinates to allow easy testing directly in cartesian
    B=make_set_B(ocean.owpps,pcc)
    return ocean,B
end
pccs=get_PccData()
owpps=get_OwppData()
@time ocn,B=main(pccs,owpps)
Tb=make_set_Tb(ocn.owpps,ocn.pcc,B,ocn)
plot_circuit(Tb[1][1])
#NOTE!!!
#pretty good results approx 3 ME diff though on HVDC though only 0.3% overall
plot_ocean(Tb[1][1])
ocn.owpps[1].mva=500
ocn.owpps[1].node.xy.x=0
ocn.owpps[1].node.xy.y=240

function make_set_Tb(owpps,pcc,B,ocn)
    Tb_basis_circuits=make_basis_of_Tb(owpps,pcc,B,ocn)
    return Tb_basis_circuits
end

function make_basis_of_Tb(owpps,pcc,B,ocn)
    hmv_circuits=Vector{Array{circuit,1}}()
    A=range(1,stop=length(owpps))
    for b in B
        hmv_circuit=circuit()
        hmv_circuit.binary=b[1]
        hmv_circuit.decimal=b[2]
        first_owp=findfirst(x->x==1,b[1])[2]
        hmv_circuit.base=owpps[first_owp]
        attached_owpps=findall(x->x==1,b[1])
        net_mva=0.0
        net_wnd=Array{wind,1}()
        for owp in attached_owpps
            push!(net_wnd,owpps[owp[2]].wnd)
            net_mva=net_mva+owpps[owp[2]].mva
            push!(hmv_circuit.owpps,deepcopy(owpps[owp[2]]))
        end
        hmv_circuit.mva=net_mva
        hmv_circuit.wnd=find_netWind(net_wnd)
        hmv_circuit.pcc=deepcopy(pcc)
        if (length(attached_owpps)==1)
            hmv_circuit.id="mvhv"*string(hmv_circuit.decimal)
            hmv_circuit=optimal_owpp2pcc(hmv_circuit,ocn)
        end
        push!(hmv_circuits,[deepcopy(hmv_circuit)])
    end
    return hmv_circuits
end








#testing
#owpps[7].node.xy.x=abs(owpps[7].node.xy.x) - for testing purposes
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
