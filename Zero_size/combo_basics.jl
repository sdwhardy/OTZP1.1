using Combinatorics, Ipopt, JuMP, Cbc
using Geodesy
N = 3
combis=reverse.(Iterators.product(fill((0,1,2),N)...))[:]
combi_pair=Array{Tuple{NTuple{N,Int64},Tuple{Int64,Int64,Int64}},1}()
#xfo_types=Array{Tuple{Int64,Int64},1}()
for comb in combis
    zeroz=count(x->x==0,comb)
    onez=count(x->x==1,comb)
    twoz=count(x->x==2,comb)
    push!(combi_pair,(comb,(zeroz,onez,twoz)))
    #zero2two=sort([zeroz,onez])
    #push!(xfo_types,((zero2two[1],zero2two[2]),))
end

a=[66,220,400]
collect(permutations(a,8))

typeof(((10,5),69.69))
unique(xfo_types)
b=[0,0,0,0,0,0,0,0]
n=2
b=zeros(Int8,1,8)
B=Array{Tuple{Array{Int8,2},Int64},1}()

PWR_tbl[2]=1
#owpps_tbl=top_hvTopos(8)
ts=collect(combinations(a))
for st in ts
    b2=deepcopy(b)
    dec=0
    for bits in st
        b2[bits+1]=1
        dec=dec+2^(bits)
    end
    push!(B,(b2,deepcopy(dec)))
end
sort!(B, by = x -> x[2])
return B
t=collect(combinations(a,n))
print(t.f[1])

function topologies(n)
    total=0
    for i=2:n
        m1=binomial(n,i)
        for j=2:i
            total=total+m1*binomial(n-i,j)
        end
    end
    return total
end
t=topologies(8)*2^4


println(t)
owpps_tbl_hv=top_hvTopos(9)
for i=1:length(t)
    println(t[i])
    dc=0
    for bt in t[i]
        dc=dc+2^bt
    end
    println(dc)
    println(top_dec2bin(dc))
    println()
end
top_dec2bin(128)
top_bin2dec([0,0,0,0,0,1,1,1])
    owpps_tbl_hv=owpps_tbl_hv[end:-1:1,end:-1:1]
    oss_system=circuit()
    oss_systems=Array{circuit,1}()


    # latitude and longitude for two cities in Sweden
    stockholm_lla = LLA(59.3293, 18.0686)
    uppsala_lla = LLA(59.8586, 17.6389)
    pcc_lla=LLA(51.32694,3.183611)
    norther_lla=LLA(51.59,2.939972)
    2.939972	51.59

    # all of Sweden uses UTM zone 33 north
    utm_sweden = UTMfromLLA(32, true, wgs84)

pcc_utm = utm_sweden(pcc_lla)

norther_utm = utm_sweden(norther_lla)

diff = pcc_utm.x - norther_utm.x, pcc_utm.y - norther_utm.y

√sum(diff.^2)/1000



minmax_range(ocean.hv_circuits)
cst_mvakm=(10+ocean.mv_circuits[3][1].osss_mog[1].xfmrs[1].costs.ttl)*1000/500
function minmax_range(circs)
    cst_min=Inf
    cst_max=0
    cst_mean=0
    for c in circs[1:1]
        cst_mvakm=0
        cst_mvakm=c[1].pcc_cbls[1].costs.cpx*1000*1.25/(c[1].pcc_cbls[1].length*250*length(c[1].owpps))
        println(cst_mvakm)
        if (cst_mvakm<cst_min)
            cst_min=deepcopy(cst_mvakm)
        end
        if (cst_mvakm>cst_max)
            cst_max=deepcopy(cst_mvakm)
        end
        cst_mean=cst_mean+cst_mvakm
    end
    println("min: "*string(cst_min))
    println("max: "*string(cst_max))
    println("mean: "*string(cst_mean/length(circs)))
end


#below functional MIP as example


for_cons=deepcopy(hv_connections)
aft_cons=deepcopy(owpp2oss_connections)
w_p=deepcopy(wnd_power)
topo=deepcopy(circ)

topo.oss=[]
topo.MVcbls=[]
topo.HVcbls=[]
topo.PCCcbls=[]
topo.pcc.xfmrs=[]
topo.mog=[]

ac_oss=ks.pac_f+ks.pac_f*ks.opx_pl*npv_years()
dc_oss=ks.pdc_h+ks.pdc_h*ks.opx_pl*npv_years()+ks.conv_d+ks.conv_d*ks.opx_co*npv_years()

in_nodes=build_input_nodes(aft_cons)
out_nodes=build_output_nodes(for_cons)
candidates=build_candidates(in_nodes,out_nodes,ocn.database,ks)

#400
println(in_nodes["220.0"][1]["cost"]+in_nodes["220.0"][2]["cost"]+out_nodes["220.0"]["cost"]+candidates["220.0220.0"][2]["cost"]+ac_oss)

#400 input
println(in_nodes["220.0"][1]["cost"]+in_nodes["400.0"][2]["cost"]+out_nodes["400.0"]["cost"]
+candidates["220.0400.0"][1]["cost"]+candidates["400.0400.0"][1]["cost"]+ac_oss)

#400 output
println(out_nodes["400.0"]["cost"]+candidates["220.0400.0"][3]["cost"]+ac_oss)

#220
println(in_nodes["220.0"][1]["cost"]+in_nodes["220.0"][2]["cost"]+in_nodes["220.0"][3]["cost"]
+in_nodes["220.0"][4]["cost"]+in_nodes["220.0"][5]["cost"]+out_nodes["220.0"]["cost"]+candidates["220.0220.0"][5]["cost"]+ac_oss)

sol["300.0"][1]["in_"*"300.0"*string(sol["300.0"][1]["mva"])]
old_ocean=load("v2.0/tempfiles/ocean/pre_solve_tnep_ocean.jld2")["ocean"]
plot_circuit(old_ocean.hv_circuits[255][1])

cbl1a=1
cbl1b=2
cbl2a=1
cbl2b=2
xb=1.5

@constraint(m,vars["1"][1]+vars["1"][2]<=1)
@constraint(m,vars["2"][1]+vars["2"][2]<=1)

@objective(m, Min,cbl1a*vars["1"][1]+cbl1b*vars["1"][2]+cbl2a*vars["2"][1]+cbl2b*vars["2"][2]+xb*(1-vars["1"][2]))
optimize!(m)

x0=JuMP.value.((vars["2"][1]))

###########################
#transformer
xfo=deepcopy(database["transformers"][string(500.0)]["offshore"])
#converter
conv=deepcopy(database["converters"][string(500.0)]["offshore"])
#platforms
plat_ac=platform()
plat_ac.mva=500.0
plat_ac.wnd=xfo.wnd
plat_dc=deepcopy(plat_ac)
plat_ac=cost_ac_platform(plat_ac,ks)
plat_dc=cost_dc_platform(plat_dc,ks)
plat_dc=adjust_base_dc_platform(plat_dc,ks)
println(plat_dc.costs.ttl+plat_ac.costs.ttl+xfo.costs.ttl+conv.costs.ttl)

plat_ac=adjust_base_ac_platform(plat_ac,ks)
println(plat_ac.costs.ttl+xfo.costs.ttl)


##############################
