using Combinatorics
using Geodesy
a=[0,1,2,3,4,5,6,7]
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
