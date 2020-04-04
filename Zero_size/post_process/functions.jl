function ppf_saveSystem(ocn,nme)
	if (isdir("Zero_size/tempfiles/data/solutions/") == false)
		mkdir("Zero_size/tempfiles/data/solutions")
    else
    end

    #save cicuits
    save("Zero_size/tempfiles/data/solutions/"*nme*".jld2", "ocean",ocn)
end

function ppf_printIt(ocean,bsf_mvhv)
    for i=1:12
        ppf_equipment_OSS_MOG(ocean,bsf_mvhv[i])
    end
	#gui()
end

function ppf_printCost(bsf_mvhv)
	for (i,c) in enumerate(bsf_mvhv)
	    println(string(i)*" - "*string(c.cost))
	end
end

function ppf_saveCircuit(oss_systems,nme)

    #check appropriate directories exist
	if (isdir("Zero_size/tempfiles/data/circuits/") == false)
		mkdir("Zero_size/tempfiles/data/circuits")
    else
    end

    #save cicuits
    save("Zero_size/tempfiles/data/circuits/circ_"*nme*".jld2", "circuits",oss_systems)
#	load("tempFiles/data/circuits/CRC_"*string(1)*string(2)*".jld2")["circuits"]
end

function ppf_testing(circ)
    total=0
    for i=1:length(circ)
        println(string(circ[i].decimal)*") Cst:"*string(circ[i].cost)*", mvC:"*string(length(circ[i].owp_MVcbls))*", hvC:"*string(length(circ[i].owp_HVcbls))*", oss:"*string(length(circ[i].osss_owp))*", mog:"*string(length(circ[i].osss_mog))*", o2oC:"*string(length(circ[i].oss2oss_cbls)))
        total=total+circ[i].cost
    end
    println(total)
end

function ppf_equipment_OSS_MOG(ocn,circ)
    println("MV cables: ")
    for mvc in circ.owp_MVcbls
        println(string(mvc.mva*mvc.num)*" MVA, "*string(mvc.elec.volt)*" KV, "*string(mvc.length)*" KM, "*string(mvc.costs.ttl)*" EURO, "*string(mvc.pth[1].num)*" - "*string(mvc.pth[length(mvc.pth)].num))
    end
    println()
    println("OSS: ")
    for mvc in circ.osss_owp
        xfmTTL=0
        println(string(mvc.node.num)*": ")
        for x in mvc.xfmrs
            xfmTTL=(xfmTTL+(x.costs.ttl))
            print(string(x.num)*" - "*string(x.mva)*" MVA, "*string(x.lv)*" KV, "*string(x.hv)*" KV, ")
        end
        println(string(mvc.base_cost+xfmTTL)*" EURO")
    end
    println()
    println()
    println("HV cables: ")
    for mvc in circ.owp_HVcbls
        println(string(mvc.mva*mvc.num)*" MVA, "*string(mvc.elec.volt)*" KV, "*string(mvc.length)*" KM, "*string(mvc.costs.ttl)*" EURO, "*string(mvc.pth[1].num)*" - "*string(mvc.pth[length(mvc.pth)].num))
    end
    println()
    println("MOG: ")
    for mvc in circ.osss_mog
        xfmTTL=0
        println(string(mvc.node.num)*": ")
        for x in mvc.xfmrs
            xfmTTL=(xfmTTL+(x.costs.ttl))
            println(string(x.num)*" - "*string(x.mva)*" MVA, "*string(x.lv)*" KV, "*string(x.hv)*" KV, ")
        end
        println(string(mvc.base_cost+xfmTTL)*" EURO")
    end
    println()
    println()
    println("O2O cables: ")
    for mvc in circ.oss2oss_cbls
        println(string(mvc.mva*mvc.num)*" MVA, "*string(mvc.elec.volt)*" KV, "*string(mvc.length)*" KM, "*string(mvc.costs.ttl)*" EURO, "*string(mvc.pth[1].num)*" - "*string(mvc.pth[length(mvc.pth)].num))
    end
    println()
    println("PCC cables: ")
    for mvc in circ.pcc_cbls
        println(string(mvc.mva*mvc.num)*" MVA, "*string(mvc.elec.volt)*" KV, "*string(mvc.length)*" KM, "*string(mvc.costs.ttl)*" EURO, "*string(mvc.pth[1].num)*" - "*string(mvc.pth[length(mvc.pth)].num))
    end
    println()
    println("PCC: ")
    xfmTTL=0
    println(string(circ.pcc.node.num)*": ")
    for x in circ.pcc.xfmrs
        xfmTTL=(xfmTTL+(x.costs.ttl))
        println(string(x.num)*" - "*string(x.mva)*" MVA, "*string(x.lv)*" KV, "*string(x.hv)*" KV, ")
    end
    println(string(circ.pcc.base_cost+xfmTTL)*" EURO")
    ppf_printOcnXY_cables_OSS_MOG(ocn,circ)
end

function ppf_printOcnXY_cables_OSS_MOG(ocn,pths)

	plotly()
	p=plot([ocn.pccs[1].node.xy.x],[ocn.pccs[1].node.xy.y],color = :green,markersize=2,seriestype=:scatter,xticks = 0:2:30,xlims=(0,30),yticks = 0:5:56,label="",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	for indx = 2:length(ocn.pccs)
		plot!(p,[ocn.pccs[indx].node.xy.x],[ocn.pccs[indx].node.xy.y],color = :green,seriestype=:scatter,label="")
	end
	xs=Float32[]
	ys=Float32[]
	tx=Tuple[]


	for (i,owpp) in enumerate(ocn.owpps)
		txt=text(string(i),15,:blue,:left)
		push!(tx,(owpp.node.xy.x,owpp.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :blue,seriestype=:scatter,label="")


	for cbl in pths.owp_MVcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
				push!(xs,pth.xy.x)
				push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :blue,linewidth = 2,label="")
	end
	for cbl in pths.owp_HVcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,linewidth = 2,label="")
	end

	for cbl in pths.oss2oss_cbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,linewidth = 2,label="")
	end
	for cbl in pths.pcc_cbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,linewidth = 2,label="")
	end
	xs=Float32[]
	ys=Float32[]
	tx=Tuple[]

	for (i,oss) in enumerate(pths.osss_owp)
		txt=text(string(""),15,:red,:right)
		push!(tx,(oss.node.xy.x,oss.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :black,seriestype=:scatter,label="")
	xs=Float32[]
	ys=Float32[]
	tx=Tuple[]
	for (i,oss) in enumerate(pths.osss_mog)
		txt=text(string(""),15,:black,:right)
		push!(tx,(oss.node.xy.x,oss.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :black,seriestype=:scatter,label="")
	p
	gui()
end

function ppf_cbl_count(ocnhv)
    cbles=0
    oss=0
    for (i,circ) in enumerate(ocnhv)

        cbles=length(circ.pcc_cbls)
		cbleso20=length(circ.oss2oss_cbls)
        oss=length(circ.osss_mog)+length(circ.osss_owp)
		#if (oss>7 && cbles>1 && cbleso20==0)
        	println(string(i)*" OSS: "*string(oss)*" PCC Cables: "*string(cbles)*" Cost: "*string(circ.cost)*" ID: "*string(circ.id))
		#end

        oss=0
		cbles=0
    end
end
