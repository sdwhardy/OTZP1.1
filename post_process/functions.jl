function printLines(ocn)
	ocn=ocean
	plotly()
	p=plot()
	for edg in ocn.discretedom.edges[2000:5000]
		x0=ocn.discretedom.nodes[edg.head].xy.x
		x1=ocn.discretedom.nodes[edg.tail].xy.x
		y0=ocn.discretedom.nodes[edg.head].xy.y
		y1=ocn.discretedom.nodes[edg.tail].xy.y
		plot!(p,[x0,x1],[y0,y1])
	end
	p
	gui()
end

function ppf_testing(ocn)
    mv=0
    hv=0
    total=0
    for i=1:length(ocn.circuits)
        println(string(ocn.circuits[i].decimal)*") Cst:"*string(ocn.circuits[i].cost)*", mvC:"*string(length(ocn.circuits[i].owp_MVcbls))*", hvC:"*string(length(ocn.circuits[i].owp_HVcbls))*", oss:"*string(length(ocn.circuits[i].osss_owp))*", mog:"*string(length(ocn.circuits[i].osss_mog))*", o2oC:"*string(length(ocn.circuits[i].oss2oss_cbls)))
        total=total+ocn.circuits[i].cost
    end
    println(total)
end


function ppf_equipment(ocn,circ)
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
    ppf_printOcnXY_cables(ocn,circ)
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


function ppf_layout_testing(ocn)


	plotly()
	p=plot([ocn.discretedom.nodes[1].xy.x],[ocn.discretedom.nodes[1].xy.y],color = :red,markersize=2,seriestype=:scatter,xticks = 0:2:30,xlims=(0,30),yticks = 0:5:56,label="",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	for indx = 2:length(ocn.discretedom.nodes)
		plot!(p,[ocn.discretedom.nodes[indx].xy.x],[ocn.discretedom.nodes[indx].xy.y],color = :red,markersize=1,seriestype=:scatter,label="")
	end

	for indx = 1:length(ocn.pccs)
		plot!(p,[ocn.pccs[indx].node.xy.x],[ocn.pccs[indx].node.xy.y],color = :green,seriestype=:scatter,label="")
	end

	for owpp in ocn.owpps
		plot!(p,[owpp.node.xy.x],[owpp.node.xy.y],color = :blue,seriestype=:scatter,label="")
	end
	p
	gui()

end

function ppf_printOcnXY_cables(ocn,pths)

	plotly()
	p=plot([ocn.pccs[1].node.xy.x],[ocn.pccs[1].node.xy.y],color = :green,markersize=2,seriestype=:scatter,xticks = 0:2:30,xlims=(0,30),yticks = 0:5:56,label="",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	for indx = 2:length(ocn.pccs)
		plot!(p,[ocn.pccs[indx].node.xy.x],[ocn.pccs[indx].node.xy.y],color = :green,seriestype=:scatter,label="")
	end
	for owpp in ocn.owpps
		plot!(p,[owpp.node.xy.x],[owpp.node.xy.y],color = :blue,seriestype=:scatter,label="")
		for nd in owpp.zone.nodes
			plot!(p,[nd.xy.x],[nd.xy.y],color = :red,markersize=2,seriestype=:scatter,label="")
		end
	end
	for ng in ocn.nogos
		for nd in ng.nodes
			plot!(p,[nd.xy.x],[nd.xy.y],color = :red,markersize=2,seriestype=:scatter,label="")
		end
	end



	for cbl in pths.owp_MVcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
				push!(xs,pth.xy.x)
				push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :green,label="")
	end
	for cbl in pths.owp_HVcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :red,label="")
	end

	for cbl in pths.oss2oss_cbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :blue,label="")
	end
	for cbl in pths.pcc_cbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.pth
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,label="")
	end
	p
	gui()
end


function ppf_printOcnXY_cables_OSS_MOG(ocn,pths)

	plotly()
	p=plot([ocn.pccs[1].node.xy.x],[ocn.pccs[1].node.xy.y],color = :green,markersize=2,seriestype=:scatter,xticks = 0:2:30,xlims=(0,30),yticks = 0:5:56,label="",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	for indx = 2:length(ocn.pccs)
		plot!(p,[ocn.pccs[indx].node.xy.x],[ocn.pccs[indx].node.xy.y],color = :green,seriestype=:scatter,label="")
	end
	for owpp in ocn.owpps
		plot!(p,[owpp.node.xy.x],[owpp.node.xy.y],color = :blue,seriestype=:scatter,label="")
		xs=Float64[]
		ys=Float64[]

		push!(xs,owpp.node.xy.x-owpp.zone.neg_width)
		push!(ys,owpp.node.xy.y+owpp.zone.pos_height)

		push!(xs,owpp.node.xy.x+owpp.zone.pos_width)
		push!(ys,owpp.node.xy.y+owpp.zone.pos_height)

		push!(xs,owpp.node.xy.x+owpp.zone.pos_width)
		push!(ys,owpp.node.xy.y-owpp.zone.neg_height)

		push!(xs,owpp.node.xy.x-owpp.zone.neg_width)
		push!(ys,owpp.node.xy.y-owpp.zone.neg_height)

		push!(xs,owpp.node.xy.x-owpp.zone.neg_width)
		push!(ys,owpp.node.xy.y+owpp.zone.pos_height)

		plot!(p,xs,ys,color = :red,linewidth = 1,linestyle=:dot,label="")
	end
	for ng in ocn.nogos
			xs=Float64[]
			ys=Float64[]
			for pnts in ng.bndryPnts
				push!(xs,pnts.xy.x)
				push!(ys,pnts.xy.y)
			end
			push!(xs,ng.bndryPnts[1].xy.x)
			push!(ys,ng.bndryPnts[1].xy.y)
			plot!(p,xs,ys,color = :red,linewidth = 1,linestyle=:dot,label="")
	end



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
		txt=text(string(i),15,:red,:right)
		push!(tx,(oss.node.xy.x,oss.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :black,seriestype=:scatter,label="")
	xs=Float32[]
	ys=Float32[]
	tx=Tuple[]
	for (i,oss) in enumerate(pths.osss_mog)
		txt=text(string(i+1),15,:black,:right)
		push!(tx,(oss.node.xy.x,oss.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :black,seriestype=:scatter,label="")
	p
	gui()
end

function ppf_printOcnGPS(ocean)
	pcc=Array{Tuple,1}()
	gen=Array{Tuple,1}()
	bnd=Array{Tuple,1}()
	nogo=Array{Tuple,1}()

	for i in ocean.pccs
		txt=text(string(i.num),12,:black,:right)
		#txt=text(string(1),12,:black,:right)
		push!(pcc,(i.node.gps.lng,i.node.gps.lat,txt))
	end

	for i in ocean.owpps
		txt=text(string(i.num),12,:red,:right)
		push!(gen,(i.node.gps.lng,i.node.gps.lat,txt))
	end

	for i in ocean.bndryPnts
		push!(bnd,(i.gps.lng,i.gps.lat))
	end
	push!(bnd,(ocean.bndryPnts[1].gps.lng,ocean.bndryPnts[1].gps.lat))

	for h in ocean.nogos
		for i in h.nodes
			push!(ng,(i.gps.lng,i.gps.lat))
		end
		push!(ng,(h.nodes[1].gps.lng,h.nodes[1].gps.lat))
	end
	os=0

	Xpcc=[x[1] for x in pcc]
	Xgen=[x[1] for x in gen]
	Xbnd=[x[1] for x in bnd]
	Xnogo=[x[1] for x in ng]

	Ypcc=[x[2] for x in pcc]
	Ygen=[x[2] for x in gen]
	Ybnd=[x[2] for x in bnd]
	Ynogo=[x[2] for x in nogo]
	xlimax=trunc(Int,findmax(findmax([Xnogo,Xbnd,Xpcc,Xgen])[1])[1])+os
	ylimax=trunc(Int,findmax(findmax([Ynogo,Ybnd,Ypcc,Ygen])[1])[1])+os
	xlimin=trunc(Int,findmin(findmin([Xnogo,Xbnd,Xpcc,Xgen])[1])[1])-os
	ylimin=trunc(Int,findmin(findmin([Ynogo,Ybnd,Ypcc,Ygen])[1])[1])-os

	plotly()
	p=plot(Xpcc,Ypcc,annotations=pcc,color = :blue,seriestype=:scatter,label="PCC",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	plot!(p,Xgen,Ygen,annotations=gen,color = :red,seriestype=:scatter,label="OWPP")
	#plot!(p,Xbnd,Ybnd,color = :green,seriestype=:scatter,label="BND")
	plot!(p,Xbnd,Ybnd,color = :green,label="")
	#plot!(p,Xnogo,Ynogo,color = :black,seriestype=:scatter,label="NOGO")
	plot!(p,Xnogo,Ynogo,color = :black,label="")
	p
end
function ppf_saveSystem(ocn,nme)
	if (isdir("tempFiles/data/solutions/") == false)
		mkdir("tempFiles/data/solutions")
    else
    end

    #save cicuits
    save("tempFiles/data/solutions/"*nme*".jld2", "ocean",ocn)
end

function ppf_saveCircuit(oss_systems,nme)

    #check appropriate directories exist
	if (isdir("tempFiles/data/circuits/") == false)
		mkdir("tempFiles/data/circuits")
    else
    end

    #save cicuits
    save("tempFiles/data/circuits/circ_"*nme*".jld2", "circuits",oss_systems)
#	load("tempFiles/data/circuits/CRC_"*string(1)*string(2)*".jld2")["circuits"]
end
########################################### depricated #######################################


#=
function ppf_printOcnXYt(ocean,pth,goal)
	Owpp_peri=Array{nogo,1}()
	pcc=Array{Tuple,1}()
	gen=Array{Tuple,1}()
	bnd=Array{Tuple,1}()
	domOwpps_peri=Array{Tuple,1}()
	ng=Array{Tuple,1}()
	domNodes=Array{Tuple,1}()
	domOwpps=Array{Tuple,1}()
	domOwppsMV=Array{Tuple,1}()
	domOwpps_MVperi=Array{Tuple,1}()


	for i in ocean.pccs
		txt=text(string(i.num),12,:black,:right)
		#txt=text(string(1),12,:black,:right)
		push!(pcc,(i.node.xy.x,i.node.xy.y,txt))
	end

	for i in ocean.owpps
		dummy_peri=nogo()
		txt=text(string(i.num),12,:red,:right)
		push!(gen,(i.node.xy.x,i.node.xy.y,txt))
		dummy_node=node()
		dummy_node.xy.x=i.node.xy.x-i.zone.neg_width
		dummy_node.xy.y=i.node.xy.y-i.zone.neg_height
		push!(dummy_peri.nodes,deepcopy(dummy_node))


		dummy_node=node()
		dummy_node.xy.x=i.node.xy.x+i.zone.pos_width
		dummy_node.xy.y=i.node.xy.y-i.zone.neg_height
		push!(dummy_peri.nodes,deepcopy(dummy_node))


		dummy_node=node()
		dummy_node.xy.x=i.node.xy.x+i.zone.pos_width
		dummy_node.xy.y=i.node.xy.y+i.zone.pos_height
		push!(dummy_peri.nodes,deepcopy(dummy_node))

		dummy_node=node()
		dummy_node.xy.x=i.node.xy.x-i.zone.neg_width
		dummy_node.xy.y=i.node.xy.y+i.zone.pos_height
		push!(dummy_peri.nodes,deepcopy(dummy_node))

		push!(Owpp_peri,deepcopy(dummy_peri))
	end

	for i in ocean.bndryPnts
		push!(bnd,(i.xy.x,i.xy.y))
	end
	push!(bnd,(ocean.bndryPnts[1].xy.x,ocean.bndryPnts[1].xy.y))

	for i in ocean.discretedom.nodes
		push!(domNodes,(i.xy.x,i.xy.y))
	end

	for i in ocean.owpps
		#i=ocean.owpps[length(ocean.owpps)]
		#i=ocean.owpps[1]
		#for j in i.zone.pnts
		#	push!(domOwpps,(j.xy.x,j.xy.y))
		#end
		for j in i.zone.nodes
			push!(domOwpps_peri,(j.xy.x,j.xy.y))
		end
		for j in i.mv_zone.nodes
			push!(domOwppsMV,(j.xy.x,j.xy.y))
		end
		#for j in i.mv_zone.pnts
		#	push!(domOwpps_MVperi,(j.xy.x,j.xy.y))
		#end
	end

	for h in ocean.nogos
		for i in h.nodes
			push!(ng,(i.xy.x,i.xy.y))
		end
		push!(ng,(h.nodes[1].xy.x,h.nodes[1].xy.y))
	end
	os=0

	Xpcc=[x[1] for x in pcc]
	Xgen=[x[1] for x in gen]
	Xbnd=[x[1] for x in bnd]
	Xnogo=[x[1] for x in ng]
	Xdom=[x[1] for x in domNodes]
	Xowpp=[x[1] for x in domOwpps]
	Xowpp_peri=[x[1] for x in domOwpps_peri]
	XowppMV=[x[1] for x in domOwppsMV]
	Xowpp_MVperi=[x[1] for x in domOwpps_MVperi]

	Ypcc=[x[2] for x in pcc]
	Ygen=[x[2] for x in gen]
	Ybnd=[x[2] for x in bnd]
	Ynogo=[x[2] for x in ng]
	Ydom=[x[2] for x in domNodes]
	Yowpp=[x[2] for x in domOwpps]
	Yowpp_peri=[x[2] for x in domOwpps_peri]
	YowppMV=[x[2] for x in domOwppsMV]
	Yowpp_MVperi=[x[2] for x in domOwpps_MVperi]

	xlimax=trunc(Int,findmax(findmax([Xnogo,Xbnd,Xpcc,Xgen])[1])[1])+os
	ylimax=trunc(Int,findmax(findmax([Ynogo,Ybnd,Ypcc,Ygen])[1])[1])+os
	xlimin=trunc(Int,findmin(findmin([Xnogo,Xbnd,Xpcc,Xgen])[1])[1])-os
	ylimin=trunc(Int,findmin(findmin([Ynogo,Ybnd,Ypcc,Ygen])[1])[1])-os

	plotly()
	p=plot(Xpcc,Ypcc,annotations=pcc,color = :green,seriestype=:scatter,label="PCC",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	plot!(p,Xdom,Ydom,color = :black,seriestype=:scatter,markersize = 1,label="")
	plot!(p,XowppMV,YowppMV,color = :yellow,seriestype=:scatter,markersize = 2,label="")
	#plot!(p,Xowpp_MVperi,Yowpp_MVperi,color = :orange,seriestype=:scatter,markersize = 2,label="")
	plot!(p,Xowpp,Yowpp,color = :blue,seriestype=:scatter,markersize = 2,label="")
	plot!(p,Xowpp_peri,Yowpp_peri,color = :red,seriestype=:scatter,markersize = 2,label="")
	plot!(p,Xgen,Ygen,annotations=gen,color = :blue,seriestype=:scatter,label="OWPP")

	#plot!(p,Xbnd,Ybnd,color = :green,seriestype=:scatter,label="BND")
	plot!(p,Xbnd,Ybnd,color = :green,label="")
	#plot!(p,Xnogo,Ynogo,color = :black,seriestype=:scatter,label="NOGO")
	plot!(p,Xnogo,Ynogo,color = :black,label="")

	for nds in Owpp_peri
		Xperi=Array{Float32,1}()
		Yperi=Array{Float32,1}()
		for nd in nds.nodes
			push!(Xperi,deepcopy(nd.xy.x))
			push!(Yperi,deepcopy(nd.xy.y))
		end
		push!(Xperi,deepcopy(nds.nodes[1].xy.x))
		push!(Yperi,deepcopy(nds.nodes[1].xy.y))
		plot!(p,Xperi,Yperi,color = :black,label="")
	end

#=
		for edge in ocean.discretedom.edges
			Xedge=Array{Float32,1}()
			Yedge=Array{Float32,1}()
			push!(Xedge,ocean.discretedom.nodes[edge.tail].xy.x)
			push!(Yedge,ocean.discretedom.nodes[edge.tail].xy.y)
			push!(Xedge,ocean.discretedom.nodes[edge.head].xy.x)
			push!(Yedge,ocean.discretedom.nodes[edge.head].xy.y)
			#push!(nm,text(string(edge.tail),5,:black,:right))
			plot!(p,Xedge,Yedge,color = :orange,label="")
		end
		for owpp in ocean.owpps
			Xedge=Array{Float32,1}()
			Yedge=Array{Float32,1}()
			for edge in owpp.node.edges
				push!(Xedge,ocean.discretedom.nodes[edge.tail].xy.x)
				push!(Yedge,ocean.discretedom.nodes[edge.tail].xy.y)
				push!(Xedge,ocean.discretedom.nodes[edge.head].xy.x)
				push!(Yedge,ocean.discretedom.nodes[edge.head].xy.y)
			end
			plot!(p,Xedge,Yedge,color = :orange,label="")
		end=#
		nd=pth
		full=pth.G_cost
		strtXY=pth.xy
		while nd.num != goal.num
			rough=full-nd.G_cost
			if(nd.parent.parent.parent.num != goal.num && nd.parent.parent.num != goal.num && nd.parent.num != goal.num)
				crow=lof_pnt2pnt_dist(nd.xy,nd.parent.parent.parent.xy)
				road=nd.G_cost-nd.parent.parent.parent.G_cost
				if (((road-crow)/crow) >= 0.08)
					println("road: "*string(road))
					println("crow: "*string(crow))
					println(string(nd.num)*" dist: "*string((road-crow)/crow)*"%")
				end
			end
			Xedge=Array{Float32,1}()
			Yedge=Array{Float32,1}()
			push!(Xedge,nd.xy.x)
			push!(Yedge,nd.xy.y)
			push!(Xedge,nd.parent.xy.x)
			push!(Yedge,nd.parent.xy.y)
			#push!(nm,text(string(edge.tail),5,:black,:right))
			#println("node: "*string(nd.num)*" is roughly "*string(rough)*" km from goal")
			#println("node: "*string(nd.num)*" is "*string(nd.H_cost)*" km from goal")
			#println(string((rough-nd.H_cost)/nd.H_cost)*"%")
			plot!(p,Xedge,Yedge,color = :blue,label="")
			nd=nd.parent
		end
	p
end

function ppf_printOcnXY(ocn,pths)

	plotly()
	p=plot([ocn.discretedom.nodes[1].xy.x],[ocn.discretedom.nodes[1].xy.y],color = :red,markersize=2,seriestype=:scatter,label="",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
	for indx = 2:length(ocn.discretedom.nodes)
		plot!(p,[ocn.discretedom.nodes[indx].xy.x],[ocn.discretedom.nodes[indx].xy.y],color = :red,markersize=2,seriestype=:scatter,label="")
	end

	for indx = 1:length(ocn.pccs)
		plot!(p,[ocn.pccs[indx].node.xy.x],[ocn.pccs[indx].node.xy.y],color = :green,seriestype=:scatter,label="")
	end

	for owpp in ocn.owpps
		plot!(p,[owpp.node.xy.x],[owpp.node.xy.y],color = :blue,seriestype=:scatter,label="")
	end

	#=for constraint in ocean.constrain.ellipses
		exs=Array{Float32,1}()
		wys_p=Array{Float32,1}()
		wys_n=Array{Float32,1}()
		for x=(constraint.x0-constraint.rx):ocean.sys.prec:(constraint.x0+constraint.rx)
			y=(constraint.ry/constraint.rx)*sqrt(constraint.rx^2-(x-constraint.x0)^2)
			push!(exs,x)
			push!(wys_n,constraint.y0-y)
			push!(wys_p,constraint.y0+y)

		end
		x=(constraint.x0+constraint.rx)
		y=(constraint.ry/constraint.rx)*sqrt(constraint.rx^2-(x-constraint.x0)^2)
		push!(exs,x)
		push!(wys_n,constraint.y0-y)
		push!(wys_p,constraint.y0+y)
		plot!(p,exs,wys_n,color = :purple,label="")
		plot!(p,exs,wys_p,color = :purple,label="")
	end=#

	#=for edge in ocean.discretedom.edges
		plot!(p,[ocean.discretedom.nodes[edge.head].xy.x,ocean.discretedom.nodes[edge.tail].xy.x],[ocean.discretedom.nodes[edge.head].xy.y,ocean.discretedom.nodes[edge.tail].xy.y],color = :blue,label="")
	end

	for owpp in ocean.owpps
		for edge in owpp.node.edges
			plot!(p,[ocean.discretedom.nodes[edge.head].xy.x,ocean.discretedom.nodes[edge.tail].xy.x],[ocean.discretedom.nodes[edge.head].xy.y,ocean.discretedom.nodes[edge.tail].xy.y],color = :blue,label="")
		end
	end=#

	for pth in pths
		nd=pth
		full=pth.G_cost
		strtXY=pth.xy
		goal=deepcopy(nd.goal)
		while nd.num != goal
			rough=full-nd.G_cost
			if(nd.parent.parent.parent.num != goal && nd.parent.parent.num != goal && nd.parent.num != goal)
				crow=lof_pnt2pnt_dist(nd.xy,nd.parent.parent.parent.xy)
				road=nd.G_cost-nd.parent.parent.parent.G_cost
				if (((road-crow)/crow) >= 0.08)
					println("road: "*string(road))
					println("crow: "*string(crow))
					println(string(nd.num)*" dist: "*string((road-crow)/crow)*"%")
				end
			end
			Xedge=Array{Float32,1}()
			Yedge=Array{Float32,1}()
			push!(Xedge,nd.xy.x)
			push!(Yedge,nd.xy.y)
			push!(Xedge,nd.parent.xy.x)
			push!(Yedge,nd.parent.xy.y)
			plot!(p,Xedge,Yedge,color = :black,label="")
			nd=nd.parent
		end
	end

	p
end

function ppf_saveSystem(ocn,nme)
	if (isdir("tempFiles/data/solutions/") == false)
		mkdir("tempFiles/data/solutions")
    else
    end

    #save cicuits
    save("tempFiles/data/solutions/"*nme*".jld2", "ocean",ocn)
end

function ppf_saveCircuit(oss_systems,first,last)
    #check appropriate directories exist
	if (isdir("tempFiles/data/circuits/") == false)
		mkdir("tempFiles/data/circuits")
    else
    end

    #save cicuits
    save("tempFiles/data/circuits/circ_"*string(first)*string(last)*".jld2", "circuits",oss_systems)
#	load("tempFiles/data/circuits/CRC_"*string(1)*string(2)*".jld2")["circuits"]
end


function ppf_loadFile()
    ab=Array{eez,1}()
    dm=Array{eez,1}()
    Ids=Array{Array{Int32,1},1}()
    obj=Array{Float32,1}()

    for i=1:18
		if isfile("../../results/setUp220kv_66kv/22066_"*string(i)*".jld") == true
	        solution=load("../../results/setUp220kv_66kv/22066_"*string(i)*".jld")
	        push!(Ids,solution["optIds"])
	        push!(obj,solution["objective"])
	        push!(dm,solution["domain"])
	        push!(ab,solution["asBuilt"])
		end
    end
    solutions=(ab,obj,Ids,dm)
    return solutions
end
=#
