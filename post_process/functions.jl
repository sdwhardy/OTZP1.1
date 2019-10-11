
function ppf_printOcnXY(ocean)
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
		for j in i.zone.pnts
			push!(domOwpps_peri,(j.xy.x,j.xy.y))
		end
		for j in i.mv_zone.pnts
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
		Xperi=Array{Float64,1}()
		Yperi=Array{Float64,1}()
		for nd in nds.nodes
			push!(Xperi,deepcopy(nd.xy.x))
			push!(Yperi,deepcopy(nd.xy.y))
		end
		push!(Xperi,deepcopy(nds.nodes[1].xy.x))
		push!(Yperi,deepcopy(nds.nodes[1].xy.y))
		plot!(p,Xperi,Yperi,color = :black,label="")
	end
	p
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
