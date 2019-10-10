
function ppf_printOcnXY(ocean)
	pcc=Array{Tuple,1}()
	gen=Array{Tuple,1}()
	bnd=Array{Tuple,1}()
	nogo=Array{Tuple,1}()
	domNodes=Array{Tuple,1}()
	domOwpps=Array{Tuple,1}()
	domOwpps_peri=Array{Tuple,1}()
	domOwppsMV=Array{Tuple,1}()
	domOwpps_MVperi=Array{Tuple,1}()

	for i in ocean.pccs
		txt=text(string(i.num),12,:black,:right)
		#txt=text(string(1),12,:black,:right)
		push!(pcc,(i.node.xy.x,i.node.xy.y,txt))
	end

	for i in ocean.owpps
		txt=text(string(i.num),12,:red,:right)
		push!(gen,(i.node.xy.x,i.node.xy.y,txt))
	end

	for i in ocean.bndryPnts
		push!(bnd,(i.xy.x,i.xy.y))
	end
	push!(bnd,(ocean.bndryPnts[1].xy.x,ocean.bndryPnts[1].xy.y))

	for i in ocean.discretedom.nodes
		push!(domNodes,(i.xy.x,i.xy.y))
	end

	#for i in ocean.owpps
		i=ocean.owpps[length(ocean.owpps)-4]
		for j in i.zone.pnts
			push!(domOwpps,(j.xy.x,j.xy.y))
		end
		for j in i.zone.periPnts
			push!(domOwpps_peri,(j.xy.x,j.xy.y))
		end
		for j in i.mv_zone.pnts
			push!(domOwppsMV,(j.xy.x,j.xy.y))
		end
		for j in i.mv_zone.periPnts
			push!(domOwpps_MVperi,(j.xy.x,j.xy.y))
		end
	#end

	for h in ocean.nogos
		for i in h.nodes
			push!(nogo,(i.xy.x,i.xy.y))
		end
		push!(nogo,(h.nodes[1].xy.x,h.nodes[1].xy.y))
	end
	os=0

	Xpcc=[x[1] for x in pcc]
	Xgen=[x[1] for x in gen]
	Xbnd=[x[1] for x in bnd]
	Xnogo=[x[1] for x in nogo]
	Xdom=[x[1] for x in domNodes]
	Xowpp=[x[1] for x in domOwpps]
	Xowpp_peri=[x[1] for x in domOwpps_peri]
	XowppMV=[x[1] for x in domOwppsMV]
	Xowpp_MVperi=[x[1] for x in domOwpps_MVperi]

	Ypcc=[x[2] for x in pcc]
	Ygen=[x[2] for x in gen]
	Ybnd=[x[2] for x in bnd]
	Ynogo=[x[2] for x in nogo]
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
	plot!(p,Xowpp_MVperi,Yowpp_MVperi,color = :orange,seriestype=:scatter,markersize = 2,label="")
	plot!(p,Xowpp,Yowpp,color = :blue,seriestype=:scatter,markersize = 2,label="")
	plot!(p,Xowpp_peri,Yowpp_peri,color = :red,seriestype=:scatter,markersize = 2,label="")
	plot!(p,Xgen,Ygen,annotations=gen,color = :blue,seriestype=:scatter,label="OWPP")

	#plot!(p,Xbnd,Ybnd,color = :green,seriestype=:scatter,label="BND")
	plot!(p,Xbnd,Ybnd,color = :green,label="")
	#plot!(p,Xnogo,Ynogo,color = :black,seriestype=:scatter,label="NOGO")
	plot!(p,Xnogo,Ynogo,color = :black,label="")
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
			push!(nogo,(i.gps.lng,i.gps.lat))
		end
		push!(nogo,(h.nodes[1].gps.lng,h.nodes[1].gps.lat))
	end
	os=0

	Xpcc=[x[1] for x in pcc]
	Xgen=[x[1] for x in gen]
	Xbnd=[x[1] for x in bnd]
	Xnogo=[x[1] for x in nogo]

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
