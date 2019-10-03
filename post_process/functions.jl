function ppf_printOcnGPS(ocean)
	pcc=Array{Tuple,1}()
	gen=Array{Tuple,1}()
	bnd=Array{Tuple,1}()
	nogo=Array{Tuple,1}()

	for i in ocean.pccs
		txt=text(string(i.num),12,:black,:right)
		#txt=text(string(1),12,:black,:right)
		push!(pcc,(i.cntr.gps.lng,i.cntr.gps.lat,txt))
	end

	for i in ocean.owpps
		txt=text(string(i.num),12,:red,:right)
		push!(gen,(i.cntr.gps.lng,i.cntr.gps.lat,txt))
	end

	for i in ocean.bndryPnts
		push!(bnd,(i.cntr.gps.lng,i.cntr.gps.lat))
	end

	for h in ocean.nogos
		for i in h.nodes
			push!(nogo,(i.cntr.gps.lng,i.cntr.gps.lat))
		end
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
	plot!(p,Xbnd,Ybnd,color = :green,seriestype=:scatter,label="BND")
	plot!(p,Xnogo,Ynogo,color = :black,seriestype=:scatter,label="BND")
	p
end

function ppf_printOcnXY(ocean)
	pcc=Array{Tuple,1}()
	gen=Array{Tuple,1}()
	bnd=Array{Tuple,1}()
	nogo=Array{Tuple,1}()

	for i in ocean.pccs
		txt=text(string(i.num),12,:black,:right)
		#txt=text(string(1),12,:black,:right)
		push!(pcc,(i.cntr.xy.x,i.cntr.xy.y,txt))
	end

	for i in ocean.owpps
		txt=text(string(i.num),12,:red,:right)
		push!(gen,(i.cntr.xy.x,i.cntr.xy.y,txt))
	end

	for i in ocean.bndryPnts
		push!(bnd,(i.cntr.xy.x,i.cntr.xy.y))
	end
	push!(bnd,(ocean.bndryPnts[1].cntr.xy.x,ocean.bndryPnts[1].cntr.xy.y))
	os=0

	for h in ocean.nogos
		for i in h.nodes
			push!(nogo,(i.cntr.xy.x,i.cntr.xy.y))
		end
		push!(nogo,(h.nodes[1].cntr.xy.x,h.nodes[1].cntr.xy.y))
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
