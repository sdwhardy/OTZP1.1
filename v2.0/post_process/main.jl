using Plots
function plot_ocean(ocn)
    plotly()
    p=plot()
    p=plot_PCCnOWPP(p,ocn)
    p
    gui()
end


function plot_circuit(circ)
    plotly()
    p=plot()
    p=plot_PCCnOWPP(p,circ)
    p=plot_OSSnLINES(p,circ)
    p
    gui()
	print_circuit_Details(circ)
end


function plot_OSSnLINES(p,circ)#eex or bus as arg
    for cbl in circ.MVcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.path
				push!(xs,pth.xy.x)
				push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :blue,linewidth = 2,label="")
	end
	for cbl in circ.HVcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.path
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,linewidth = 2,label="")
	end

	for cbl in circ.O2Ocbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.path
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,linewidth = 2,label="")
	end
	for cbl in circ.PCCcbls
		xs=Float32[]
		ys=Float32[]
		for pth in cbl.path
			push!(xs,pth.xy.x)
			push!(ys,pth.xy.y)
		end
		plot!(p,xs,ys,color = :black,linewidth = 2,label="")
	end
	xs=Float32[]
	ys=Float32[]
	tx=Tuple[]

	for (i,oss) in enumerate(circ.oss)
		txt=text(string(""),15,:red,:right)
		push!(tx,(oss.node.xy.x,oss.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :black,seriestype=:scatter,label="")
	xs=Float32[]
	ys=Float32[]
	tx=Tuple[]
	for (i,oss) in enumerate(circ.mog)
		txt=text(string(""),15,:black,:right)
		push!(tx,(oss.node.xy.x,oss.node.xy.y,txt))
	end
	xs=[x[1] for x in tx]
	ys=[x[2] for x in tx]
	plot!(p,xs,ys,annotation=tx,color = :black,seriestype=:scatter,label="")
    return p
end



function plot_PCCnOWPP(p,ocn)#eex or bus as arg

    plot!(p,[ocn.pcc.node.xy.x],[ocn.pcc.node.xy.y],color = :green,markersize=2,seriestype=:scatter,label="",xaxis = ("km", font(15, "Courier")),yaxis = ("km", font(15, "Courier")))
    plot!(p,[ocn.pcc.node.xy.x],[ocn.pcc.node.xy.y],color = :green,seriestype=:scatter,label="")
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
    return p
end



function print_hvac_canditate_mog2pcc_result(ac_sys)
    #platform
    println("Offshore platform:")
    println("CAPEX: "*string(ac_sys["plat_aft_ac"].costs.cpx))
    println("OPEX: "*string(ac_sys["plat_aft_ac"].costs.cm))
    println("TTl: "*string(ac_sys["plat_aft_ac"].costs.ttl))
    println()
    #offshore transformer
    println("Offshore Xfo: ")
    println("CAPEX: "*string(ac_sys["xfo_aft"].costs.cpx_p+ac_sys["xfo_aft"].costs.cpx_i))
    println("OPEX: "*string(ac_sys["xfo_aft"].costs.cm+ac_sys["xfo_aft"].costs.eens+ac_sys["xfo_aft"].costs.tlc))
    println("TTL: "*string(ac_sys["xfo_aft"].costs.ttl))
    println()
    #cable
    println("Cable:")
    if (ac_sys["cable"].mpc_ac==true)
        println("Mid point compensated cable: All but total is for 1/2 the cable required.")
    end
    println("CAPEX: "*string(ac_sys["cable"].costs.cpx_p+ac_sys["cable"].costs.cpx_i))
    println("S.G.: "*string(ac_sys["cable"].costs.sg))
    println("Compensation: "*string(ac_sys["cable"].costs.qc))
    println("OPEX: "*string(ac_sys["cable"].costs.rlc+ac_sys["cable"].costs.cm+ac_sys["cable"].costs.eens))
    println("TTL: "*string(ac_sys["cable"].costs.grand_ttl))
    println()
    #onshore transformer
    println("Onshore Xfo:")
    println("CAPEX: "*string(ac_sys["xfo_for"].costs.cpx_p+ac_sys["xfo_for"].costs.cpx_i))
    println("OPEX: "*string(ac_sys["xfo_for"].costs.cm+ac_sys["xfo_for"].costs.eens+ac_sys["xfo_for"].costs.tlc))
    println("TTL: "*string(ac_sys["xfo_for"].costs.ttl))
    println()
    #totals
    println("System:")
    println("CAPEX: "*string(ac_sys["cable"].costs.cpx_p+ac_sys["cable"].costs.cpx_i+ac_sys["cable"].costs.sg
    +ac_sys["cable"].costs.qc+ac_sys["xfo_aft"].costs.cpx_p+ac_sys["xfo_aft"].costs.cpx_i+ac_sys["xfo_for"].costs.cpx_p
    +ac_sys["xfo_for"].costs.cpx_i+ac_sys["plat_aft_ac"].costs.cpx))
    println("OPEX: "*string(ac_sys["cable"].costs.rlc+ac_sys["cable"].costs.cm+ac_sys["cable"].costs.eens
    +ac_sys["xfo_for"].costs.cm+ac_sys["xfo_for"].costs.tlc+ac_sys["xfo_for"].costs.eens
    +ac_sys["xfo_aft"].costs.cm+ac_sys["xfo_aft"].costs.tlc+ac_sys["xfo_aft"].costs.eens
    +ac_sys["plat_aft_ac"].costs.cm))
    println("TTL: "*string(ac_sys["cable"].costs.grand_ttl+ac_sys["xfo_aft"].costs.ttl+ac_sys["xfo_for"].costs.ttl
    +ac_sys["plat_aft_ac"].costs.ttl))
end



function print_hvdc_canditate_mog2pcc_result(dc_sys)
    #platform
    println("Offshore platform:")
    println("CAPEX: "*string(dc_sys["plat_aft_ac"].costs.cpx+dc_sys["plat_aft_dc"].costs.cpx))
    println("OPEX: "*string(dc_sys["plat_aft_ac"].costs.cm+dc_sys["plat_aft_dc"].costs.cm))
    println("TTl: "*string(dc_sys["plat_aft_ac"].costs.ttl+dc_sys["plat_aft_dc"].costs.ttl))
    println()
    #offshore transformer
    println("Offshore Xfo: ")
    println("CAPEX: "*string(dc_sys["xfo_aft"].costs.cpx_p+dc_sys["xfo_aft"].costs.cpx_i))
    println("OPEX: "*string(dc_sys["xfo_aft"].costs.cm+dc_sys["xfo_aft"].costs.eens+dc_sys["xfo_aft"].costs.tlc))
    println("TTL: "*string(dc_sys["xfo_aft"].costs.ttl))
    println()
    #offshore Converter
    println("Offshore Converter:")
    println("CAPEX: "*string(dc_sys["conv_aft"].costs.cpx))
    println("OPEX: "*string(dc_sys["conv_aft"].costs.tlc+dc_sys["conv_aft"].costs.eens+dc_sys["conv_aft"].costs.cm))
    println("TTL: "*string(dc_sys["conv_aft"].costs.ttl))
    println()
    #cable
    println("Cable:")
    println("CAPEX: "*string(dc_sys["cable"].costs.cpx_p+dc_sys["cable"].costs.cpx_i))
    println("OPEX: "*string(dc_sys["cable"].costs.rlc+dc_sys["cable"].costs.cm+dc_sys["cable"].costs.eens))
    println("TTL: "*string(dc_sys["cable"].costs.ttl))
    println()
    #onshore transformer
    println("Onshore Xfo:")
    println("CAPEX: "*string(dc_sys["xfo_for"].costs.cpx_p+dc_sys["xfo_for"].costs.cpx_i))
    println("OPEX: "*string(dc_sys["xfo_for"].costs.cm+dc_sys["xfo_for"].costs.eens+dc_sys["xfo_for"].costs.tlc))
    println("TTL: "*string(dc_sys["xfo_for"].costs.ttl))
    println()
    #onshore Converter
    println("Onshore Converter:")
    println("CAPEX: "*string(dc_sys["conv_for"].costs.cpx))
    println("OPEX: "*string(dc_sys["conv_for"].costs.tlc+dc_sys["conv_for"].costs.eens+dc_sys["conv_for"].costs.cm))
    println("TTL: "*string(dc_sys["conv_for"].costs.ttl))
    println()
    #totals
    println("System:")
    println("CAPEX: "*string(dc_sys["cable"].costs.cpx_p+dc_sys["cable"].costs.cpx_i+dc_sys["conv_aft"].costs.cpx
    +dc_sys["conv_for"].costs.cpx+dc_sys["xfo_aft"].costs.cpx_p+dc_sys["xfo_aft"].costs.cpx_i+dc_sys["xfo_for"].costs.cpx_p
    +dc_sys["xfo_for"].costs.cpx_i+dc_sys["plat_aft_ac"].costs.cpx+dc_sys["plat_aft_dc"].costs.cpx))
    println("OPEX: "*string(dc_sys["cable"].costs.rlc+dc_sys["cable"].costs.cm+dc_sys["cable"].costs.eens
    +dc_sys["conv_aft"].costs.tlc+dc_sys["conv_aft"].costs.eens+dc_sys["conv_aft"].costs.cm
    +dc_sys["conv_for"].costs.tlc+dc_sys["conv_for"].costs.eens+dc_sys["conv_for"].costs.cm
    +dc_sys["xfo_for"].costs.cm+dc_sys["xfo_for"].costs.tlc+dc_sys["xfo_for"].costs.eens
    +dc_sys["xfo_aft"].costs.cm+dc_sys["xfo_aft"].costs.tlc+dc_sys["xfo_aft"].costs.eens
    +dc_sys["plat_aft_ac"].costs.cm+dc_sys["plat_aft_dc"].costs.cm))
    println("TTL: "*string(dc_sys["cable"].costs.ttl+dc_sys["conv_aft"].costs.ttl+
    dc_sys["conv_for"].costs.ttl+dc_sys["xfo_aft"].costs.ttl+dc_sys["xfo_for"].costs.ttl
    +dc_sys["plat_aft_ac"].costs.ttl+dc_sys["plat_aft_dc"].costs.ttl))
end



function print_circuit_Details(circ)
	println()
	print("MVcbls: ")
	for mv_cbl in circ.MVcbls
		print(mv_cbl.length)
		print("km - ")
		print(mv_cbl.elec.volt)
		print("kV - ")
		print(mv_cbl.num)
		print(" - ")
		print(mv_cbl.size)
		print("mm - ")
		print(mv_cbl.elec.mva)
		print("mva |")
        print(mv_cbl.costs.grand_ttl)
		print("|ME")
		println()
    end
	println()
	print("HVcbls: ")
    for hv_cbl in circ.HVcbls
		print(hv_cbl.length)
		print("km - ")
		print(hv_cbl.elec.volt)
		print("kV - ")
		print(hv_cbl.num)
		print(" - ")
		print(hv_cbl.size)
		print("mm - ")
		print(hv_cbl.elec.mva)
		print("mva |")
		print(hv_cbl.costs.grand_ttl)
		print("|ME")
		println()
		if (hv_cbl.mpc_ac==true)
			print("Mid Point Compensation: ")
			print(hv_cbl.plat.acdc)
			print(" - ")
			print(hv_cbl.plat.mva)
			print("mva |")
            print(hv_cbl.plat.costs.ttl)
			print("|ME ")
			println()
		end
    end
	println()
	print("O2Ocbls: ")
    for hv_cbl in circ.O2Ocbls
		print(hv_cbl.length)
		print("km - ")
		print(hv_cbl.elec.volt)
		print("kV - ")
		print(hv_cbl.num)
		print(" - ")
		print(hv_cbl.size)
		print("mm - ")
		print(hv_cbl.elec.mva)
		print("mva |")
		print(hv_cbl.costs.grand_ttl)
		print("|ME ")
		println()
		if (hv_cbl.mpc_ac==true)
			print("Mid Point Compensation: ")
			print(hv_cbl.plat.acdc)
			print(" - ")
			print(hv_cbl.plat.mva)
			print("mva |")
            print(hv_cbl.plat.costs.ttl)
			print("|ME ")
			println()
		end
    end
	println()
	print("PCCcbls: ")
    for hv_cbl in circ.PCCcbls
		print(hv_cbl.length)
		print("km - ")
		print(hv_cbl.elec.volt)
		print("kV - ")
		print(hv_cbl.num)
		print(" - ")
		print(hv_cbl.size)
		print("mm - ")
		print(hv_cbl.elec.mva)
		print("mva |")
		print(hv_cbl.costs.grand_ttl)
		print("|ME ")
		println()
		if (hv_cbl.mpc_ac==true)
			print("Mid Point Compensation: ")
			print(hv_cbl.plat.acdc)
			print(" - ")
			print(hv_cbl.plat.mva)
			print("mva |")
            print(hv_cbl.plat.costs.ttl)
			print("|ME ")
			println()
		end
    end

    #total transformer and converter costs
	println()
	print("MOGs: ")
    for ss in circ.mog
		println()
		print("transformers: ")
        for xfo in ss.xfmrs
			print(xfo.num)
			print(" - ")
			print(xfo.elec.mva)
			print("mva |")
            print(xfo.costs.ttl)
			print("|ME ")
			println()
        end
		println()
		print("Converters: ")
        for cv in ss.conv
			print(cv.mva)
			print("mva |")
            print(cv.costs.ttl)
			print("|ME ")
			println()
        end
        #platforms per ss
		println()
		print("platforms: ")
        for plt in ss.plat
			print(plt.acdc)
			print(" - ")
			print(plt.mva)
			print("mva |")
            print(plt.costs.ttl)
			print("|ME ")
			println()
        end
    end
    #totals transformers and converters and platform costs
	println()
	print("OSSs: ")
    for ss in circ.oss
        #transformers per ss
		println()
		print("transformers: ")
        for xfo in ss.xfmrs
			print(xfo.num)
			print(" - ")
			print(xfo.elec.mva)
			print("mva |")
            print(xfo.costs.ttl)
			print("|ME ")
			println()
        end
		println()
        #converters per ss
		print("converters: ")
        for cv in ss.conv
			print(cv.mva)
			print("mva |")
            print(cv.costs.ttl)
			print("|ME ")
			println()
        end
		println()
        #platforms per ss
		print("platforms: ")
        for plt in ss.plat
			print(plt.acdc)
			print(" - ")
			print(plt.mva)
			print("mva |")
            print(plt.costs.ttl)
			print("|ME ")
			println()
        end
    end
    #transformers per pcc
	println()
	println("PCC: ")
	println("transformers: ")
    for xfo in circ.pcc.xfmrs
		print(xfo.num)
		print(" - ")
		print(xfo.elec.mva)
		print("mva |")
		print(xfo.costs.ttl)
		print("|ME ")
		println()
    end
	println()
    #converters per pcc
	print("Converters: ")
    for cv in circ.pcc.conv
		print(cv.mva)
		print("mva |")
		print(cv.costs.ttl)
		print("|ME ")
		println()
    end
	println()
	println()
	println("Total MVA: "*string(circ.mva))
	println("Total ME: "*string(circ.cost))
end
