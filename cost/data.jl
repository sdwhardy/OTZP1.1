#=
In this file all input cost data is set
=#

#Sets the exchange rate
#original input data is in 2008 British Pounds
function cstD_xchg()
    return 1.16
end

#Sets values of all cost factors discribed below
#the structure of object ks is described in file cst_structure.jl
function cstD_cfs()
    sht = XLSX.readxlsx("layout//data.xlsx")["cost_data"]
    ks=cstS_ks()#create an instance of object ks
    ks.FC_ac=sht["B1"]#fixed AC cost
    ks.FC_dc=sht["B2"]#fixed DC cost
    ks.dc=sht["B3"]#penalization factor for different than 2 xfrms
    ks.f_ct=sht["B4"]#generating plant variable cost
    ks.p_ct=sht["B5"]#substructure variable cost
    ks.c_ct=sht["B6"]#hvdc converter variable cost
    ks.Qc_oss=sht["B7"]#M£/MVAr
    ks.Qc_pcc=sht["B8"]#M£/MVAr
    ks.life=sht["B9"]#lifetime of wind farm
    ks.T_op=sht["B10"]#Operational lifetime in hours
    ks.E_op=sht["B11"]#Energy price £/Wh
    ks.cf=sht["B12"]#Capitalization factor
    return ks
end



#################################################### Depricated functions below ####################################

#Sets Cost factor for cable compensation
#bin is a binary variable that which is true if a cable connects 2 OSS and false if from OSS to PCC
#ks is a set of cost factors set within this file
#=function cstD_QC(bin,ks)
    p2e=cstD_xchg()
    #offshore(OSS) to onshore(PCC) connection
    ks.Qc_oss=0.025*p2e#M£/MVAr
    ks.Qc_pcc=0.015*p2e#M£/MVAr
    #offshore(OSS) to offshore(OSS) connection
    if bin == true
        ks.Qc_pcc=ks.Qc_oss
    end
    return nothing
end=#


##############################
########## PCCs ##############
#=function lod_Belgium()
    return true
end

function lod_pccGps()
    pcc=Array{Tuple,1}()
    if lod_Belgium()==true
    ##### Belgium #####
        push!(pcc,(2.939692,51.239737))
        push!(pcc,(3.183611,51.32694))
    else
        ################################## India ###################################
        #Dui
        push!(pcc,(71.0,20.71))
        #Barbarkot (Ultratech)
        push!(pcc,(71.399491,20.866641))
    end
    return pcc
end
##############################
########## Gens ##############

function lod_gensGps()
    c=Array{Tuple,1}()
    wnd=Array{String,1}()
    p=Array{Float32,1}()
##### Belgium #####
    if lod_Belgium()==true
        #Norther
        push!(c,(3.015833,51.52806))
        push!(p,250.0)
        push!(wnd,"Norther")
        #Thornton
        push!(c,((2.97+2.919972)/2,(51.56+51.53997)/2))
        push!(p,250.0)
        push!(wnd,"Thornton")
        #Rentel
        push!(c,(2.939972,51.59))
        push!(p,250.0)
        push!(wnd,"Rentel")
        #Northwind
        push!(c,(2.900972,51.61897))
        push!(p,250.0)
        push!(wnd,"Northwind")
        #Seastar
        push!(c,(2.859972,51.63))
        push!(p,250.0)
        push!(wnd,"Seastar")
        #Nobelwind/Belwind
        push!(c,((2.819972+2.799972)/2,(51.664+51.67)/2))
        push!(p,250.0)
        push!(wnd,"Nobelwind")
        #Northwester
        push!(c,(2.757,51.68597))
        push!(p,250.0)
        push!(wnd,"Northwester")
        #Mermaid
        push!(c,(2.74,51.71997))
        push!(p,250.0)
        push!(wnd,"Mermaid")
    else
        #A5
        push!(c,(71.375,20.53))
        push!(p,500.0)
        push!(wnd,"A5")
        #A1
        push!(c,(71.4,20.58))
        push!(p,500.0)
        push!(wnd,"A1")
        #A2
        push!(c,(71.53,20.55))
        push!(p,500.0)
        push!(wnd,"A2")
        #A3
        push!(c,(71.5,20.64))
        push!(p,500.0)
        push!(wnd,"A3")
        #A4
        push!(c,(71.435,20.53))
        push!(p,500.0)
        push!(wnd,"A4")
        #A6
        push!(c,(71.525,20.45))
        push!(p,500.0)
        push!(wnd,"A6")
        #A7
        push!(c,(71.44,20.475))
        push!(p,500.0)
        push!(wnd,"A7")
        #A8
        push!(c,(71.36,20.475))
        push!(p,500.0)
        push!(wnd,"A8")
        #A9
        push!(c,(71.27,20.475))
        push!(p,500.0)
        push!(wnd,"A9")
        #A10
        push!(c,(71.15,20.47))
        push!(p,500.0)
        push!(wnd,"A10")
        #A11
        push!(c,(71.11,20.39))
        push!(p,500.0)
        push!(wnd,"A11")
        #A12
        push!(c,(71.28,20.415))
        push!(p,500.0)
        push!(wnd,"A12")
        #A13
        push!(c,(71.25,20.365))
        push!(p,500.0)
        push!(wnd,"A13")
        #A14
        push!(c,(71.425,20.415))
        push!(p,500.0)
        push!(wnd,"A14")
        #A15
        push!(c,(71.37,20.365))
        push!(p,500.0)
        push!(wnd,"A15")
        #A16
        push!(c,(71.5,20.365))
        push!(p,500.0)
        push!(wnd,"A16")
        #A17
        push!(c,(71.5,20.32))
        push!(p,500.0)
        push!(wnd,"A17")
        #A18
        push!(c,(71.35,20.32))
        push!(p,500.0)
        push!(wnd,"A18")
        #A19
        push!(c,(71.45,20.27))
        push!(p,500.0)
        push!(wnd,"A19")
    end
    return c,p,wnd
end

#Grid PU mva
function lod_cnceMva()
    if lod_Belgium()==true
        mva=250.0
    else
        mva=500.0
    end
    return mva
end

#set onshore transmission voltage
function lod_pccKv()
    return 220.0
end
=#
