
# Reading data in from Excel file
function get_PccData()
    df = DataFrame(XLSX.readtable("v2.0/layout/input_data/data.xlsx", "pcc_data")...)
    pccs=Array{bus,1}()
    for index=1:length(df[:, 1])
        dummy_bus=bus()
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.kV=df.kv[index]
        push!(pccs,deepcopy(dummy_bus))
    end
    return pccs
end
#**
#**#reads OWPP data in excel file
function get_OwppData()
    df = DataFrame(XLSX.readtable("v2.0/layout/input_data/data.xlsx", "owpp_data")...)
    owpps=Array{bus,1}()
    for index=1:length(df[:, 1])
        dummy_bus=bus()
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.name=df.name[index]
        dummy_bus.mva=df.mva[index]
        push!(owpps,deepcopy(dummy_bus))
    end
    get_WindData(owpps)
    return owpps
end
#**Reads Wind data from excel file
function get_WindData(owpps)
    df = DataFrame(XLSX.readtable("v2.0/layout/input_data/data.xlsx", "wind_data")...)
    for owpp in owpps
        wnd=wndF_wndPrf([getproperty(df,Symbol(owpp.name))])
        owpp.wnd=deepcopy(wnd)
    end
end
#**
#**reads system data from Excel file
function get_SystemData()
    sht = XLSX.readxlsx("v2.0/layout/input_data/data.xlsx")["sys_data"]
    sys=system()
    sys.nogoNum=sht["B1"]
    sys.prec=sht["B2"]
    sys.mwPerKm=sht["B3"]
    sys.mvCl=sht["B4"]
    return sys
end
