#Imports the Corwind data and calculates load loss/ constrained energy
function wndF_wndPrf(profiles)
    wp=Array{Float64,1}()
    dummy=Array{Float64,1}()
    dummy=profiles[1]
    for i=2:length(profiles)
        dummy=dummy+profiles[i]
    end
    wp=dummy./length(profiles)

    wnd=wind()
    mx=findmax(wp)[1]
    ord=reverse(sort(wp))./mx
    is=Array{Int,1}()
    for i=1:length(ord)
            push!(is,i)
    end
    wndF_conEng([is ord],wnd,mx)#constraint energy calc
    wndF_ldLss([is ord], wnd)#calculates loss factor
    wnd
    return wnd
end

#Imports the Corwind data and calculates load loss/ constrained energy
function wndF_wndPrf_pu(profiles)
    wnd=wind()
    wnd.ce=profiles[1].ce
    wnd.pu=profiles[1].pu
    wnd.delta=profiles[1].delta
    wnd.lf=profiles[1].lf
    for i=2:length(profiles)
        wnd.ce=wnd.ce+profiles[i].ce
        wnd.pu=wnd.pu+profiles[i].pu
        wnd.delta=wnd.delta+profiles[i].delta
        wnd.lf=wnd.lf+profiles[i].lf
    end
    wnd.ce=(wnd.ce)./length(profiles)
    wnd.pu=(wnd.pu)./length(profiles)
    wnd.delta=(wnd.delta)./length(profiles)
    wnd.lf=(wnd.lf)./length(profiles)
    wnd
    return wnd
end


#=function wndF_wndPrf(nmes)
    prof = DataFrame(XLSX.readtable("layout//data.xlsx", "wind_data")...)

    wp=Array{Float64,1}()
    dummy=Array{Float64,1}()
    dummy=getproperty(prof,Symbol(nmes[1]))
    for i=2:length(nmes)
        i=2
        dummy=dummy+getproperty(prof,Symbol(nmes[i]))
    end
    wp=dummy./length(nmes)

    wnd=wind()
    mx=findmax(wp)[1]
    ord=reverse(sort(wp))./mx
    is=Array{Int,1}()
    for i=1:length(ord)
            push!(is,i)
    end
    wndF_conEng([is ord],wnd,mx)#constraint energy calc
    wndF_ldLss([is ord], wnd)#calculates loss factor
    wnd
    return wnd
end=#

#Calculates the loss factor associated with the wind profile
function wndF_ldLss(div, wind)
  wind.lf=(sum(div[:,2]))*0.85/length(div[:,2])#saves loss factor, 0.85 is wake penalization
  #loss factor/llf formula ref: Guidelines on the calculation and use of loss factors Te Mana Hiko Electricity Authority
  #0.85 for wake effect from Evaluation of the wind direction uncertainty and its impact on wake modeling at the Horns Rev offshore wind farm
#M. Gaumond  P.‐E. Réthoré  S. Ott  A. Peña  A. Bechmann  K. S. Hansen
  llf=0.0
  for pu in div[:,2]
    llf=llf+(pu*0.85)^2
    #llf=llf+(pu)^2
  end
  wind.delta=llf/length(div[:,2])#saves load loss factor
end

#calc for constrained energy
function wndF_conEng(graph,wind,max)
#create sized arrays
    B=zeros(length(graph[:,2]),2)
    conENG=zeros(length(graph[:,2]),2)
    p_div=polyfit(graph[:,1],graph[:,2],3)#make a polynomial approximation
    area=zeros(length(graph[:,2]),1)
    integral=polyint(p_div)#set up the integral
    area=(polyval(integral,graph[:,1])-(graph[:,1].*graph[:,2]))#take integral to find area under curve
    x_axis=reverse(graph[:,2],2)
    y_axis=reverse(area[:,1],2)#reverse x and y axis
    B=[x_axis y_axis]
    conENG=sortslices(B,dims=1)#sort by x axis
    wind.ce=conENG[:,2]
    wind.pu=conENG[:,1]#store pu and constraind energy in wind object
    return nothing
end
