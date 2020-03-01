function wndF_wndPrf(profiles)
    wp=Array{Float32,1}()
    dummy=Array{Float32,1}()
    dummy=profiles[1]
    for i=2:length(profiles)
        dummy=dummy+profiles[i]
    end
    wp=dummy./length(profiles)

    wnd=wind()
    mx=findmax(wp)[1]
    ord=reverse(sort(wp))./mx
    wndF_conEng(ord,wnd)#constraint energy calc
    wndF_ldLss(ord, wnd)#calculates loss factor
    wnd
    return wnd
end
#**loss factor
function wndF_ldLss(div, wind)
  wind.lf=(sum(div))*0.85/length(div)#saves loss factor, 0.85 is wake penalization
  #loss factor/llf formula ref: Guidelines on the calculation and use of loss factors Te Mana Hiko Electricity Authority
  #0.85 for wake effect from Evaluation of the wind direction uncertainty and its impact on wake modeling at the Horns Rev offshore wind farm
  #M. Gaumond  P.‐E. Réthoré  S. Ott  A. Peña  A. Bechmann  K. S. Hansen
  llf=0.0
  for pu in div
    llf=llf+(pu*0.85)^2
  end
  wind.delta=llf/length(div)#saves load loss factor
end

#**Finds Constraind energy
function wndF_conEng(graph,wnd)
    #create sized arrays
    ce=Float32[]
    for hr=1:length(graph)
        smPu=0
        for pu=1:hr
            smPu=smPu+graph[pu]
        end
        smPu=smPu-(graph[hr]*hr)
        push!(ce,smPu)
    end
    wnd.ce=deepcopy(ce)
    wnd.pu=deepcopy(graph)
    return nothing
end
