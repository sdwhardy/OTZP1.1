
################################################################################
#################################### EENS ######################################
function eensF_eqp_eens(eqp, S, ks, wp)
    cpt_tbl=eensF_eqp_cpt(eqp,S)#make capacity probability table
    eens_all=[]#Create eens array
    eens=0.0

    for i=1:length(cpt_tbl[:,1])#loop through rows of cpt
        ratio_curt=cpt_tbl[i,1]/S#find PU curtailment ratio
        diff=wp.pu.-ratio_curt#find closest PU of wind power series to PU curtail ratio
        #i_min=argmin(sqrt.((diff[:]).^2))diff[5240]
        i_min=argmin(abs.((diff[:])))
#i_min=argmin(diff[:])
        if ratio_curt>=1#check if curt ratio is at or above full power and set ce=curtailed energy to 0
            ce=wp.ce[1]
        elseif ratio_curt<=0#check if curt ratio is at or below zero power and set ce to max
            ce=wp.ce[length(wp.ce)]
        elseif i_min < length(diff) && diff[i_min]<0#if curt ratio is a mid point interpolate ce
            ce=eensF_intPole(ratio_curt,wp.pu[i_min],wp.pu[i_min+1],wp.ce[i_min],wp.ce[i_min+1])
        elseif i_min > 1 && diff[i_min]>0
            ce=eensF_intPole(ratio_curt,wp.pu[i_min-1],wp.pu[i_min],wp.ce[i_min-1],wp.ce[i_min])
        else#if exact match occurs
            ce=wp.ce[i_min]
        end

        push!(eens_all, ce*S*cpt_tbl[i,2])#multiply PU curtailed energy with max power and availability, then store
    end
    eens=sum(eens_all)*ks.life*ks.E_op#sum all eens and multiply by cost factors
    return eens
end


#The calculation of equipment level capacity probability table **
function eensF_eqp_cpt(eqp,S)
    #Calculate failure rate for entire length if cable
    if typeof(eqp)==typeof(cbl())
        eqp.reliability.fr=(eqp.reliability.fr/100.0)*eqp.length
    end
    #Calculate Availability of eqiupment
    A_eqp=1.0/(1.0+eqp.reliability.fr*(eqp.reliability.mttr*30.0*24.0/8760.0))
    #Create combinatorial matrix of 0s and 1s
    clms=trunc(Int,eqp.num)
    rows=trunc(Int, 2.0^clms)
    empty_tbl=eensF_blankTbl(rows,clms)
    #Create blank power and availability tables
    PWR_tbl=zeros(Float32,rows,1)
    AVL_tbl=ones(Float32,rows,1)
    #Set powers and availabilities by looping through the CPT
    for k=1:clms
        for j=1:rows
        #if 1 the equipment is functional and the power is added to total
        #the availability is multiplied
            if trunc(Int,empty_tbl[j,k])==1
                AVL_tbl[j]=AVL_tbl[j]*A_eqp
                PWR_tbl[j]=min(S,PWR_tbl[j]+eqp.mva)
        #if 0 the equipment is broken and no power is transmitted
            else
                AVL_tbl[j]=AVL_tbl[j]*(1-A_eqp)
            end
          end
      end
    #all unique power levels are extracted
    tbl_c1=unique(PWR_tbl)
    tbl_c2=zeros(Float32,length(tbl_c1),1)
    for k=1:length(tbl_c1)
      for j=1:length(AVL_tbl)
    #Availabilities are summed for common power levels
          if PWR_tbl[j]==tbl_c1[k]
              tbl_c2[k]=tbl_c2[k]+AVL_tbl[j]
          end
      end
    end
    #Checks if probability sums to 1 else error is thrown
    if sum(tbl_c2) > 1.0001 || sum(tbl_c2) < 0.9999
        println("sum is: "*string(sum(tbl_c2)))
        error("probability does not sum to 1")
    elseif maximum(tbl_c1) > S && minimum(tbl_c1) > 1
        error("power is not correct")
    else
      return [tbl_c1 tbl_c2]
    end
end

#creates a blank capacity probability table **
function eensF_blankTbl(rows,clms)
    XFM_CBL=trunc.(Int8,zeros(rows,clms))
#=create all combinations ie
transpose(
11101000
11010100
10110010
)
=#
    round=1
    k=1
    multi=1
    while round<=clms
      while k<rows
        while k<=(multi*2^(clms-round))
          XFM_CBL[k,round]=1
          k=k+1
        end
        multi=multi+2
        k=k+2^(clms-round)
      end
      round=round+1
      k=1
      multi=1
    end
    return XFM_CBL
end

#linearly interpolates 2 points of graph **
function eensF_intPole(true_x,min_x,max_x,min_y,max_y)
    slope=(max_y-min_y)/(max_x-min_x)
    b=min_y-slope*min_x
    true_y=slope*true_x+b
    return true_y
end
