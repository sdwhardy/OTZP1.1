#=ocn=deepcopy(ocean)
i=46
nsb=bn0[1]
circ=hv_orig=#

#kept each loop
function opt_set()
    return 1
end
#number returned as aft circuits
function opt_setCompPart()
    return 1
end
function opt_compoundOSS_single(ocn,i,max_cst)
    println("rank: "*string(i))
    #find binary and location of 1s of loop index
    i_bn=top_dec2bin(i)
    bn0=findall(x->x==1,i_bn)
    #only process if 3 or more OWPPs are connected
    if (length(bn0)>2)
        #copy originals
        hv_orig=deepcopy(ocn.hv_circuits[i][1])
        mv_orig=deepcopy(ocn.mv_circuits[i][1])
        hv_bsf=deepcopy(ocn.hv_circuits[i][1])
        mv_bsf=deepcopy(ocn.mv_circuits[i][1])
        #println("bn0: "*string(bn0))
        for nsb in bn0
            #println("NSB: "*string(nsb))
            #println("nsb: "*string(nsb))
            nsb_sum_bn,i_min_sum_bn=opt_NSBandi(hv_orig.binary,nsb)
            bn_i=findall(x->x==1,i_min_sum_bn)
            #println(string(top_bin2dec(nsb_sum_bn))*"nsb_sum_bn: "*string(nsb_sum_bn)*" i_min_sum_bn: "*string(i_min_sum_bn)*"- "*string(top_bin2dec(i_min_sum_bn)))
            if (length(bn_i)>1)
                #find mv/hv topology combo
                #ppf_equipment_OSS_MOG(ocean,forwardMV)
                i_min_sum_int=round(Int32,top_bin2dec(i_min_sum_bn))
                #forwardHV, hv_lv=opt_mvhv_layout(deepcopy(hv_orig),nsb_sum_bn,i_min_sum_int,ocn)
                forwardHV=opt_copyBaseEquipment(deepcopy(hv_orig),nsb_sum_bn,ocn)
                #forwardMV, mv_lv=opt_mvhv_layout(deepcopy(mv_orig),nsb_sum_bn,i_min_sum_int,ocn)
                forwardMV=opt_copyBaseEquipment(deepcopy(mv_orig),nsb_sum_bn,ocn)
                #decompose circuits
                hv_bd,ocn.hv_circuits[i]=opt_decomposeLO_single(forwardHV, deepcopy(hv_bsf), i_min_sum_bn, ocn.hv_circuits, ocn,max_cost)
                mv_bd,ocn.mv_circuits[i]=opt_decomposeLO_single(forwardMV, deepcopy(mv_bsf), i_min_sum_bn, ocn.mv_circuits, ocn,max_cost)
                if (hv_bd.cost<hv_bsf.cost)
                    #println("NSB3: "*string(nsb))
                    hv_bsf=deepcopy(hv_bd)
                end
                if (mv_bd.cost<mv_bsf.cost)
                    #println("NSB2: "*string(nsb))
                    mv_bsf=deepcopy(mv_bd)
                end
            end
        end
    end

    return ocn
end

function opt_decomposeLO_single(forward, full, aft_bn, circs, ocn)
    Q_2bd=opt_buildDaBDQ(aft_bn, circs, ocn)
    #println(length(Q_2bd))
    #sort!(Q_2bd, by = v -> v.cost, rev=false)
    #println("Q_2bd: "*string(length(Q_2bd)))
    if (length(Q_2bd)>0)
        for q_set in Q_2bd[1:length(Q_2bd)]
            #println("q owpps: "*string(length(q_set.owpps)))
            #println("f owpps: "*string(length(forward.owpps)))
            single_parent=circuit()
            single_parent=deepcopy(forward)
            single_parent.id=lof_combineIDs(single_parent.id,q_set.id)
            #println("SP cost0: "*string(single_parent.cost))
            pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_copyWindProfiles(ocn.owpps,forward)
            if (length(q_set.owpps)>1)
                single_parent=opt_compoundCbls(q_set,single_parent,ocn)
                pHv,wHv,p132,w132,p220,w220,p400,w400=opt_combineCblPW(single_parent,q_set,pHv,wHv,p132,w132,p220,w220,p400,w400)
            else
                single_parent,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_str8Oss2Oss(q_set,single_parent,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
            end
            single_parent.osss_mog[1]=opt_mogXfmrs(single_parent.osss_mog[1],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,single_parent.pcc_cbls[length(single_parent.pcc_cbls)].elec.volt)
            opt_ttlMvCirc(single_parent)
            #nsb_sum_int=round(Int32,top_bin2dec(nsb_sum_bn))

            single_parent=opt_circOpt(single_parent,ocn)
            #push!(circs[round(Int,single_parent.decimal)],deepcopy(single_parent))
            insert_and_dedup!(circs[round(Int,single_parent.decimal)],deepcopy(single_parent))
            circs[round(Int,single_parent.decimal)]=opt_unicOnly(circs[round(Int,single_parent.decimal)])
            if (single_parent.cost<full.cost)
                println("SP cost1: "*string(single_parent.cost)*" full cost: "*string(full.cost))
                full=deepcopy(single_parent)
            end
            #println("circs: "*string(length(circs[round(Int,single_parent.decimal)])))
        end
    end
    return deepcopy(full), circs[round(Int,full.decimal)]
end
#ocn=ocean
#nsb=bn0[1]
function opt_compoundOSS(ocn)
    #=all_circsMV=Array{Array{circuit, 1}, 1}()
    #all_circsHV=Array{Array{circuit, 1}, 1}()
    for (i,c) in enumerate(ocn.hv_circuits)
        tp=circuit[]
        push!(tp,c)
        push!(all_circsHV,deepcopy(tp))
    end
    for (i,c) in enumerate(ocn.mv_circuits)
        tp=circuit[]
        push!(tp,c)
        push!(all_circsMV,deepcopy(tp))
    end
    mv=-1
    hv=1=#
#i=length(ocn.hv_circuits)
    for i=7:1:length(ocn.hv_circuits)
        println("rank: "*string(i))
        #find binary and location of 1s of loop index
        i_bn=top_dec2bin(i)
        bn0=findall(x->x==1,i_bn)
        #only process if 3 or more OWPPs are connected
        if (length(bn0)>2)
            #copy originals
            hv_orig=deepcopy(ocn.hv_circuits[i][1])
            mv_orig=deepcopy(ocn.mv_circuits[i][1])
            hv_bsf=deepcopy(ocn.hv_circuits[i][1])
            mv_bsf=deepcopy(ocn.mv_circuits[i][1])
            #println("bn0: "*string(bn0))
            for nsb in bn0
                #println("NSB: "*string(nsb))
                #println("nsb: "*string(nsb))
                nsb_sum_bn,i_min_sum_bn=opt_NSBandi(hv_orig.binary,nsb)
                bn_i=findall(x->x==1,i_min_sum_bn)
                #println(string(top_bin2dec(nsb_sum_bn))*"nsb_sum_bn: "*string(nsb_sum_bn)*" i_min_sum_bn: "*string(i_min_sum_bn)*"- "*string(top_bin2dec(i_min_sum_bn)))
                if (length(bn_i)>1)
                    #find mv/hv topology combo
                    #ppf_equipment_OSS_MOG(ocean,forwardMV)
                    i_min_sum_int=round(Int32,top_bin2dec(i_min_sum_bn))
                    #forwardHV, hv_lv=opt_mvhv_layout(deepcopy(hv_orig),nsb_sum_bn,i_min_sum_int,ocn)
                    forwardHV=opt_copyBaseEquipment(deepcopy(hv_orig),nsb_sum_bn,ocn)
                    #forwardMV, mv_lv=opt_mvhv_layout(deepcopy(mv_orig),nsb_sum_bn,i_min_sum_int,ocn)
                    forwardMV=opt_copyBaseEquipment(deepcopy(mv_orig),nsb_sum_bn,ocn)
                    #decompose circuits
                    hv_bd,ocn.hv_circuits[i]=opt_decomposeLO(forwardHV, deepcopy(hv_bsf), i_min_sum_bn, ocn.hv_circuits, ocn)
                    mv_bd,ocn.mv_circuits[i]=opt_decomposeLO(forwardMV, deepcopy(mv_bsf), i_min_sum_bn, ocn.mv_circuits, ocn)
                    if (hv_bd.cost<hv_bsf.cost)
                        #println("NSB3: "*string(nsb))
                        hv_bsf=deepcopy(hv_bd)
                    end
                    if (mv_bd.cost<mv_bsf.cost)
                        #println("NSB2: "*string(nsb))
                        mv_bsf=deepcopy(mv_bd)
                    end
                end
            end
            #optimize circuits
            #=mv_bsf=opt_circOpt(mv_bsf,ocn)
            hv_bsf=opt_circOpt(hv_bsf,ocn)
            #update if cheaper

            if (hv_bsf.cost<ocn.hv_circuits[i].cost)
                println(string(ocn.hv_circuits[i].cost)*" -HV- "*string(hv_bsf.cost))
                ocn.hv_circuits[i]=deepcopy(hv_bsf)
            end

            if (mv_bsf.cost<ocn.mv_circuits[i].cost)
                println(string(ocn.mv_circuits[i].cost)*" -MV- "*string(mv_bsf.cost))
                ocn.mv_circuits[i]=deepcopy(mv_bsf)
            end
            #=if (ocn.mv_circuits[i].oyas.chichi[2]==0 || ocn.mv_circuits[i].oyas.haha[2]==0)
                ocn.mv_circuits[i].oyas=mv_orig.oyas
                println(string(ocn.mv_circuits[i].decimal)*" oya: "*string(ocn.mv_circuits[i].oyas))
            end
            if (ocn.hv_circuits[i].oyas.chichi[2]==0 || ocn.hv_circuits[i].oyas.haha[2]==0)
                ocn.hv_circuits[i].oyas=hv_orig.oyas
                println(string(ocn.hv_circuits[i].decimal)*" oya: "*string(ocn.hv_circuits[i].oyas))
            end=#=#
        end
    end
    return ocn
end
#=
forward=forwardHV
full=deepcopy(hv_bsf)
aft_bn=i_min_sum_bn
circs=ocn.hv_circuits
=#
#decompose selected layout
#q_set=Q_2bd[1]
#ppf_equipment_OSS_MOG(ocean,q_set)
#ppf_equipment_OSS_MOG(ocean,forward)
#ppf_equipment_OSS_MOG(ocean,full)
function opt_decomposeLO(forward, full, aft_bn, circs, ocn)
    Q_2bd=opt_buildDaBDQ(aft_bn, circs, ocn)
    #println(length(Q_2bd))
    #sort!(Q_2bd, by = v -> v.cost, rev=false)
    if (length(Q_2bd)>0)
        for q_set in Q_2bd[1:length(Q_2bd)]
            #println("q owpps: "*string(length(q_set.owpps)))
            #println("f owpps: "*string(length(forward.owpps)))
            single_parent=circuit()
            single_parent=deepcopy(forward)
            single_parent.id=lof_combineIDs(single_parent.id,q_set.id)
            #println("SP cost0: "*string(single_parent.cost))
            pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_copyWindProfiles(ocn.owpps,forward)
            if (length(q_set.owpps)>1)
                single_parent=opt_compoundCbls(q_set,single_parent,ocn)
                pHv,wHv,p132,w132,p220,w220,p400,w400=opt_combineCblPW(single_parent,q_set,pHv,wHv,p132,w132,p220,w220,p400,w400)
            else
                single_parent,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_str8Oss2Oss(q_set,single_parent,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
            end
            single_parent.osss_mog[1]=opt_mogXfmrs(single_parent.osss_mog[1],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,single_parent.pcc_cbls[length(single_parent.pcc_cbls)].elec.volt)
            opt_ttlMvCirc(single_parent)
            #nsb_sum_int=round(Int32,top_bin2dec(nsb_sum_bn))

            single_parent=opt_circOpt(single_parent,ocn)
            #push!(circs[round(Int,single_parent.decimal)],deepcopy(single_parent))
            insert_and_dedup!(circs[round(Int,single_parent.decimal)],deepcopy(single_parent))
            #circs[round(Int,single_parent.decimal)]=opt_unicOnly(circs[round(Int,single_parent.decimal)])
            if (single_parent.cost<full.cost)
                println("SP cost1: "*string(single_parent.cost)*" full cost: "*string(full.cost))
                full=deepcopy(single_parent)
            end
            #println("circs: "*string(length(circs[round(Int,single_parent.decimal)])))
        end
    end
    return deepcopy(full), circs[round(Int,full.decimal)]
end


function opt_buildDaBDQ(aft_bn, circs,ocn)
    #println(string(top_bin2dec(aft_bn))*" :aft_bn: "*string(aft_bn))
    aftbn_1s=findall(x->x==1,aft_bn)
    Q_2bd=Array{Array{Tuple{Int32,Int32},1},1}()
    if (length(aftbn_1s)>1)
        Q_2bd=opt_factorBn(aft_bn,ocn,aftbn_1s)
    end
    return Q_2bd
end

#aft_bn[5]=1
function opt_factorBn(aft_bn,ocn,aftbn_1s)
    tbl_mvCircs, tbl_hvCircs, aftmid_indx=opt_factTble(aft_bn,ocn,aftbn_1s)
    mv_fset,mv_c=opt_rollUp_partial(tbl_mvCircs, aftmid_indx)
    hv_fset,hv_c=opt_rollUp_partial(tbl_hvCircs, aftmid_indx)
    mvhv_c=Array{Array{circuit, 1}, 1}()
    for i=1:length(hv_c)
        push!(mvhv_c,optSysCompare2(mv_c[i],hv_c[i]))
        #println(mv_c[i])
        #println(hv_c[i])
    end
    mvhv_fset,mvhv_c=opt_rollUp_partial(mvhv_c, aftmid_indx)
    #mvhv_fset=comboCompareSing(mvhv_fset)
    toReturnTemp=optSysCompare2(mv_fset,hv_fset)
    toReturnTemp=keep_unic(toReturnTemp)
    toReturn=optSysCompare2(toReturnTemp,mvhv_fset)
    toReturn=keep_unic(toReturn)
    return toReturn
end

function opt_factTble(aft_bn,ocn,aftbn_1s)
    base_Int=round(Int32,2^(aftbn_1s[1]-1))
    aftbn_0s=findall(x->x==0,aft_bn)
    aftbn_dec=round(Int32,top_bin2dec(aft_bn))
    base_Int_comp=aftbn_dec-base_Int
    aftmid_dec=round((aftbn_dec/2))
    tbl_bn=top_hvTopos(aftbn_1s[length(aftbn_1s)])
    aftmid_indx=0
    set=Bool
    set=false
    tbl_bnR=Array{Array{Int8,1},1}()
    tbl_bnMns=Array{Array{Int8,1},1}()
    tbl_Ints=Array{Int32,1}()
    tbl_mvCircs=Array{Array{circuit, 1}, 1}()
    tbl_hvCircs=Array{Array{circuit, 1}, 1}()

    for row_bn=length(tbl_bn[:,1]):-1:1
        push!(tbl_bnR,reverse!(tbl_bn[row_bn,:]))
    end

    for (row_i,row_bn) in enumerate(tbl_bnR)
        aftbn_out1s=findall(x->x==1,row_bn)
        row_Int=round(Int32,top_bin2dec(row_bn))
        keep=true
        if (row_Int==base_Int || row_Int==base_Int_comp)
            keep=false
        end
        for po1 in aftbn_out1s
            for po0 in aftbn_0s
                if (po1==po0)
                    keep=false
                    break
                end
            end
            if (keep==false)
                break
            end
        end
        if (keep==true)
            #push!(tbl_bnMns,deepcopy(row_bn))

            if (row_Int>aftmid_dec && set==false)
                aftmid_indx=deepcopy(length(tbl_mvCircs)+1)
                set=true
            end
            #push!(tbl_Ints,deepcopy(row_Int))
            push!(tbl_mvCircs,deepcopy(ocn.mv_circuits[row_Int]))
            push!(tbl_hvCircs,deepcopy(ocn.hv_circuits[row_Int]))
        end
    end
    return tbl_mvCircs, tbl_hvCircs, aftmid_indx
end

#=circs=tbl_mvCircs
indx0=1
crc0=circs[indx0]
crc1=crc0[1]
crc2=circs[indx0+2]
crc3=crc2[1]
rnk_i=5=#
function keep_unic(bsf)
    unic=circuit[]
    eps=1
    for (i0,c) in enumerate(bsf[1:length(bsf)-1])
        crcCopy=false
        for (i1,d) in enumerate(bsf[i0+1:length(bsf)])
            if (c.cost-eps<d.cost && c.cost+eps>d.cost)
                if (length(bsf[i0].osss_mog)==length(d.osss_mog) && length(bsf[i0].osss_owp)==length(d.osss_owp) && length(bsf[i0].owp_MVcbls)==length(d.owp_MVcbls) && length(bsf[i0].owp_HVcbls)==length(d.owp_HVcbls) && length(bsf[i0].oss2oss_cbls)==length(d.oss2oss_cbls) && length(bsf[i0].pcc_cbls)==length(d.pcc_cbls))
                    num_array0=id_breakDown_Check(deepcopy(bsf[i0].id))
                    num_array1=id_breakDown_Check(deepcopy(d.id))
                    #println(string(num_array0)*" - "*string(num_array1))
                    if (num_array0==num_array1)
                        #println("the same!")
                        crcCopy=true
                        break
                    end
                end
            end
        end
        if (crcCopy==false)
            #println("stored!")
            push!(unic,deepcopy(c))
        end
    end
    if (length(unic)==0)
        push!(unic,bsf[1])
    end
    return unic
end

function id_breakDown_Check(eyeD)
    tst1=split(eyeD, "_")
    num_array=Int64[]
    for tst in tst1
        tInt=split(tst, "#")
        #println(tInt[1])
        push!(num_array,parse(Int,tInt[1]))
    end
    a=sort!(num_array)
    return a
end

function opt_rollUp_partial(circs, aftmid_indx)
    complete_systems=circuit[]
    for crc in circs[length(circs)]
        #push!(complete_systems, deepcopy(crc))
        insert_and_dedup!(complete_systems, deepcopy(crc))
        complete_systems=opt_unicOnly_CompPart(complete_systems)
    end
    base_1s=findall(x->x==1,circs[length(circs)][1].binary)
    for (indx0,crc0) in enumerate(circs[1:aftmid_indx])
        bn0=findall(x->x==1,crc0[length(crc0)].binary)
        for (indx1,crc1) in enumerate(crc0)
            for crc2 in circs[indx0+1:length(circs)]

                bn1=findall(x->x==1,crc2[length(crc2)].binary)
                cmb, bn01=opt_willCombine(bn0,bn1,length(crc2[length(crc2)].binary))
                if (cmb)
                    for (indx3,crc3) in enumerate(crc2)

                        #println(string(cmb)*" bn01: "*string(bn01)*" bn0: "*string(bn0)*" bn1: "*string(bn1))
                        rnk_i=0
                        bn01_1s=findall(x->x==1,bn01)
                        rnk=floor(Int64,top_bin2dec(bn01))
                        for (i_cp,c_pf) in enumerate(circs)
                            if (c_pf[1].decimal==rnk)
                                rnk_i=deepcopy(i_cp)
                                break
                            end
                        end
                        if (rnk_i != 0)
                            if (length(bn01_1s)==length(base_1s))
                                #push!(complete_systems, deepcopy(opt_replaceRnk(circs[rnk_i][1],crc1,crc3)))
                                insert_and_dedup!(complete_systems, deepcopy(opt_replaceRnk(circs[rnk_i][1],crc1,crc3)))
                                complete_systems=opt_unicOnly_CompPart(complete_systems)
                                #=if (length(complete_systems)>2*opt_set())
                                    complete_systems=complete_systems[1:2*opt_set()]
                                end=#

                            else
                                #push!(circs[rnk_i], deepcopy(opt_replaceRnk(circs[rnk_i][1],crc1,crc3)))
                                insert_and_dedup!(circs[rnk_i], deepcopy(opt_replaceRnk(circs[rnk_i][1],crc1,crc3)))
                                circs[rnk_i]=opt_unicOnly(circs[rnk_i])
                                #=if (length(circs[rnk_i])>2*opt_set())
                                    circs[rnk_i]=circs[rnk_i][1:2*opt_set()]
                                end=#

                            end
                        end
                    end
                end

            end
        end
        #println("Circuit Number: "*string(indx0))
    end
    #=sort!(complete_systems, by = v -> v.cost, rev=false)
    complete_systems=opt_unicOnly(complete_systems)
    #=if (length(complete_systems)>10)
        complete_systems=complete_systems[1:10]
    end=#
    for (i,c) in enumerate(circs)
        sort!(circs[i], by = v -> v.cost, rev=false)
        circs[i]=opt_unicOnly(circs[i])
    end=#
    #=for (i,c) in enumerate(circs)
        if (length(c)>10)
            circs[i]=c[1:10]
        else
        end
    end=#
    complete_systems=keep_unic(complete_systems)
    return complete_systems, circs
end

#circs=hv_fset
#c=circs[1]
function opt_unicOnly_comSys(circs)
    if (length(circs)>1)
        i=1
        while (i < length(circs))
            str0=deepcopy(circs[i].id)
            str1=deepcopy(circs[i+1].id)
            tst0=split(str0, "_")
            tst1=split(str1, "_")
            a=sort!(tst0)
            b=sort!(tst1)
            i=i+1
            if (a==b)
                deleteat!(circs, i)
                i=i-1
            end
        end
        #circs=keep_unic(circs)
    else
    end
    return circs
end

function opt_unicOnly(circs)
    if (length(circs)>1)
        #=i=1
        while (i < length(circs))
            str0=deepcopy(circs[i].id)
            str1=deepcopy(circs[i+1].id)
            tst0=split(str0, "_")
            tst1=split(str1, "_")
            a=sort!(tst0)
            b=sort!(tst1)
            i=i+1
            if (a==b)
                deleteat!(circs, i)
                i=i-1
            end
            if (i==opt_set())
                break
            end
        end=#
        circs=keep_unic(circs)
        if (length(circs)>opt_set())
            lm=opt_set()
        else
            lm=length(circs)
        end
        circs=circs[1:lm]
    else
    end
    return circs
end


function opt_unicOnly_CompPart(circs)
    if (length(circs)>1)
        i=1
        while (i < length(circs))
            str0=deepcopy(circs[i].id)
            str1=deepcopy(circs[i+1].id)
            tst0=split(str0, "_")
            tst1=split(str1, "_")
            a=sort!(tst0)
            b=sort!(tst1)
            i=i+1
            if (a==b)
                deleteat!(circs, i)
                i=i-1
            end
            if (i==opt_setCompPart())
                break
            end
        end
        #circs=keep_unic(circs)
        if (length(circs)>opt_setCompPart())
            lm=opt_setCompPart()
        else
            lm=length(circs)
        end
        circs=circs[1:lm]
    else
    end
    return circs
end


#=
function opt_decomposeLO_old(forward, full, aft_bn, circs, ocn)
    Q_2bd=opt_buildDaBDQ(aft_bn, circs, ocn)
    if (length(Q_2bd)>0)
        println("length of Q_2bd: "*string(length(Q_2bd)))
        println(string(top_bin2dec(aft_bn))*" - "*string(aft_bn)*" - Broken down into: ")
        for q_set in Q_2bd
            println(string(q_set))
            single_parent=deepcopy(forward)
            pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_copyWindProfiles(ocn.owpps,forward)
            for q in q_set
                orphan,parnt=lof_oyaCirc(q[1],q[2],ocn)
                if (orphan==false)
                    if (length(parnt.owpps)>1)
                        single_parent=opt_compoundCbl(parnt,single_parent,ocn)
                        pHv,wHv,p132,w132,p220,w220,p400,w400=opt_combineCblPW(single_parent,parnt,pHv,wHv,p132,w132,p220,w220,p400,w400)
                    else
                        single_parent,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_str8Oss2Oss(parnt,single_parent,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
                    end
                else
                    println("Error orphan should not exist here!!!")
                end
            end
            single_parent.osss_mog[1]=opt_mogXfmrs(single_parent.osss_mog[1],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,single_parent.pcc_cbls[length(single_parent.pcc_cbls)].elec.volt)
            opt_ttlMvCirc(single_parent)
            #nsb_sum_int=round(Int32,top_bin2dec(nsb_sum_bn))
            if (single_parent.cost<full.cost)
                single_parent.oyas=full.oyas
                #println(string(single_parent.cost)*" - "*string(full.cost)*" - "*string(single_parent.oyas))
                full=deepcopy(single_parent)
            end
        end
    end
    return full
end
#build Break-down Queue
function opt_buildDaBDQ_old(aft_bn, circs)
    println(string(top_bin2dec(aft_bn))*" :aft_bn: "*string(aft_bn))
    aftbn_1s=findall(x->x==1,aft_bn)
    Q_2bd=Array{Array{Tuple{Int32,Int32},1},1}()
    #println("aftbn_1s: "*string(aftbn_1s))
    if (length(aftbn_1s)>2)
        #println("length(aftbn_1s)>2")
        aft_int=round(Int32,top_bin2dec(aft_bn))
        println(circs[aft_int].oyas)
        if ((circs[aft_int].oyas.chichi[2] != 0) && (circs[aft_int].oyas.haha[2] != 0))

            #initialize Q_2bd
            push!(Q_2bd,[circs[aft_int].oyas.chichi,circs[aft_int].oyas.haha])
            #println("Q_2bd length: "*string(length(Q_2bd)))
            q_set_i=1
            while (q_set_i <= length(Q_2bd))
                #rintln("length: "*string(length(Q_2bd)))
                q_set=Q_2bd[q_set_i]
                #println("q_set: "*string(q_set))
                for (iq,q) in enumerate(q_set)
                    q_bn=top_dec2bin(q[1])
                    qbn_1s=findall(x->x==1,q_bn)
                    #println("qbn_1s: "*string(qbn_1s))
                    if (length(qbn_1s)>2)
                        #println("(length(qbn_1s)>2)")
                        if ((circs[q[1]].oyas.chichi[2] != 0) && (circs[q[1]].oyas.haha[2] != 0))
                            #println("((circs[q[1]].oyas.chichi[2] != 0) && (circs[q[1]].oyas.haha[2] != 0))")
                            q_new=Array{Tuple{Int32,Int32},1}()
                            for (iq_not_q,q_not_q) in enumerate(q_set)
                                if (iq_not_q!=iq)
                                    push!(q_new,deepcopy(q_not_q))
                                end
                            end
                            push!(q_new,circs[q[1]].oyas.chichi)
                            push!(q_new,circs[q[1]].oyas.haha)
                            push!(Q_2bd,deepcopy(q_new))
                        end
                    end
                end
                q_set_i=q_set_i+1
            end#while
        end#if
    end#if
    return Q_2bd
end
=#
function lof_oyaCirc(chichi,kv,ocn)
    orphan=false
    chichiBD=circuit()
    if (kv>0)
        #HV chichi
        chichiBD=deepcopy(ocn.hv_circuits[chichi])
    elseif (kv<0)
        #MV chichi
        chichiBD=deepcopy(ocn.mv_circuits[chichi])
    else
        #no chichi
        println("No parent component found for: "*string(chichi))
        orphan=true
    end
    return orphan, chichiBD
end

function lof_oyasCircs(q1_int, circs, ocn)
    orphan=false
    chichiBD=circuit()
    hahaBD=circuit()
    if (circs[q1_int].oyas.chichi[2]>0)
        #HV chichi
        chichiBD=deepcopy(ocn.hv_circuits[circs[q1_int].oyas.chichi[1]])
    elseif (circs[q1_int].oyas.chichi[2]<0)
        #MV chichi
        chichiBD=deepcopy(ocn.mv_circuits[circs[q1_int].oyas.chichi[1]])
    else
        #no chichi
        println("No chichi component found for: "*string(q1_int))
        orphan=true
    end
    if (circs[q1_int].oyas.haha[2]>0)
        #HV haha
        hahaBD=deepcopy(ocn.hv_circuits[circs[q1_int].oyas.haha[1]])
    elseif (circs[q1_int].oyas.haha[2]<0)
        #MV chichi
        hahaBD=deepcopy(ocn.mv_circuits[circs[q1_int].oyas.haha[1]])
    else
        #no haha
        println("No haha component found for: "*string(q1_int))
        orphan=true
    end
    return orphan, chichiBD, hahaBD
end
#find mv/hv topology combo
function opt_mvhv_layout(circ,nsb_sum_bn,i_min_sum_int,ocn)
    #take base information
    forward_base=opt_copyBaseEquipment(circ,nsb_sum_bn,ocn)
    aft_mv=deepcopy(ocn.mv_circuits[i_min_sum_int][1])
    aft_hv=deepcopy(ocn.hv_circuits[i_min_sum_int][1])
    cmp_mv=opt_compoundCbl(aft_mv,deepcopy(forward_base),ocn)
    cmp_hv=opt_compoundCbl(aft_hv,deepcopy(forward_base),ocn)
    pMvb,pHvb,p132b,p220b,p400b,wMvb,wHvb,w132b,w220b,w400b=opt_copyWindProfiles(ocn.owpps,forward_base)
    pHvmv,wHvmv,p132mv,w132mv,p220mv,w220mv,p400mv,w400mv=opt_combineCblPW(cmp_mv,aft_mv,deepcopy(pHvb),deepcopy(wHvb),deepcopy(p132b),deepcopy(w132b),deepcopy(p220b),deepcopy(w220b),deepcopy(p400b),deepcopy(w400b))
    pHvhv,wHvhv,p132hv,w132hv,p220hv,w220hv,p400hv,w400hv=opt_combineCblPW(cmp_hv,aft_hv,deepcopy(pHvb),deepcopy(wHvb),deepcopy(p132b),deepcopy(w132b),deepcopy(p220b),deepcopy(w220b),deepcopy(p400b),deepcopy(w400b))
    cmp_mv.osss_mog[1]=opt_mogXfmrs(cmp_mv.osss_mog[1],pMvb,pHvmv,p132mv,p220mv,p400mv,wMvb,wHvmv,w132mv,w220mv,w400mv,ocn.finance,cmp_mv.pcc_cbls[length(cmp_mv.pcc_cbls)].elec.volt)
    cmp_hv.osss_mog[1]=opt_mogXfmrs(cmp_hv.osss_mog[1],pMvb,pHvhv,p132hv,p220hv,p400hv,wMvb,wHvhv,w132hv,w220hv,w400hv,ocn.finance,cmp_hv.pcc_cbls[length(cmp_hv.pcc_cbls)].elec.volt)
    opt_ttlMvCirc(cmp_mv)
    opt_ttlMvCirc(cmp_hv)
    nsb_sum_int=round(Int32,top_bin2dec(nsb_sum_bn))
    #cmp_mv.oyas.chichi=(nsb_sum_int,mhv)
    #cmp_mv.oyas.haha=(i_min_sum_int,-1)
    #cmp_hv.oyas.chichi=(nsb_sum_int,mhv)
    #cmp_hv.oyas.haha=(i_min_sum_int,1)
    #cmp_hv=opt_circOpt(cmp_hv,ocn)
    #cmp_mv=opt_circOpt(cmp_mv,ocn)
    if (cmp_mv.cost<cmp_hv.cost)
        cmp=cmp_mv
        #aft=aft_mv
    else
        cmp=cmp_hv
        #aft=aft_hv
    end
    return forward_base
end

#optimize an individual circuit
function opt_circOpt(crc,ocn)
    longest_cable=lof_pnt2pnt_dist(ocn.owpps[length(ocn.owpps)].node.xy,ocn.pccs[length(ocn.pccs)].node.xy)
    c_o=deepcopy(crc)
    #println(string(c_o.decimal)*") mvC:"*string(length(c_o.owp_MVcbls))*", hvC:"*string(length(c_o.owp_HVcbls))*", oss:"*string(length(c_o.osss_owp))*", mog:"*string(length(c_o.osss_mog))*", o2oC:"*string(length(c_o.oss2oss_cbls)))
    new_coords=opt_reAdjust_oss(crc,ocn.owpps[1].mv_zone,ocn.sys.mvCl,10e-7)
    crc=opt_reAdjust_cbls(crc,new_coords,ocn,longest_cable)
    if ((c_o.cost)<crc.cost)
        #println("Error: re-adjustment failed, circuit: ")
        #print("Initial: "*string(c_o.cost))
        #print(" - Adjusted: "*string(crc.cost))
        new_coords=opt_reAdjust_oss(c_o,ocn.owpps[1].mv_zone,ocn.sys.mvCl,10e-6)
        crc=opt_reAdjust_cbls(crc,new_coords,ocn,longest_cable)
        #print(" - re-adjusted: "*string(crc.cost))
        #println()
        #println(string(c_o.decimal)*") mvC:"*string(length(c_o.owp_MVcbls))*", hvC:"*string(length(c_o.owp_HVcbls))*", oss:"*string(length(c_o.osss_owp))*", mog:"*string(length(c_o.osss_mog))*", o2oC:"*string(length(c_o.oss2oss_cbls)))
        if ((c_o.cost)<crc.cost)
            crc=deepcopy(c_o)
            #println("kept original layout!")
        end
    end
    return crc
end

function opt_NSBandi(binry,nsb)
    nsb_sum_bn=Int8[]
    i_nsb_sum_bn=Int8[]
    for (posi,bt) in enumerate(binry)
        if (posi <= nsb)
            push!(nsb_sum_bn,bt)
            push!(i_nsb_sum_bn,0)
        elseif (posi > nsb)
            push!(nsb_sum_bn,0)
            push!(i_nsb_sum_bn,bt)
        end
    end
    return nsb_sum_bn,i_nsb_sum_bn
end

function opt_copyBaseEquipment(crc,base_bin,ocn)
    oss_system=circuit()
    #base_bin=top_dec2bin(base_dec)
    #base quantities
    oss_system.binary=crc.binary
    oss_system.decimal=crc.decimal
    oss_system.owpps=crc.owpps
    oss_system.base_owp=crc.base_owp
    oss_system.oss_wind=crc.oss_wind
    oss_system.oss_mva=crc.oss_mva
    oss_system.id=crc.id

    #Take base OWP,PCC connections and MOG without transformers
    push!(oss_system.osss_mog,deepcopy(crc.osss_mog[1]))
    oss_system.osss_mog[1].xfmrs=xfo[]
    push!(oss_system.pcc_cbls,deepcopy(crc.pcc_cbls[1]))
    mog=oss_system.osss_mog[1].node.num
    oss_mog=Int32[]
    for (owp_i,bn) in enumerate(base_bin)
        if (bn==1)
            owp=ocn.owpps[owp_i].node.num
            for mv_c in crc.owp_MVcbls
                if (mv_c.pth[1].num==owp)
                    push!(oss_system.owp_MVcbls,deepcopy(mv_c))
                    if (mv_c.pth[length(mv_c.pth)].num!=mog)
                        #copy numbers of nodes at end of MV that aren't the terminous MOG
                        push!(oss_mog,deepcopy(mv_c.pth[length(mv_c.pth)].num))
                    end
                end
            end
        end
    end

    for hv_c in crc.owp_HVcbls
        for oss_num in oss_mog
            if (hv_c.pth[1].num==oss_num)
                push!(oss_system.owp_HVcbls,deepcopy(hv_c))
                if (hv_c.pth[length(hv_c.pth)].num != mog)
                    push!(oss_mog,deepcopy(hv_c.pth[length(hv_c.pth)].num))
                end
            end
        end
    end

    for mvhv_oss in crc.osss_owp
        for oss_num in oss_mog
            if (mvhv_oss.node.num==oss_num)
                push!(oss_system.osss_owp,deepcopy(mvhv_oss))
            end
        end
    end

    for mvhv_mog in crc.osss_mog[2:length(crc.osss_mog)]
        for oss_num in oss_mog
            if (mvhv_mog.node.num==oss_num)
                push!(oss_system.osss_mog,deepcopy(mvhv_mog))
            end
        end
    end
    for o2o_c in crc.oss2oss_cbls
        for oss_num in oss_mog
            if (o2o_c.pth[1].num==oss_num)
                push!(oss_system.oss2oss_cbls,deepcopy(o2o_c))
            end
        end
    end

    return oss_system
end

function opt_copyWindProfiles(owps,oss_system)
    pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_InitPW()
    oss_mog=Int32[]
    for c_mv in oss_system.owp_MVcbls
        if (c_mv.pth[length(c_mv.pth)].num == oss_system.osss_mog[1].node.num)
            push!(pMv,owps[c_mv.pth[1].num].mva)
            push!(wMv,owps[c_mv.pth[1].num].wnd)
        else
            push!(oss_mog,deepcopy(c_mv.pth[1].num))
        end
    end

    for (hv_i,c_hv) in enumerate(oss_system.owp_HVcbls)
        if (c_hv.elec.volt==oss_system.pcc_cbls[1].elec.volt)# this is not updated yet!!!!!!
            push!(pHv,owps[oss_mog[hv_i]].mva)
            push!(wHv,owps[oss_mog[hv_i]].wnd)
        elseif (c_hv.elec.volt==132)
            push!(p132,owps[oss_mog[hv_i]].mva)
            push!(w132,owps[oss_mog[hv_i]].wnd)
        elseif (c_hv.elec.volt==220)
            push!(p220,owps[oss_mog[hv_i]].mva)
            push!(w220,owps[oss_mog[hv_i]].wnd)
        elseif (c_hv.elec.volt==400)
            push!(p400,owps[oss_mog[hv_i]].mva)
            push!(w400,owps[oss_mog[hv_i]].wnd)
        else
        end
    end
    return pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400
end
#circ=q_set
#oss_system=single_parent

function opt_compoundCbls(circ,oss_system,ocn)

    for cb in circ.owp_MVcbls
        push!(oss_system.owp_MVcbls,cb)
    end

    for cb in circ.owp_HVcbls
        push!(oss_system.owp_HVcbls,cb)
    end

    for cb in circ.oss2oss_cbls
        push!(oss_system.oss2oss_cbls,cb)
    end

    for xm in circ.osss_owp
        push!(oss_system.osss_owp,xm)
    end

    for xm in circ.osss_mog
        push!(oss_system.osss_mog,xm)
    end

    #find oss to mog connections
    for pcc_c in circ.pcc_cbls
        con_bus=bus()
        for os in circ.osss_owp
            if (os.node.num==pcc_c.pth[1].num)
                con_bus=deepcopy(os)
                break
            end
        end
        for os in circ.osss_mog
            if (os.node.num==pcc_c.pth[1].num)
                con_bus=deepcopy(os)
                break
            end
        end
        pth=deepcopy(lof_pnt2pnt_dist(con_bus.node.xy,oss_system.osss_mog[1].node.xy))
        #println(string(con_bus.node.num)*" con_bus.mva: "*string(con_bus.mva)*" size: "*string(pcc_c.size)*" num: "*string(pcc_c.num))
        oss2ossCbl=cstF_nextSizeDown(pth,con_bus.mva,pcc_c.elec.volt,con_bus.wnd,ocn.finance,pcc_c.size,pcc_c.num,ocn.eqp_data)
        push!(oss2ossCbl.pth,con_bus.node)
        push!(oss2ossCbl.pth,oss_system.osss_mog[1].node)
        push!(oss_system.oss2oss_cbls,oss2ossCbl)
    end
    return oss_system
end


function opt_compoundCbl(circ,oss_system,ocn)

    for cb in circ.owp_MVcbls
        push!(oss_system.owp_MVcbls,cb)
    end

    for cb in circ.owp_HVcbls
        push!(oss_system.owp_HVcbls,cb)
    end

    for cb in circ.oss2oss_cbls
        push!(oss_system.oss2oss_cbls,cb)
    end

    for xm in circ.osss_owp
        push!(oss_system.osss_owp,xm)
    end

    for xm in circ.osss_mog
        push!(oss_system.osss_mog,xm)
    end

    #find oss to mog connections
    pth=deepcopy(lof_pnt2pnt_dist(circ.osss_mog[1].node.xy,oss_system.osss_mog[1].node.xy))
    oss2ossCbl=cstF_nextSizeDown(pth,circ.oss_mva,circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt,circ.oss_wind,ocn.finance,circ.pcc_cbls[length(circ.pcc_cbls)].size,circ.pcc_cbls[length(circ.pcc_cbls)].num,ocn.eqp_data)
    push!(oss2ossCbl.pth,circ.osss_mog[1].node)
    push!(oss2ossCbl.pth,oss_system.osss_mog[1].node)
    push!(oss_system.oss2oss_cbls,oss2ossCbl)
    return oss_system
end


function optSysCompare2(a1,a2)
    for (i,cs) in enumerate(a2)
        insert_and_dedup!(a1, deepcopy(cs))
        a1=opt_unicOnly_CompPart(a1)
    end

    return deepcopy(a1)
end


function optSysCompare(MVC,HVC_orig)
    HVC=deepcopy(HVC_orig)
    for (i0,hc) in enumerate(HVC)
        for (i,mc) in enumerate(MVC[i0])
            insert_and_dedup!(HVC[i0], deepcopy(mc))
            HVC[i0]=opt_unicOnly(HVC[i0])
            #=if (length(HVC[i])>opt_set())
                lm2=opt_set()
            else
                lm2=length(HVC[i])
            end
            HVC[i]=HVC[i][1:lm2]=#
        end
    end
    return HVC
end

insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x,by = v -> v.cost), [x]); v)

#circs=big_setsHV
function opt_rollUp(ocn, circs)
    complete_systems=circuit[]
    for c in circs[length(circs)]
        insert_and_dedup!(complete_systems, deepcopy(c))
        complete_systems=opt_unicOnly_comSys(complete_systems)
    end
    for (indx0,crc0) in enumerate(circs[1:floor(Int64,length(circs)/2)])
        lm0=0
        #=if (length(crc0)>opt_set())
            lm0=opt_set()
        else
            lm0=length(crc0)
        end=#
        for (indx1,crc1) in enumerate(crc0[1:length(crc0)])
            bn0=findall(x->x==1,crc1.binary)
            for (indx2,crc2) in enumerate(circs[indx0+1:length(circs)])
                #=lm1=0
                if (length(crc2)>opt_set())
                    lm1=opt_set()
                else
                    lm1=length(crc2)
                end=#
                for (indx3,crc3) in enumerate(crc2[1:length(crc2)])
                    bn1=findall(x->x==1,crc3.binary)
                    cmb, bn01=opt_willCombine(bn0,bn1,length(crc1.binary))
                    if (cmb)
                        rnk=floor(Int64,top_bin2dec(bn01))
                        if (rnk==length(circs))
                            insert_and_dedup!(complete_systems, deepcopy(opt_replaceRnk(circs[rnk][1],crc1,crc3)))
                            complete_systems=opt_unicOnly_comSys(complete_systems)
                        else
                            insert_and_dedup!(circs[rnk], deepcopy(opt_replaceRnk(circs[rnk][1],crc1,crc3)))
                            circs[rnk]=opt_unicOnly(circs[rnk])
                        end
                    end
                end
            end
            end
        println("Circuit Number: "*string(indx0))
    end
    complete_systems=keep_unic(complete_systems)
    #complete_systems=deepcopy(opt_sortFllSys(complete_systems))
    return complete_systems, circs
end

function opt_rollUp_firstInLine(ocn, circs)
    complete_systems=circuit[]
    for c in circs[length(circs)]
        insert_and_dedup!(complete_systems, deepcopy(c))
        complete_systems=opt_unicOnly_comSys(complete_systems)
    end
    for (indx0,crc0) in enumerate(circs[1:floor(Int64,length(circs)/2)])
        lm0=0
        #=if (length(crc0)>opt_set())
            lm0=opt_set()
        else
            lm0=length(crc0)
        end=#
        for (indx1,crc1) in enumerate(crc0[1:1])
            bn0=findall(x->x==1,crc1.binary)
            for (indx2,crc2) in enumerate(circs[indx0+1:length(circs)])
                #=lm1=0
                if (length(crc2)>opt_set())
                    lm1=opt_set()
                else
                    lm1=length(crc2)
                end=#
                for (indx3,crc3) in enumerate(crc2[1:1])
                    bn1=findall(x->x==1,crc3.binary)
                    cmb, bn01=opt_willCombine(bn0,bn1,length(crc1.binary))
                    if (cmb)
                        rnk=floor(Int64,top_bin2dec(bn01))
                        if (rnk==length(circs))
                            insert_and_dedup!(complete_systems, deepcopy(opt_replaceRnk(circs[rnk][1],crc1,crc3)))
                            complete_systems=opt_unicOnly_comSys(complete_systems)
                        else
                            insert_and_dedup!(circs[rnk], deepcopy(opt_replaceRnk(circs[rnk][1],crc1,crc3)))
                            circs[rnk]=opt_unicOnly(circs[rnk])
                        end
                    end
                end
            end
            end
        println("Circuit Number: "*string(indx0))
    end
    #complete_systems=deepcopy(opt_sortFllSys(complete_systems))
    return complete_systems, circs
end

function combineAndrank(bfsmv,bfshv,bfsmhv)
    for c in bfsmv
        insert_and_dedup!(bfsmhv,deepcopy(c))
    end
    for c in bfshv
        insert_and_dedup!(bfsmhv,deepcopy(c))
    end
    bfsmhv=keep_unic(bfsmhv)
    return bfsmhv
end

function combineAndrank_old(bfsmv,bfshv,bfsmhv)
    bfsmhv=comboCompare(bfsmhv,bfshv)
    bfsmhv=comboCompare(bfsmhv,bfsmv)
    bfsmhv=deepcopy(opt_sortFllSys(bfsmhv))
    return bfsmhv
end

function check_ids(bfsmv,bfshv)
    eps=0.00001
    uni=circuit[]
    for (cv_i,cv) in enumerate(bfsmv)
        if (cv[1].cost<bfshv[cv_i][1].cost+eps && cv[1].cost>bfshv[cv_i][1].cost-eps)
            if (length(cv[1].osss_mog)==length(bfshv[cv_i][1].osss_mog) && length(cv[1].osss_owp)==length(bfshv[cv_i][1].osss_owp) && length(cv[1].owp_MVcbls)==length(bfshv[cv_i][1].owp_MVcbls) && length(cv[1].owp_HVcbls)==length(bfshv[cv_i][1].owp_HVcbls) && length(cv[1].oss2oss_cbls)==length(bfshv[cv_i][1].oss2oss_cbls) && length(cv[1].pcc_cbls)==length(bfshv[cv_i][1].pcc_cbls))
                bfshv[cv_i][1].id=string(cv_i)*"#mh"
                bfsmv[cv_i][1].id=string(cv_i)*"#mh"
            end
        end
    end
    return bfsmv,bfshv
end


function combineAndrank_id(bfsmv,bfshv,bfsmhv)
    bfsmhv=comboCompare_id(bfsmhv,bfshv)
    bfsmhv=comboCompare_id(bfsmhv,bfsmv)
    bfsmhv=deepcopy(opt_sortFllSys(bfsmhv))
    return bfsmhv
end




function comboCompareSing(scs)
    eps=0.00001
    uni=circuit[]
    for (cv_i,cv) in enumerate(scs)
        crcCopy=false
        for sc in scs[cv_i+1:length(scs)]
            if (cv.cost<sc.cost+eps && cv.cost>sc.cost-eps)
                if (length(cv.osss_mog)==length(sc.osss_mog) && length(cv.osss_owp)==length(sc.osss_owp) && length(cv.owp_MVcbls)==length(sc.owp_MVcbls) && length(cv.owp_HVcbls)==length(sc.owp_HVcbls) && length(cv.oss2oss_cbls)==length(sc.oss2oss_cbls) && length(cv.pcc_cbls)==length(sc.pcc_cbls))
                    crcCopy=true
                    break
                end
            end
        end
        if (crcCopy==false)
            push!(uni,deepcopy(cv))
        end
    end
    return uni
end

function comboCompare(scs,bfs)
    eps=0.5
    for cv in bfs
        crcCopy=false
        for sc in scs
            if (cv.cost<sc.cost+eps && cv.cost>sc.cost-eps)
                if (length(cv.osss_mog)==length(sc.osss_mog) && length(cv.osss_owp)==length(sc.osss_owp) && length(cv.owp_MVcbls)==length(sc.owp_MVcbls) && length(cv.owp_HVcbls)==length(sc.owp_HVcbls) && length(cv.oss2oss_cbls)==length(sc.oss2oss_cbls) && length(cv.pcc_cbls)==length(sc.pcc_cbls))
                    crcCopy=true
                end
            end
        end
        if (crcCopy==false)
            push!(scs,deepcopy(cv))
        end
    end
    return scs
end


function comboCompare_id(scs,bfs)
    for cv in bfs
        crcCopy=false
        for sc in scs
            if (cv.id == sc.id)
                crcCopy=true
            end
        end
        if (crcCopy==false)
            push!(scs,deepcopy(cv))
        end
    end
    return scs
end

#**
function opt_sortFllSys(complete_systems)
    c_is=Array{Tuple{Float64,Int64},1}()
    sorted_circs=circuit[]
    for (i,sys) in enumerate(complete_systems)
        push!(c_is,deepcopy((sys.cost,i)))
    end
    sort!(c_is, by = x -> x[1])
    for cst in c_is
        push!(sorted_circs,deepcopy(complete_systems[cst[2]]))
    end
    return sorted_circs
end

#l=length(circs[indx1].binary)
function opt_willCombine(bn0,bn1,l)
    cmb=true
    bn01=zeros(Int8,l)
    for one0 in bn0
        for one1 in bn1
            if (one0==one1)
                cmb=false
            end
        end
    end
    if (cmb==true)
        for one0 in bn0
            bn01[one0]=1
        end
        for one1 in bn1
            bn01[one1]=1
        end
    end
    return cmb, bn01
end


function opt_replaceRnk(circ,crc0,crc1)
    crc01=deepcopy(circ)
    crc01.osss_owp=deepcopy(crc0.osss_owp)
    crc01.osss_mog=deepcopy(crc0.osss_mog)
    crc01.owp_MVcbls=deepcopy(crc0.owp_MVcbls)
    crc01.owp_HVcbls=deepcopy(crc0.owp_HVcbls)
    crc01.oss2oss_cbls=deepcopy(crc0.oss2oss_cbls)
    crc01.pcc_cbls=deepcopy(crc0.pcc_cbls)
    crc01.cost=deepcopy(crc0.cost+crc1.cost)
    crc01.cost=deepcopy(crc0.cost+crc1.cost)
    crc01.id=crc0.id*"|_"*crc1.id
    for o_owp in crc1.osss_owp
        push!(crc01.osss_owp,deepcopy(o_owp))
    end
    for o_mog in crc1.osss_mog
        push!(crc01.osss_mog,deepcopy(o_mog))
    end
    for mv_c in crc1.owp_MVcbls
        push!(crc01.owp_MVcbls,deepcopy(mv_c))
    end
    for hv_c in crc1.owp_HVcbls
        push!(crc01.owp_HVcbls,deepcopy(hv_c))
    end
    for o2o_c in crc1.oss2oss_cbls
        push!(crc01.oss2oss_cbls,deepcopy(o2o_c))
    end
    for pcc_c in crc1.pcc_cbls
        push!(crc01.pcc_cbls,deepcopy(pcc_c))
    end
    return crc01
end
