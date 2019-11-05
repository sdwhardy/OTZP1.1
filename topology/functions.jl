#=
ocn=ocean
epty=top_mvTopos(ocn)
ocn.owpps[1].mv_zone.pos_height
=#

function top_mvTopos(owpps)
    #find all combinations
    clms=trunc(Int,length(owpps))
    rows=trunc(Int, 2.0^clms)
    empty_tbl=eensF_blankTbl(rows,clms)

    #remove anything less than a pair
    lngth=length(empty_tbl[:,1])
    indx=1
    while indx <= lngth
        if (sum(empty_tbl[indx,:]) < 1)
            empty_tbl=empty_tbl[1:size(empty_tbl,1) .!= indx,: ]
            indx=indx-1
        end
        indx=indx+1
        lngth=length(empty_tbl[:,1])
    end

    #remove those where the MV range does not overlap
    lngth=length(empty_tbl[:,1])
    indx=1
    while indx <= lngth
        closestOwpp=findfirst(isequal(1),empty_tbl[indx,:])
        farthestOwpp=findlast(isequal(1),empty_tbl[indx,:])
        ymx0=owpps[closestOwpp].node.xy.y+owpps[closestOwpp].mv_zone.pos_height
        ymn0=owpps[closestOwpp].node.xy.y-owpps[closestOwpp].mv_zone.neg_height
        xmx0=owpps[closestOwpp].node.xy.x+owpps[closestOwpp].mv_zone.pos_width
        xmn0=owpps[closestOwpp].node.xy.x-owpps[closestOwpp].mv_zone.neg_width
        ymx1=owpps[farthestOwpp].node.xy.y+owpps[farthestOwpp].mv_zone.pos_height
        ymn1=owpps[farthestOwpp].node.xy.y-owpps[farthestOwpp].mv_zone.neg_height
        xmx1=owpps[farthestOwpp].node.xy.x+owpps[farthestOwpp].mv_zone.pos_width
        xmn1=owpps[farthestOwpp].node.xy.x-owpps[farthestOwpp].mv_zone.neg_width
        if ((ymx1>ymn0 && ymx0>ymn1 && xmx1>xmn0 && xmx0>xmn1) != true)
            empty_tbl=empty_tbl[1:size(empty_tbl,1) .!= indx,: ]
            indx=indx-1
        end
        indx=indx+1
        lngth=length(empty_tbl[:,1])
    end

    return empty_tbl
end
############################################################################
#=
bn=[1,1,1,0,1,1]
dec=cir_bin2dec(bn)
bn=cir_dec2bin(dec)
=#
function top_bin2dec(bn)
    dec=0.0
    for (i,bt) in enumerate(bn)
        dec=dec+bt*(2^(i-1))
    end
    return dec
end
############################################################################
function top_dec2bin(dec)
    bn=Int8[]
    while (dec != 0)
        push!(bn,mod(dec,2))
        dec=floor(Int64, dec/2)
    end
    return bn
end
############################################################################
