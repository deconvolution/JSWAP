## indices
import ..ParallelStencil: INDICES, WITHIN_DOC
ix, iy, iz = INDICES[1], INDICES[2], INDICES[3]
ixi, iyi, izi = :($ix+1), :($iy+1), :($iz+1);
## 2nd-order accuracy for first-order derivatives
ix_1=:($ix+1);
iy_1=:($iy+1);
iz_1=:($iz+1);

macro d_xi_1(A::Symbol)  esc(:( $A[$ix_1+1,$iyi ,$izi ] - $A[$ix_1  ,$iyi ,$izi ] )) end
macro d_yi_1(A::Symbol)  esc(:( $A[$ixi ,$iy_1+1,$izi ] - $A[$ixi ,$iy_1  ,$izi ] )) end
macro d_zi_1(A::Symbol)  esc(:( $A[$ixi ,$iyi ,$iz_1+1] - $A[$ixi ,$iyi ,$iz_1  ] )) end
## 2th-order central first-order derivatives
ix_c_4=:($ix+3);
iy_c_4=:($iy+3);
iz_c_4=:($iz+3);

macro d_xc_4(A::Symbol)
    esc(:( 1/12*$A[$ix_c_4-2,$iyi ,$izi]-2/3*$A[$ix_c_4-1,$iyi ,$izi]+
    2/3*$A[$ix_c_4+1,$iyi ,$izi]-1/12*$A[$ix_c_4+2,$iyi ,$izi]))
end

macro d_yc_4(A::Symbol)
    esc(:( 1/12*$A[$ixi,$iy_c_4-2,$izi]-3/2*$A[$ixi,$iy_c_4-1 ,$izi ]+
     2/3*$A[$ixi,$iy_c_4+1,$izi]-1/12*$A[$ixi,$iy_c_4+2 ,$izi ]))
end

macro d_zc_4(A::Symbol)
    esc(:( 1/12*$A[$ixi ,$iyi ,$iz_c_4-2]-3/2*$A[$ixi ,$iyi ,$iz_c_4-1]+
     3/2*$A[$ixi ,$iyi ,$iz_c_4+1]-1/12*$A[$ixi ,$iyi ,$iz_c_4+2] ))
end


## 8th-order accuracy for first-order derivatives
ix_8=:($ix+3);
iy_8=:($iy+3);
iz_8=:($iz+3);

macro d_xi_8(A::Symbol)
    esc(:(1225/1024*($A[$ix_8+1,$iyi,$izi]-$A[$ix_8,$iyi,$izi])-
    245/3072*($A[$ix_8+2,$iyi,$izi]-$A[$ix_8-1,$iyi,$izi])+
    49/5120*($A[$ix_8+3,$iyi,$izi]-$A[$ix_8-2,$iyi,$izi])-
    5/7186*($A[$ix_8+4,$iyi,$izi]-$A[$ix_8-3,$iyi,$izi])));
end

macro d_yi_8(A::Symbol)
    esc(:(1225/1024*($A[$ixi,$iy_8+1,$izi]-$A[$ixi,$iy_8,$izi])-
    245/3072*($A[$ixi,$iy_8+2,$izi]-$A[$ixi,$iy_8-1,$izi])+
    49/5120*($A[$ixi,$iy_8+3,$izi]-$A[$ixi,$iy_8-2,$izi])-
    5/7168*($A[$ixi,$iy_8+4,$izi]-$A[$ixi,$iy_8-3,$izi])));
end
macro d_zi_8(A::Symbol)
    esc(:(1225/1024*($A[$ixi,$iyi,$iz_8+1]-$A[$ixi,$iyi,$iz_8])-
    245/3072*($A[$ixi,$iyi,$iz_8+2]-$A[$ixi,$iyi,$iz_8-1])+
    49/5120*($A[$ixi,$iyi,$iz_8+3]-$A[$ixi,$iyi,$iz_8-2])-
    5/7168*($A[$ixi,$iyi,$iz_8+4]-$A[$ixi,$iyi,$iz_8-3])));
end
## 12th-order accuracy for first-order derivatives
ix_12=:($ix+5);
iy_12=:($iy+5);
iz_12=:($iz+5);

macro d_xi_12(A::Symbol)
    esc(:(160083/131072*($A[$ix_12+1,$iyi,$izi]-$A[$ix_12,$iyi,$izi])-
    12705/131072*($A[$ix_12+2,$iyi,$izi]-$A[$ix_12-1,$iyi,$izi])+
    22869/1310720*($A[$ix_12+3,$iyi,$izi]-$A[$ix_12-2,$iyi,$izi])-
    5445/1835008*($A[$ix_12+4,$iyi,$izi]-$A[$ix_12-3,$iyi,$izi])+
    847/2359296*($A[$ix_12+5,$iyi,$izi]-$A[$ix_12-4,$iyi,$izi])-
    63/2883584*($A[$ix_12+6,$iyi,$izi]-$A[$ix_12-5,$iyi,$izi])
    ));
end

macro d_yi_12(A::Symbol)
    esc(:(160083/131072*($A[$ixi,$iy_12+1,$izi]-$A[$ixi,$iy_12,$izi])-
    12705/131072*($A[$ixi,$iy_12+2,$izi]-$A[$ixi,$iy_12-1,$izi])+
    22869/1310720*($A[$ixi,$iy_12+3,$izi]-$A[$ixi,$iy_12-2,$izi])-
    5445/1835008*($A[$ixi,$iy_12+4,$izi]-$A[$ixi,$iy_12-3,$izi])+
    847/2359296*($A[$ixi,$iy_12+5,$izi]-$A[$ixi,$iy_12-4,$izi])-
    63/2883584*($A[$ixi,$iy_12+6,$izi]-$A[$ixi,$iy_12-5,$izi])
    ));
end
macro d_zi_12(A::Symbol)
    esc(:(160083/131072*($A[$ixi,$iyi,$iz_12+1]-$A[$ixi,$iyi,$iz_12])-
    12705/131072*($A[$ixi,$iyi,$iz_12+2]-$A[$ixi,$iyi,$iz_12-1])+
    22869/1310720*($A[$ixi,$iyi,$iz_12+3]-$A[$ixi,$iyi,$iz_12-2])-
    5445/1835008*($A[$ixi,$iyi,$iz_12+4]-$A[$ixi,$iyi,$iz_12-3])+
    847/2359296*($A[$ixi,$iyi,$iz_12+5]-$A[$ixi,$iyi,$iz_12-4])-
    63/2883584*($A[$ixi,$iyi,$iz_12+6]-$A[$ixi,$iyi,$iz_12-5])
    ));
end
