## indices
import ..ParallelStencil: INDICES, WITHIN_DOC
ix, iy, iz = INDICES[1], INDICES[2], INDICES[3]
## 2nd-order accuracy for first-order derivatives
ix_1=:($ix+1);
iy_1=:($iy+1);
iz_1=:($iz+1);

macro   dx_1(A::Symbol)  esc(:( $A[$ix_1+1,$iy ,$iz ] - $A[$ix_1  ,$iy ,$iz ] )) end
macro   dy_1(A::Symbol)  esc(:( $A[$ix ,$iy_1+1,$iz ] - $A[$ix ,$iy_1  ,$iz ] )) end
macro   dz_1(A::Symbol)  esc(:( $A[$ix ,$iy ,$iz_1+1] - $A[$ix ,$iy ,$iz_1  ] )) end
## 8th-order accuracy for first-order derivatives
ix_8=:($ix+3);
iy_8=:($iy+3);
iz_8=:($iz+3);

macro   dx_8(A::Symbol)
    esc(:(1225/1024*($A[$ix_8+1,$iy,$iz]-$A[$ix_8,$iy,$iz])-
    245/3072*($A[$ix_8+2,$iy,$iz]-$A[$ix_8-1,$iy,$iz])+
    49/5120*($A[$ix_8+3,$iy,$iz]-$A[$ix_8-2,$iy,$iz])-
    5/7186*($A[$ix_8+4,$iy,$iz]-$A[$ix_8-3,$iy,$iz])));
end

macro   dy_8(A::Symbol)
    esc(:(1225/1024*($A[$ix,$iy_8+1,$iz]-$A[$ix,$iy_8,$iz])-
    245/3072*($A[$ix,$iy_8+2,$iz]-$A[$ix,$iy_8-1,$iz])+
    49/5120*($A[$ix,$iy_8+3,$iz]-$A[$ix,$iy_8-2,$iz])-
    5/7168*($A[$ix,$iy_8+4,$iz]-$A[$ix,$iy_8-3,$iz])));
end
macro   dz_8(A::Symbol)
    esc(:(1225/1024*($A[$ix,$iy,$iz_8+1]-$A[$ix,$iy,$iz_8])-
    245/3072*($A[$ix,$iy,$iz_8+2]-$A[$ix,$iy,$iz_8-1])+
    49/5120*($A[$ix,$iy,$iz_8+3]-$A[$ix,$iy,$iz_8-2])-
    5/7168*($A[$ix,$iy,$iz_8+4]-$A[$ix,$iy,$iz_8-3])));
end
## 12th-order accuracy for first-order derivatives
ix_12=:($ix+5);
iy_12=:($iy+5);
iz_12=:($iz+5);

macro   dx_12(A::Symbol)
    esc(:(160083/131072*($A[$ix_12+1,$iy,$iz]-$A[$ix_12,$iy,$iz])-
    12705/131072*($A[$ix_12+2,$iy,$iz]-$A[$ix_12-1,$iy,$iz])+
    22869/1310720*($A[$ix_12+3,$iy,$iz]-$A[$ix_12-2,$iy,$iz])-
    5445/1835008*($A[$ix_12+4,$iy,$iz]-$A[$ix_12-3,$iy,$iz])+
    847/2359296*($A[$ix_12+5,$iy,$iz]-$A[$ix_12-4,$iy,$iz])-
    63/2883584*($A[$ix_12+6,$iy,$iz]-$A[$ix_12-5,$iy,$iz])
    ));
end

macro   dy_12(A::Symbol)
    esc(:(160083/131072*($A[$ix,$iy_12+1,$iz]-$A[$ix,$iy_12,$iz])-
    12705/131072*($A[$ix,$iy_12+2,$iz]-$A[$ix,$iy_12-1,$iz])+
    22869/1310720*($A[$ix,$iy_12+3,$iz]-$A[$ix,$iy_12-2,$iz])-
    5445/1835008*($A[$ix,$iy_12+4,$iz]-$A[$ix,$iy_12-3,$iz])+
    847/2359296*($A[$ix,$iy_12+5,$iz]-$A[$ix,$iy_12-4,$iz])-
    63/2883584*($A[$ix,$iy_12+6,$iz]-$A[$ix,$iy_12-5,$iz])
    ));
end
macro   dz_12(A::Symbol)
    esc(:(160083/131072*($A[$ix,$iy,$iz_12+1]-$A[$ix,$iy,$iz_12])-
    12705/131072*($A[$ix,$iy,$iz_12+2]-$A[$ix,$iy,$iz_12-1])+
    22869/1310720*($A[$ix,$iy,$iz_12+3]-$A[$ix,$iy,$iz_12-2])-
    5445/1835008*($A[$ix,$iy,$iz_12+4]-$A[$ix,$iy,$iz_12-3])+
    847/2359296*($A[$ix,$iy,$iz_12+5]-$A[$ix,$iy,$iz_12-4])-
    63/2883584*($A[$ix,$iy,$iz_12+6]-$A[$ix,$iy,$iz_12-5])
    ));
end
## for curvilinear
ix_1=:($ix+1);
iy_1=:($iy+1);
iz_1=:($iz+1);

macro   cur1(A::Symbol)
    esc(:( ($A[$ix_1,$iy ,$iz_1+1 ] + $A[$ix_1-1  ,$iy ,$iz_1-1 ] -
    $A[$ix_1,$iy ,$iz_1+1 ] - $A[$ix_1-1  ,$iy ,$iz_1-1 ])/4 ))
end

macro   cur3(A::Symbol)
    esc(:( ($A[$ix,$iy_1 ,$iz_1+1 ] + $A[$ix  ,$iy_1-1 ,$iz_1-1 ] -
    $A[$ix,$iy_1 ,$iz_1+1 ] - $A[$ix  ,$iy_1-1 ,$iz_1-1 ])/4 ))
end
