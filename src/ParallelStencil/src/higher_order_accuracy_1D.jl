## indices
import ..ParallelStencil: INDICES, WITHIN_DOC
const ix = INDICES[1];
const ixi = :($ix+1);
## 2nd-order accuracy for first-order derivatives
ix_1=:($ix+1);
macro dx_1(A::Symbol)  esc(:($A[$ix_1+1]-$A[$ix_1])) end
## 8th-order accuracy for first-order derivatives
ix_8=:($ix+3);

macro dx_8(A::Symbol)
    esc(:(1225/1024*($A[$ix_8+1]-$A[$ix_8])-
    245/3072*($A[$ix_8+2]-$A[$ix_8-1])+
    49/5120*($A[$ix_8+3]-$A[$ix_8-2])-
    5/7186*($A[$ix_8+4]-$A[$ix_8-3]) ));
end
## 12th-order accuracy for first-order derivatives
ix_12=:($ix+5);

macro dx_12(A::Symbol)
    esc(:(160083/131072*($A[$ix_12+1]-$A[$ix_12])-
    12705/131072*($A[$ix_12+2]-$A[$ix_12-1])+
    22869/1310720*($A[$ix_12+3]-$A[$ix_12-2])-
    5445/1835008*($A[$ix_12+4]-$A[$ix_12-3])+
    847/2359296*($A[$ix_12+5]-$A[$ix_12-4])-
    63/2883584*($A[$ix_12+6]-$A[$ix_12-5]) ));
end
