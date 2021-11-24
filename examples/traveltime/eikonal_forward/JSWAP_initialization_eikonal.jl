## import packages
using JSWAP
Threads.nthreads()
## input struct
mutable struct input3
    h
    nx
    ny
    nz
    X
    Y
    Z
    v
    r1
    r2
    r3
    s1
    s2
    s3
    r1t
    r2t
    r3t
    s1t
    s2t
    s3t
    path
end
input2=input3(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
