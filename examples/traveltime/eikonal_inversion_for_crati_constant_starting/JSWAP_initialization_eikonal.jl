## import packages
using JSWAP,MATLAB,FileIO
Threads.nthreads()
## input struct
mutable struct input3
    T
    R_true
    R_cal
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
    write_rec
    n_iteration
    max_gradient
    fu
end
input2=input3(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
