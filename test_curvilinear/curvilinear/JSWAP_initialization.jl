## import packages
## Using ParallelStencil
include("../ParallelStencil/ParallelStencil.jl");
## utilities
include("..//utilities/utilities.jl");
## using packages
using Random,MAT,Plots,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,CSV,
.ParallelStencil,.ParallelStencil.FiniteDifferences3D
## input struct
mutable struct input3
    dt
    dx
    dy
    dz
    nt
    nx
    ny
    nz
    X
    Y
    Z
    lambda
    mu
    rho
    inv_Qa
    r1
    r2
    r3
    s1
    s2
    s3
    src1
    src2
    src3
    srcp
    r1t
    r2t
    r3t
    Rm
    s1t
    s2t
    s3t
    lp
    nPML
    Rc
    PML_active
    plot_interval
    wavefield_interval
    path
    path_pic
    path_model
    path_wavefield
    path_rec
    Kmax
end
input2=input3(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
