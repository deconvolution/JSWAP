## Initialization of JSWAP
include("./JSWAP_initialization.jl");
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
## Specify where the input file is
path_to_input=["./input_template.jl"];
## Run solvers
for I=1:length(path_to_input)
    include(path_to_input[I]);
    JSWAP.CPU_3D.JSWAP_CPU_3D_forward_isotropic_solver(nt=input2.nt,
    nx=input2.nx,
    ny=input2.ny,
    nz=input2.nz,
    dt=input2.dt,
    dx=input2.dx,
    dy=input2.dy,
    dz=input2.dz,
    X=input2.X,
    Y=input2.Y,
    Z=input2.Z,
    lambda=input2.lambda,
    mu=input2.mu,
    rho=input2.rho,
    inv_Qa=input2.inv_Qa,
    s1=input2.s1,
    s2=input2.s2,
    s3=input2.s3,
    s1t=input2.s1t,
    s2t=input2.s2t,
    s3t=input2.s3t,
    r1=input2.r1,
    r2=input2.r2,
    r3=input2.r3,
    r1t=input2.r1t,
    r2t=input2.r2t,
    r3t=input2.r3t,
    path=input2.path,
    wavefield_interval=input2.wavefield_interval,
    plot_interval=input2.plot_interval,
    M11=input2.M11,
    M22=input2.M22,
    M33=input2.M33,
    M23=input2.M23,
    M13=input2.M13,
    M12=input2.M12,
    lp=input2.lp,
    Rc=input2.Rc,
    nPML=input2.nPML,
    PML_active=input2.PML_active);
end
