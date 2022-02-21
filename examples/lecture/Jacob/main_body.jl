## import packages
using JSWAP,DelimitedFiles,Statistics
ti=JSWAP.TimerOutput();
Threads.nthreads()
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
## Run solvers
include("./input_template.jl");
@JSWAP.timeit ti "simulation" JSWAP.CPU_3D.JSWAP_CPU_3D_forward_isotropic_solver(nt=nt,
nx=nx,
ny=ny,
nz=nz,
dt=dt,
dx=dx,
dy=dy,
dz=dz,
X=X,
Y=Y,
Z=Z,
lambda=lambda,
mu=mu,
rho=rho,
inv_Qa=inv_Qa,
s1=s1,
s2=s2,
s3=s3,
s1t=s1t,
s2t=s2t,
s3t=s3t,
r1=r1,
r2=r2,
r3=r3,
r1t=r1t,
r2t=r2t,
r3t=r3t,
path=path,
wavefield_interval=wavefield_interval,
plot_interval=plot_interval,
src1=src1,
src2=src2,
src3=src3,
srcp=srcp,
lp=lp,
Rc=Rc,
nPML=nPML,
PML_active=PML_active);
