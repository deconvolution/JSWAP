## import packages
using JSWAP
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
JSWAP.CPU_3D.JSWAP_CPU_3D_forward_isotropic_solver(nt=nt,
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
M11=M11,
M22=M22,
M33=M33,
M23=M23,
M13=M13,
M12=M12,
lp=lp,
Rc=Rc,
nPML=nPML,
PML_active=PML_active);
## load synthetic recordings
synthetic_rec_location="./input_campi_3Hz_sensitivity_constant_model/rec";
R1=JSWAP.readmat(string(synthetic_rec_location,"/rec_1.mat"),"data");
R2=JSWAP.readmat(string(synthetic_rec_location,"/rec_2.mat"),"data");
R3=JSWAP.readmat(string(synthetic_rec_location,"/rec_3.mat"),"data");
RP=JSWAP.readmat(string(synthetic_rec_location,"/rec_p.mat"),"data");
## load true recordings
true_rec_location="../campi 3Hz/main_body/rec";
tR1=JSWAP.readmat(string(true_rec_location,"/rec_1.mat"),"data");
tR2=JSWAP.readmat(string(true_rec_location,"/rec_2.mat"),"data");
tR3=JSWAP.readmat(string(true_rec_location,"/rec_3.mat"),"data");
tP=JSWAP.readmat(string(true_rec_location,"/rec_p.mat"),"data");
## adjoint simulation
i=5;
misR1=R1-tR1;
misR2=R2-tR2;
misR3=R3-tR3;
misP=0*(RP-tP);

rev_misR1=reverse(misR1,dims=1);
rev_misR2=reverse(misR2,dims=1);
rev_misR3=reverse(misR3,dims=1);
rev_misP=reverse(misP,dims=1);

JSWAP.CPU_3D.JSWAP_CPU_3D_adjoint_isotropic_solver(dt=dt,
dx=dx,
dy=dy,
dz=dz,
nt=nt,
nx=nx,
ny=ny,
nz=nz,
X=X,
Y=Y,
Z=Z,
r1=r1,
r2=r2,
r3=r3,
s1=r1[i],
s2=r2[i],
s3=r3[i],
src1=rev_misR1[:,i],
src2=rev_misR2[:,i],
src3=rev_misR3[:,i],
srcp=rev_misP[:,i],
r1t=r1t[i],
r2t=r2t[i],
r3t=r3t[i],
s1t=r1t,
s2t=r2t,
s3t=r3t,
lp=lp,
nPML=nPML,
Rc=Rc,
PML_active=PML_active,
lambda=lambda,
mu=mu,
rho=rho,
plot_interval=0,
wavefield_interval=wavefield_interval,
path="./adjoint_wavefield");
## wavefield correlation
path_forward_wavefield=string("./input_campi_3Hz_sensitivity_constant_model/wavefield");
path_adjoint_wavefield=string("./adjoint_wavefield/wavefield_adjoint");
JSWAP.CPU_3D.JSWAP_CPU_3D_compute_sensitivity(path="./sensitivity_kernel/",
nt=nt,
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
path_forward_wavefield=path_forward_wavefield,
path_adjoint_wavefield=path_adjoint_wavefield,
wavefield_interval=50);
