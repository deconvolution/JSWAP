## dimensions
# Time increment
dt=10.0^-3;
# dx
dx=10.0;
# dy
dy=10.0;
# dz
dz=10.0;
# number of time steps
nt=30;
# nx
nx=50;
# ny
ny=50;
# nz
nz=50;
# 3D true coordinate X, Y and Z
Y,X,Z=JSWAP.meshgrid(1:500,1:500,1:500);
## material properties
lambda=ones(50,50,50)*10^9*1.0;
mu=ones(50,50,50)*10^9*1.0;
rho=ones(50,50,50)*1000.0;
inv_Qa=ones(50,50,50)*10.0^-4;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
r1=zeros(Int64,1,2);
r1[:]=[20 30];
# receiver grid location y
r2=zeros(Int64,1,2);
r2[:]=[20 30];
# receiver grid location z
r3=zeros(Int64,1,2);
r3 .=15;
# source grid location x
s1=zeros(Int64,1,1);
s1[:] .=30;
# source grid location y
s2=zeros(Int64,1,1);
s2[:] .=30;
# source grid location z
s3=zeros(Int64,1,1);
s3[:] .=30;
# source signal x
src1=zeros(nt,1);
# source signal y
src2=zeros(nt,1);
# source signal z
src3=reshape(rickerWave(20.0,10.0^-3*1.0,30,2),nt,1);
# source signal pressure
srcp=zeros(nt,1);
# receiver true location x
r1t=r1*dx;
# receiver true location y
r2t=r2*dy;
# receiver true location z
r3t=r3*dz;
# activate receivers
Rm=ones(length(r3),4);
# source true location x
s1t=s1*dx;
# source true location y
s2t=s2*dx;
# source true location z
s3t=s3*dx;
## PML
# PML layers
lp=10;
# PML power
nPML=2;
# PML theorecital coefficient
Rc=.01;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
PML_active=[1 1 1 1 1 1];
## plot
path=nothing;
# plot interval
plot_interval=nothing;
# wavefield interval
wavefield_interval=nothing;
