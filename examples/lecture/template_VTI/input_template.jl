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
nt=1000;
# nx
nx=100;
# ny
ny=100;
# nz
nz=100;
# 3D true coordinate X, Y and Z
Y,X,Z=JSWAP.meshgrid(1:10:1000,1:10:1000,1:10:1000);
## material properties
C11=ones(100,100,100)*10^9*1.506;
C12=ones(100,100,100)*10^9*(1.506-2*.4);
C13=ones(100,100,100)*10^9*.164;
C14=ones(100,100,100)*10^9*0;
C15=ones(100,100,100)*10^9*0;
C16=ones(100,100,100)*10^9*0;
C22=ones(100,100,100)*10^9*1.506;
C23=ones(100,100,100)*10^9*.164;
C24=ones(100,100,100)*10^9*0;
C25=ones(100,100,100)*10^9*0;
C26=ones(100,100,100)*10^9*0;
C33=ones(100,100,100)*10^9*1.084;
C34=ones(100,100,100)*10^9*0;
C35=ones(100,100,100)*10^9*0;
C36=ones(100,100,100)*10^9*0;
C44=ones(100,100,100)*10^9*.312;
C45=ones(100,100,100)*10^9*0;
C46=ones(100,100,100)*10^9*0;
C55=ones(100,100,100)*10^9*0.312;
C56=ones(100,100,100)*10^9*0;
C66=ones(100,100,100)*10^9*.4;

rho=ones(100,100,100)*1000.0;
inv_Qa=ones(100,100,100)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
r1=zeros(Int32,1,50);
r1[:]=30:79;
# receiver grid location y
r2=zeros(Int32,1,50);
r2[:]=30:79;
# receiver grid location z
r3=zeros(Int32,1,50);
r3 .=15;
# source grid location x
s1=zeros(Int32,1,1);
s1[:] .=51;
# source grid location y
s2=zeros(Int32,1,1);
s2[:] .=51;
# source grid location z
s3=zeros(Int32,1,1);
s3[:] .=51;
# source signal 1
src1=zeros(nt,1);
# source signal 2
src2=zeros(nt,1);
# source signal 3
src3=reshape(rickerWave(10,10^-3*1.0,1000,2),nt,1)*1;
# source signal pressure
srcp=reshape(rickerWave(10,10^-3*1.0,1000,2),nt,1)*0;
# receiver true location x
r1t=r1*dx;
# receiver true location y
r2t=r2*dy;
# receiver true location z
r3t=r3*dz;
# source true location x
s1t=s1*dx;
# source true location y
s2t=s2*dx;
# source true location z
s3t=s3*dx;
## PML
# PML layers
lp=15;
# PML power
nPML=2;
# PML theorecital coefficient
Rc=.01;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
PML_active=[1 1 1 1 1 1];
## plot
# path for storage. This must be the master branch of the following pathes.
path=p3;
# plot interval
plot_interval=100;
# wavefield interval
wavefield_interval=nothing;
