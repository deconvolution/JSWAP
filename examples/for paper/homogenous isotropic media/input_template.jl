## dimensions
# Time increment
dt=10.0^-3;
# dx
dx=15.0;
# dy
dy=15.0;
# dz
dz=15.0;
# number of time steps
nt=1000;
# nx
nx=80;
# ny
ny=80;
# nz
nz=90;
# 3D true coordinate X, Y and Z
Y,X,Z=JSWAP.meshgrid((1:80)*dx,(1:80)*dy,(1:90)*dz);
## material properties
lambda=ones(80,80,90)*10^9*1.0;
mu=ones(80,80,90)*10^9*1.0;
rho=ones(80,80,90)*1000.0;
inv_Qa=ones(80,80,90)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
r1=zeros(Int32,1,1);
r1[:] .=60;
# receiver grid location y
r2=zeros(Int32,1,1);
r2[:] .=50;
# receiver grid location z
r3=zeros(Int32,1,1);
r3[:] .=55;
# source grid location x
s1=zeros(Int32,1,1);
s1[:] .=40;
# source grid location y
s2=zeros(Int32,1,1);
s2[:] .=40;
# source grid location z
s3=zeros(Int32,1,1);
s3[:] .=40;
# moment tensor source
M11=zeros(nt,1);
M22=zeros(nt,1);
M33=zeros(nt,1);
M23=zeros(nt,1);
M13=zeros(nt,1);
M12=zeros(nt,1);
freq=10;
M11[:]=1*rickerWave(freq,dt,nt,2);
M22[:]=-1*rickerWave(freq,dt,nt,2);
M33[:]=0*rickerWave(freq,dt,nt,2);
M23[:]=0*rickerWave(freq,dt,nt,2);
M13[:]=0*rickerWave(freq,dt,nt,2);
M12[:]=0*rickerWave(freq,dt,nt,2);

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
lp=10;
# PML power
nPML=2;
# PML theorecital coefficient
Rc=.1;
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
