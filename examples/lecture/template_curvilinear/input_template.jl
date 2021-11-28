## read material properties
data=JSWAP.readmat("./example.mat","data");
## dimensions
# Time increment
dt=10.0^-3;
# dx
dx=data["dx"];
# dy
dy=data["dy"];
# dz
dz=data["dz"];
# number of time steps
nt=1000;
# nx
nx=floor(Int64,data["nx"]);
# ny
ny=floor(Int64,data["ny"]);
# nz
nz=floor(Int64,data["nz"]);
# 3D true coordinate X, Y and Z
X=data["X"];
Y=data["Y"];
K=data["K"];
Kmax=data["Kmax"];
## material properties
lambda=data["lambda"];
mu=data["mu"];
rho=data["rho"];
inv_Qa=ones(nx,ny,nz)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
r1=zeros(Int64,1,1);
r1[:] .=data["r1"];
# receiver grid location y
r2=zeros(Int64,1,1);
r2[:] .=data["r2"];
# receiver grid location z
r3=zeros(Int64,1,1);
r3[:] .=data["r3"];
# source grid location x
s1=zeros(Int64,1,1);
s1[:] .=data["s1"];
# source grid location y
s2=zeros(Int64,1,1);
s2[:] .=data["s2"];
# source grid location z
s3=zeros(Int64,1,1);
s3[:] .=data["s3"];
# Directional source
src1=zeros(nt,1);
src2=zeros(nt,1);
src3=zeros(nt,1);
srcp=zeros(nt,1);
freq=15;
src1[:]=0*rickerWave(freq,dt,nt,2);
src2[:]=0*rickerWave(freq,dt,nt,2);
src3[:]=1*rickerWave(freq,dt,nt,2);
srcp[:]=0*rickerWave(freq,dt,nt,2);

r1t=zeros(Int64,1,1);
r2t=zeros(Int64,1,1);
r3t=zeros(Int64,1,1);
s1t=zeros(Int64,1,1);
s2t=zeros(Int64,1,1);
s3t=zeros(Int64,1,1);
# receiver true location x
r1t[:] .=data["r1t"];
# receiver true location y
r2t[:] .=data["r2t"];
# receiver true location z
r3t[:] .=data["r3t"];
# source true location x
s1t[:] .=data["s1t"];
# source true location y
s2t[:] .=data["s2t"];
# source true location z
s3t[:] .=data["s3t"];
## PML
# PML layers
lp=10;
# PML power
nPML=2;
# PML theorecital coefficient
Rc=.001;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
PML_active=[1 1 1 1 1 0];
## plot
# path for storage. This must be the master branch of the following pathes.
path=p3;
# plot interval
plot_interval=100;
# wavefield interval
wavefield_interval=nothing;
