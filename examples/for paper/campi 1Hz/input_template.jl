## read stiffness and density
data=JSWAP.readmat("./campi_model_random.mat","campi_model");
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
nt=20000;
# nx
nx=convert(Int32,data["nx"]);;
# ny
ny=convert(Int32,data["ny"]);;
# nz
nz=convert(Int32,data["nz"]);;
# 3D true coordinate X, Y and Z
X=data["X"];
Y=data["Y"];
Z=data["Z"];
## material properties
lambda=data["C"]["C13"];
mu=data["C"]["C55"];
rho=data["C"]["rho"];
inv_Qa=ones(nx,ny,nz)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
r1=zeros(Int32,1,length(data["r3"]));
r1[:]=data["r1"];
# receiver grid location y
r2=zeros(Int32,size(r1));
r2[:]=data["r2"];
# receiver grid location z
r3=zeros(Int32,size(r1));
r3[:]=data["r3"];
# source grid location x
s1=zeros(Int32,1,length(data["s3"]));
s1[:] .=data["s1"];
# source grid location y
s2=zeros(Int32,1,length(data["s3"]));
s2[:] .=data["s2"];
# source grid location z
s3=zeros(Int32,1,length(data["s3"]));
s3[:] .=data["s3"];;
# moment tensor source
M11=zeros(nt,1);
M22=zeros(nt,1);
M33=zeros(nt,1);
M23=zeros(nt,1);
M13=zeros(nt,1);
M12=zeros(nt,1);
M11[:]=rickerWave(1,dt,nt,2);
M22[:]=-rickerWave(1,dt,nt,2);
M33[:]=rickerWave(1,dt,nt,2)*0.0;
M23[:]=rickerWave(1,dt,nt,2)*0.0;
M13[:]=rickerWave(1,dt,nt,2)*0.0;
M12[:]=rickerWave(1,dt,nt,2)*0.0;

# receiver true location x
r1t=minimum(X) .+r1*dx;
# receiver true location y
r2t=minimum(Y) .+r2*dy;
# receiver true location z
r3t=maximum(Z) .-r3*dz;
# source true location x
s1t=minimum(X) .+s1*dx;
# source true location y
s2t=minimum(Y) .+s2*dy;
# source true location z
s3t=minimum(Z) .-s3*dz;
## PML
# PML layers
lp=10;
# PML power
nPML=2;
# PML theorecital coefficient
Rc=.1;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
PML_active=[1 1 1 1 0 1];
## plot
# path for storage. This must be the master branch of the following pathes.
path=p3;
# plot interval
plot_interval=100;
# wavefield interval
wavefield_interval=nothing;
