## import packages
using MAT,Plots,Dates,TimerOutputs,WriteVTK,DataFrames,CSV,
Statistics,DelimitedFiles,ProgressMeter,JSWAP
Threads.nthreads()
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
## read stiffness and density
data=JSWAP.readmat("./Campi/campi_model.mat","campi_model");
##
nx=convert(Int16,data["nx"]);
ny=convert(Int16,data["ny"]);
nz=convert(Int16,data["nz"]);

X=data["X"];
Y=data["Y"];
Z=data["Z"];

dx=data["dx"];
dy=data["dy"];
dz=data["dz"];

nt=2*10^4;
ns=nt;
dt=10.0^-3;

mutable struct C2
    lambda
    mu
    inv_Qa
    rho
end

rho=data["C"]["rho"];
lambda=data["C"]["C13"];
mu=data["C"]["C55"];

#lambda[:,:,1:30] .=0;
#mu[:,:,1:30] .=0;

inv_Qa=ones(nx,ny,nz)*0;

C=C2(lambda,
mu,
inv_Qa,
rho,
);
lambda=mu=rho=nothing;
## define model parameters

# PML layers
lp=15;

# xminus,xplus,yminus,yplus,zminus,zplus
pml_active=[1,1,1,1,0,1];

# PML coefficient, usually 2
nPML=2;

# Theoretical coefficient, more PML layers, less R
# Empirical values
# lp=[10,20,30,40]
# R=[.1,.01,.001,.0001]
Rc=.01;
## read source and receiver locations
r1t=zeros(Float32,1,length(data["r3"]));
r2t=copy(r1t);
r3t=copy(r1t);

r1=zeros(Int32,1,size(r1t,2));
r2=copy(r1);
r3=copy(r1);

r1[:]=data["r1"];
r2[:]=data["r2"];
r3[:]=data["r3"];

r1t=minimum(X) .+ dx .*r1;
r2t=minimum(Y) .+ dy .*r2;
r3t=maximum(Z) .- dz .*r3;
## source
# source location
# source location grid
s_s1=zeros(Int32,1,1);
s_s2=copy(s_s1);
s_s3=copy(s_s1);
s_s1[:] .=data["s1"];
s_s2[:] .=data["s2"];
s_s3[:] .=data["s3"];

# source locations true
s_s1t=minimum(X) .+ dx .*s_s1;
s_s2t=minimum(Y) .+ dy .*s_s2;
s_s3t=maximum(Z) .- dz .*s_s3;

# magnitude
M=2.7;
# source frequency [Hz]
freq=3;

# source signal
singles=JSWAP.rickerWave(freq,dt,nt,M);

# give source signal to x direction
s_src1=zeros(Float32,nt,1);

s_src2=copy(s_src1);
# give source signal to z direction
s_src3=copy(s_src1);

s_srcp=copy(s_src1);

s_src1[:] .=0;

s_src2[:] .=0;


s_src3[:] .=0;

s_srcp[:]=singles;

## plot
# point interval in time steps, 0 = no plot
plot_interval=200;
# save wavefield
wavefield_interval=0;
## mute some receiver components
Rm=ones(nt,length(r3),3);
## initialize seismograms
R1=zeros(Float32,nt,length(r3));
R2=copy(R1);
R3=copy(R1);
P=copy(R1);
data=nothing;
## implement solver
#for source_code=1:length(s_s3)
source_code=1;
global s1,s3,s1t,s3t,src1,src3,source_type,v1,v3,R1,R3,P,path,data,
path_pic,path_model,path_wavefield,path_rec
# source locations
s1=s_s1[source_code];
s2=s_s2[source_code];
s3=s_s3[source_code];

# path for this source
path=string(p3,"/source_code_",
(source_code));

path_pic=string(path,"/pic");
path_model=string(path,"/model");
path_wavefield=string(path,"/wavefield");
path_rec=string(path,"/rec");

s1=s_s1[source_code];
s2=s_s2[source_code];
s3=s_s3[source_code];

s1t=s_s1t[source_code];
s2t=s_s2t[source_code];
s3t=s_s3t[source_code];

src1=s_src1[:,source_code];
src2=s_src2[:,source_code];
src3=s_src3[:,source_code];
srcp=s_srcp[:,source_code];

# pass parameters to solver
v1,v2,v3,R1,R2,R3,P=JSWAP.CPU_3D.isotropic_forward_solver(dt,dx,dy,dz,nt,
nx,ny,nz,X,Y,Z,r1,r2,r3,s1,s2,s3,src1,src2,src3,srcp,
r1t,r2t,r3t,
Rm,
s1t,s2t,s3t,
lp,nPML,Rc,pml_active,
C,
plot_interval,
wavefield_interval,
path,
path_pic,
path_model,
path_wavefield,
path_rec);
