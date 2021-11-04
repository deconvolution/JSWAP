## dimensions
data=JSWAP.readmat("./example.mat","data");
# Time increment
input2.dt=10.0^-3;
# dx
input2.dx=data["dx"];
# dy
input2.dy=data["dy"];
# dz
input2.dz=data["dz"];
# number of time steps
input2.nt=800;
# nx
input2.nx=round(Int32,data["nx"]);
# ny
input2.ny=round(Int32,data["ny"]);
# nz
input2.nz=round(Int32,data["nz"]);

# topography
input2.Kmax=data["Kmax"];
# input2.Kmax[:] .=900;
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z=JSWAP.meshgrid((1:input2.nx)*input2.dx,(1:input2.ny)*input2.dy,(1:input2.nz)*input2.dz);
## material properties
input2.lambda=data["lambda"];
input2.mu=data["mu"];
input2.rho=data["rho"];
input2.inv_Qa=ones(input2.nx,input2.ny,input2.nz)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
input2.r1=zeros(Int32,1,1);
input2.r1[:] .=data["r1"];;
# receiver grid location y
input2.r2=zeros(Int32,1,1);
input2.r2[:] .=data["r2"];
# receiver grid location z
input2.r3=zeros(Int32,1,1);
input2.r3[:] .=data["r3"];
# source grid location x
input2.s1=zeros(Int32,1,1);
input2.s1[:] .=data["s1"];
# source grid location y
input2.s2=zeros(Int32,1,1);
input2.s2[:] .=data["s2"];
# source grid location z
input2.s3=zeros(Int32,1,1);
input2.s3[:] .=data["s3"];
# point source
freq=15;
input2.src1=zeros(input2.nt,1);
input2.src2=zeros(input2.nt,1);
input2.src3=zeros(input2.nt,1);
input2.srcp=zeros(input2.nt,1);
input2.src1[:]=0*rickerWave(freq,input2.dt,input2.nt,2);
input2.src2[:]=0*rickerWave(freq,input2.dt,input2.nt,2);
input2.src3[:]=1*rickerWave(freq,input2.dt,input2.nt,2);
input2.srcp[:]=0*rickerWave(freq,input2.dt,input2.nt,2);


# receiver true location x
input2.r1t=input2.r1*input2.dx;
# receiver true location y
input2.r2t=input2.r2*input2.dy;
# receiver true location z
input2.r3t=input2.r3*input2.dz;
# activate receivers
input2.Rm=ones(length(input2.r3),4);
# source true location x
input2.s1t=input2.s1*input2.dx;
# source true location y
input2.s2t=input2.s2*input2.dx;
# source true location z
input2.s3t=input2.s3*input2.dx;
## PML
# PML layers
input2.lp=20;
# PML power
input2.nPML=2;
# PML theorecital coefficient
input2.Rc=.001;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
input2.PML_active=[1 1 1 1 1 0];
## plot
# path for storage. This must be the master branch of the following pathes.
input2.path=p3;
# path for wavefield .vtk
input2.path_pic=string(input2.path,"/pic");
# path for model
input2.path_model=string(input2.path,"/model");
# path for wavefield .mat
input2.path_wavefield=string(input2.path,"/wavefield");
# path for recordings
input2.path_rec=string(input2.path,"/rec");
# plot interval
input2.plot_interval=100;
# wavefield interval
input2.wavefield_interval=0;
