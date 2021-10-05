## read stiffness and density
data=JSWAP.readmat("./m.mat","data");
## dimensions
# Time increment
input2.dt=10.0^-3;
# dx
input2.dx=data["dx"];
# dy
input2.dy=data["dy"];
# dz
input2.dz=data["dz"];
# number of time steps
input2.nt=20000;
# nx
input2.nx=convert(Int32,data["nx"]);;
# ny
input2.ny=convert(Int32,data["ny"]);;
# nz
input2.nz=convert(Int32,data["nz"]);;
# 3D true coordinate X, Y and Z
input2.X=data["X"];
input2.Y=data["Y"];
input2.Z=data["Z"];
## material properties
vp=data["vp"];
vs=data["vs"];
rho=data["rho"];
input2.lambda=rho .*vp.^2 .- 2 .*rho .*vs.^2;
input2.mu=rho .*vs.^2;
input2.rho=data["rho"];
input2.inv_Qa=ones(input2.nx,input2.ny,input2.nz)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
input2.r1=zeros(Int32,1,size(data["r"],1));
input2.r1[:]=data["r"][:,1];
# receiver grid location y
input2.r2=zeros(Int32,size(input2.r1));
input2.r2[:]=data["r"][:,2];
# receiver grid location z
input2.r3=zeros(Int32,size(input2.r1));
input2.r3[:]=data["r"][:,3];
# source grid location x
input2.s1=zeros(Int32,1,size(data["s"],1));
input2.s1[:] .=data["s"][1];
# source grid location y
input2.s2=zeros(Int32,1,size(data["s"],1));
input2.s2[:] .=data["s"][2];
# source grid location z
input2.s3=zeros(Int32,1,size(data["s"],1));
input2.s3[:] .=data["s"][3];

# moment tensor source
freq=.5;
input2.M11=zeros(input2.nt,1);
input2.M22=zeros(input2.nt,1);
input2.M33=zeros(input2.nt,1);
input2.M23=zeros(input2.nt,1);
input2.M13=zeros(input2.nt,1);
input2.M12=zeros(input2.nt,1);
input2.M11[:]=rickerWave(freq,input2.dt,input2.nt,1)*data["moment_tensor"][1,1];
input2.M22[:]=rickerWave(freq,input2.dt,input2.nt,1)*data["moment_tensor"][2,2];
input2.M33[:]=rickerWave(freq,input2.dt,input2.nt,1)*data["moment_tensor"][3,3];
input2.M23[:]=rickerWave(freq,input2.dt,input2.nt,1)*data["moment_tensor"][2,3];
input2.M13[:]=rickerWave(freq,input2.dt,input2.nt,1)*data["moment_tensor"][1,3];
input2.M12[:]=rickerWave(freq,input2.dt,input2.nt,1)*data["moment_tensor"][1,2];


# receiver true location x
input2.r1t=zeros(Float64,1,size(data["r"],1));
input2.r1t[:]=data["rt"][:,1];
# receiver true location y
input2.r2t=zeros(Float64,size(input2.r1));
input2.r2t[:]=data["rt"][:,2];
# receiver true location z
input2.r3t=zeros(Float64,size(input2.r1));
input2.r3t[:]=data["rt"][:,3];
# source true location x
input2.s1t=zeros(Float64,1,size(data["s"],1));
input2.s1t[:] .=data["st"][1];
# source true location y
input2.s2t=zeros(Float64,1,size(data["s"],1));
input2.s2t[:] .=data["st"][2];
# source true location z
input2.s3t=zeros(Float64,1,size(data["s"],1));
input2.s3t[:] .=data["st"][3];
## mute?
input2.Rm=ones(length(input2.r3),4);
## PML
# PML layers
input2.lp=10;
# PML power
input2.nPML=2;
# PML theorecital coefficient
input2.Rc=.1;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
input2.PML_active=[1 1 1 1 0 1];
## plot
# path for storage. This must be the master branch of the following pathes.
input2.path="./normal";
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
