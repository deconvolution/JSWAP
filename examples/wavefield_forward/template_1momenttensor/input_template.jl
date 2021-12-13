## dimensions
# Time increment
input2.dt=10.0^-3;
# dx
input2.dx=15.0;
# dy
input2.dy=15.0;
# dz
input2.dz=15.0;
# number of time steps
input2.nt=1000;
# nx
input2.nx=80;
# ny
input2.ny=80;
# nz
input2.nz=90;
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z=JSWAP.meshgrid((1:80)*input2.dx,(1:80)*input2.dy,(1:90)*input2.dz);
## material properties
input2.lambda=ones(80,80,90)*10^9*1.0;
input2.mu=ones(80,80,90)*10^9*1.0;
input2.rho=ones(80,80,90)*1000.0;
input2.inv_Qa=ones(80,80,90)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
input2.r1=zeros(Int32,1,1);
input2.r1[:] .=50;
# receiver grid location y
input2.r2=zeros(Int32,1,1);
input2.r2[:] .=40;
# receiver grid location z
input2.r3=zeros(Int32,1,1);
input2.r3[:] .=45;
# source grid location x
input2.s1=zeros(Int32,1,1);
input2.s1[:] .=30;
# source grid location y
input2.s2=zeros(Int32,1,1);
input2.s2[:] .=30;
# source grid location z
input2.s3=zeros(Int32,1,1);
input2.s3[:] .=30;
# moment tensor source
input2.M11=zeros(input2.nt,1);
input2.M22=zeros(input2.nt,1);
input2.M33=zeros(input2.nt,1);
input2.M23=zeros(input2.nt,1);
input2.M13=zeros(input2.nt,1);
input2.M12=zeros(input2.nt,1);
freq=10;
input2.M11[:]=1*rickerWave(freq,input2.dt,input2.nt,2);
input2.M22[:]=-1*rickerWave(freq,input2.dt,input2.nt,2);
input2.M33[:]=0*rickerWave(freq,input2.dt,input2.nt,2);
input2.M23[:]=0*rickerWave(freq,input2.dt,input2.nt,2);
input2.M13[:]=0*rickerWave(freq,input2.dt,input2.nt,2);
input2.M12[:]=0*rickerWave(freq,input2.dt,input2.nt,2);

# receiver true location x
input2.r1t=input2.r1*input2.dx;
# receiver true location y
input2.r2t=input2.r2*input2.dy;
# receiver true location z
input2.r3t=input2.r3*input2.dz;
# source true location x
input2.s1t=input2.s1*input2.dx;
# source true location y
input2.s2t=input2.s2*input2.dx;
# source true location z
input2.s3t=input2.s3*input2.dx;
## PML
# PML layers
input2.lp=10;
# PML power
input2.nPML=2;
# PML theorecital coefficient
input2.Rc=.001;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
input2.PML_active=[1 1 1 1 1 1];
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
input2.wavefield_interval=nothing;