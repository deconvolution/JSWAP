## dimensions
# Time increment
input2.dt=10.0^-3;
# dx
input2.dx=10.0;
# dy
input2.dy=10.0;
# dz
input2.dz=10.0;
# number of time steps
input2.nt=1000;
# nx
input2.nx=80;
# ny
input2.ny=80;
# nz
input2.nz=90;
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z=JSWAP.meshgrid(1:80,1:80,1:90);
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
input2.r1=zeros(Int32,1,50);
input2.r1[:]=30:79;
# receiver grid location y
input2.r2=zeros(Int32,1,50);
input2.r2[:]=30:79;
# receiver grid location z
input2.r3=zeros(Int32,1,50);
input2.r3 .=15;
# source grid location x
input2.s1=zeros(Int32,1,1);
input2.s1[:] .=60;
# source grid location y
input2.s2=zeros(Int32,1,1);
input2.s2[:] .=50;
# source grid location z
input2.s3=zeros(Int32,1,1);
input2.s3[:] .=30;
# source signal x
input2.src1=zeros(input2.nt,1);
# source signal y
input2.src2=zeros(input2.nt,1);
# source signal z
input2.src3=reshape(rickerWave(2,10^-3*.1,1000,2),input2.nt,1);
# source signal pressure
input2.srcp=zeros(input2.nt,1);
# receiver true location x
input2.r1t=input2.r1*input2.dx;
# receiver true location y
input2.r2t=input2.r2*input2.dy;
# receiver true location z
input2.r3t=input2.r3*input2.dz;
# activate receivers
input2.Rm=ones(input2.nt,length(input2.r3),4);
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
input2.Rc=.01;
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
input2.wavefield_interval=0;
