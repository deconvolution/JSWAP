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
input2.nx=100;
# ny
input2.ny=100;
# nz
input2.nz=100;
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z=JSWAP.meshgrid(1:10:1000,1:10:1000,1:10:1000);
## material properties
input2.C11=ones(100,100,100)*10^9*1.506;
input2.C12=ones(100,100,100)*10^9*(1.506-2*.4);
input2.C13=ones(100,100,100)*10^9*.164;
input2.C14=ones(100,100,100)*10^9*0;
input2.C15=ones(100,100,100)*10^9*0;
input2.C16=ones(100,100,100)*10^9*0;
input2.C22=ones(100,100,100)*10^9*1.506;
input2.C23=ones(100,100,100)*10^9*.164;
input2.C24=ones(100,100,100)*10^9*0;
input2.C25=ones(100,100,100)*10^9*0;
input2.C26=ones(100,100,100)*10^9*0;
input2.C33=ones(100,100,100)*10^9*1.084;
input2.C34=ones(100,100,100)*10^9*0;
input2.C35=ones(100,100,100)*10^9*0;
input2.C36=ones(100,100,100)*10^9*0;
input2.C44=ones(100,100,100)*10^9*.312;
input2.C45=ones(100,100,100)*10^9*0;
input2.C46=ones(100,100,100)*10^9*0;
input2.C55=ones(100,100,100)*10^9*0.312;
input2.C56=ones(100,100,100)*10^9*0;
input2.C66=ones(100,100,100)*10^9*.4;

input2.rho=ones(100,100,100)*1000.0;
input2.inv_Qa=ones(100,100,100)*0.0;
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
input2.s1[:] .=51;
# source grid location y
input2.s2=zeros(Int32,1,1);
input2.s2[:] .=51;
# source grid location z
input2.s3=zeros(Int32,1,1);
input2.s3[:] .=51;
# source signal x
input2.src1=zeros(input2.nt,1);
# source signal y
input2.src2=zeros(input2.nt,1);
# source signal z
input2.src3=reshape(rickerWave(10,10^-3*1.0,1000,2),input2.nt,1)*1;
# source signal pressure
input2.srcp=reshape(rickerWave(10,10^-3*1.0,1000,2),input2.nt,1)*0;
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
input2.lp=15;
# PML power
input2.nPML=2;
# PML theorecital coefficient
input2.Rc=.01;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
input2.PML_active=[1 1 1 1 1 1];
## plot
# path for storage. This must be the master branch of the following pathes.
input2.path="./P-source";
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
