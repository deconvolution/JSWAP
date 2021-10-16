## read stiffness and density
## dimensions
# dx
input2.h=10;
# nx
input2.nx=100;
# ny
input2.ny=90;
# nz
input2.nz=80;
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z=JSWAP.meshgrid(1:input2.ny,1:input2.nx,1:input2.nz);
## material properties
input2.v=ones(input2.nx,input2.ny,input2.nz)*100;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
input2.r1=zeros(Int32,1,3);
input2.r1[:]=[3,5,10];
# receiver grid location y
input2.r2=zeros(Int32,1,3);
input2.r2[:]=[3,5,8];
# receiver grid location z
input2.r3=zeros(Int32,1,3);
input2.r3[:]=[4,5,6];
# source grid location x
input2.s1=zeros(Int32,1,1);
input2.s1[:] .=50;
# source grid location y
input2.s2=zeros(Int32,1,1);
input2.s2[:] .=50;
# source grid location z
input2.s3=zeros(Int32,1,1);
input2.s3[:] .=30;

# receiver true location x
input2.r1t=minimum(input2.X) .+input2.r1*input2.h;
# receiver true location y
input2.r2t=minimum(input2.Y) .+input2.r2*input2.h;
# receiver true location z
input2.r3t=maximum(input2.Z) .-input2.r3*input2.h;
# source true location x
input2.s1t=minimum(input2.X) .+input2.s1*input2.h;
# source true location y
input2.s2t=minimum(input2.Y) .+input2.s2*input2.h;
# source true location z
input2.s3t=minimum(input2.Z) .-input2.s3*input2.h;
## plot
# path for storage. This must be the master branch of the following pathes.
input2.path="./source_1";
