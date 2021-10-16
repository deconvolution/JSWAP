## read recordings
input2.R_true=JSWAP.readmat("./source_1_true/recording.mat","data");
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
input2.Y,input2.X,input2.Z=JSWAP.meshgrid((1:input2.ny)*input2.h,
(1:input2.nx)*input2.h,(1:input2.nz)*input2.h);
## material properties
input2.v=ones(input2.nx,input2.ny,input2.nz)*100;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
b,a=meshgrid(2:10:input2.ny,2:10:input2.nx);
a=reshape(a,length(a),1);
b=reshape(b,length(b),1);
# receiver grid location x
input2.r1=zeros(Int32,1,length(a));
input2.r1[:]=a;
# receiver grid location y
input2.r2=zeros(Int32,1,length(a));
input2.r2[:]=b;
# receiver grid location z
input2.r3=zeros(Int32,1,length(a));
input2.r3[:]=ones(1,length(a))*3;
# source grid location x
input2.s1=zeros(Int32,1,1);
input2.s1[:] .=50;
# source grid location y
input2.s2=zeros(Int32,1,1);
input2.s2[:] .=50;
# source grid location z
input2.s3=zeros(Int32,1,1);
input2.s3[:] .=70;

# receiver true location x
input2.r1t=minimum(input2.X) .+input2.r1*input2.h;
# receiver true location y
input2.r2t=minimum(input2.Y) .+input2.r2*input2.h;
# receiver true location z
input2.r3t=minimum(input2.Z) .+input2.r3*input2.h;
# source true location x
input2.s1t=minimum(input2.X) .+input2.s1*input2.h;
# source true location y
input2.s2t=minimum(input2.Y) .+input2.s2*input2.h;
# source true location z
input2.s3t=minimum(input2.Z) .+input2.s3*input2.h;
## plot
# path for storage. This must be the master branch of the following pathes.
input2.path="./source_1";

input2.write_rec=0;
