## dimensions
nx=80;
ny=80;
nz=80;
h=100;
# grid true location
Y,X,Z=meshgrid(h:h:h*ny,h:h:h*nx,h:h:h*nz);
# velocity
v=zeros(nx,ny,nz);
v[:] .=5000;

# source grid location
s1=zeros(Int64,1,1);
s2=zeros(Int64,1,1);
s3=zeros(Int64,1,1);

s1[:] .=20;
s2[:] .=60;
s3[:] .=70;

# source true location
s1t=h*s1;
s2t=h*s2;
s3t=h*s3;

# receiver grid location
r1=zeros(Int64,1,5);
r2=zeros(Int64,1,5);
r3=zeros(Int64,1,5);

r1[:]=[10 20 30 40 50];
r2[:]=[20 30 40 50 60];
r3[:]=[3 3 3 3 3];

# receiver true location
r1t=h*r1;
r2t=h*r2;
r3t=h*r3;

# path to save output
path=p3;

# write traveltime and recording?
write_t=1;
