## read data
nx=80;
ny=80;
nz=80;

h=100;
Y,X,Z=meshgrid(h:h:h*ny,h:h:h*nx,h:h:h*nz);

v=zeros(nx,ny,nz);
v[:] .=5000;

s1=zeros(Int64,1,1);
s2=zeros(Int64,1,1);
s3=zeros(Int64,1,1);

s1[:] .=20;
s2[:] .=60;
s3[:] .=70;

s1t=h*s1;
s2t=h*s2;
s3t=h*s3;

r1=zeros(Int64,1,5);
r2=zeros(Int64,1,5);
r3=zeros(Int64,1,5);

r1[:]=[10 20 30 40 50];
r2[:]=[20 30 40 50 60];
r3[:]=[3 3 3 3 3];

r1t=h*r1;
r2t=h*r2;
r3t=h*r3;

path=p3;

write_t=1;
