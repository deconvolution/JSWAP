using JSWAP
nx=50;
ny=100;
nz=80;
h=1;

v=zeros(nx,ny,nz);
for i=1:nz
    v[:,:,i] .=800;
end

v[:,:,1:9] .=340;
##
s1=25;
s2=50;
s3=40;

s1t=h*s1;
s2t=h*s2;
s3t=h*s3;

r1=[2,4,6];
r2=[3,4,5];
r3=[3,3,3];
r1t=h*r1;
r2t=h*r2;
r3t=h*r3;
## compute distance to the source
Y3D,X3D,Z3D=JSWAP.meshgrid(1:ny,1:nx,1:nz);
Y,X,Z=JSWAP.meshgrid(2:ny-1,2:nx-1,2:nz-1);
Z=reshape(Z,(nx-2)*(ny-2)*(nz-2),1);
Y=reshape(Y,(nx-2)*(ny-2)*(nz-2),1);
X=reshape(X,(nx-2)*(ny-2)*(nz-2),1);
dis_s=(X .-s1).^2+(Y .-s2).^2+(Z .-s3).^2;

tt=mapslices(sortperm,dis_s,dims=(1));

X=X[tt];
Y=Y[tt];
Z=Z[tt];

path="./source_1";
##
T,T_obs=JSWAP.eikonal.acoustic_eikonal_forward(nx,ny,nz,h,v,s1,s2,s3,s1t,s2t,s3t,r1,r2,r3,r1t,r2t,r3t,path);
##
plot(heatmap(z=reshape(T[:,s2,:],nx,nz)))
