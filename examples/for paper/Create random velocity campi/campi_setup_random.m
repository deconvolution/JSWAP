%% model set up for event
%{
This is a routine of creating random media based on the velocity model of Campi Flegrei.

Additional packages:
deg2utm, https://de.mathworks.com/matlabcentral/fileexchange/10915-deg2utm
vtkwrite, https://de.mathworks.com/matlabcentral/fileexchange/47814-vtkwrite-exports-various-2d-3d-data-to-paraview-in-vtk-file-format
%}
%% Dimensions
close all;
clear all;
%% add grids
nx_a_minus=20;
nx_a_plus=15;
ny_a_minus=15;
ny_a_plus=15;
nz_a_plus=15;
nz_a_minus=15;
%% read velocity
%{
% The velocity model is obtained from Battaglia et al. (2008).

Battaglia, J., Zollo, A., Virieux, J., & Dello Iacono, D. (2008). Merging active
570 and passive data sets in traveltime tomography: The case study of campi
571 flegrei caldera (southern italy) [Journal Article]. Geophysical Prospecting,
572 56 (4), 555-573. 
%}
M=table2array(readtable('./modvPS.txt'));

X0=unique(M(:,1));
Y0=flip(unique(M(:,2)),1);
Z0=unique(M(:,3));

nx0=length(X0);
ny0=length(Y0);
nz0=length(Z0);


X1=flip(permute(reshape(M(:,1),[nz0,ny0,nx0]),[3,2,1]),2);
Y1=flip(permute(reshape(M(:,2),[nz0,ny0,nx0]),[3,2,1]),2);
Z1=flip(permute(reshape(M(:,3),[nz0,ny0,nx0]),[3,2,1]),2);

X2=unique(X1);
Y2=unique(Y1);
Z2=unique(Z1);

nx=nx0+nx_a_minus+nx_a_plus;
ny=ny0+ny_a_minus+ny_a_plus;
nz=nz0+nz_a_minus+nz_a_plus;

X=zeros(nx,ny,nz);
Y=zeros(nx,ny,nz);
Z=zeros(nx,ny,nz);

dX=X2(2)-X2(1);
dY=Y2(2)-Y2(1);
dZ=Z2(2)-Z2(1);

X0=zeros(nx,1);
Y0=zeros(ny,1);
Z0=zeros(nz,1);

X0(:)=min(X2)-nx_a_minus*dX:dX:max(X2)+nx_a_plus*dX;
Y0(:)=min(Y2)-ny_a_minus*dY:dY:max(Y2)+ny_a_plus*dY;
Z0(:)=max(Z2)+nz_a_plus*dZ:-dZ:min(Z2)-nz_a_minus*dZ;

[Y,X,Z]=meshgrid(Y0,X0,Z0);

vp2=flip(permute(reshape(M(:,4),[nz0,ny0,nx0]),[3,2,1]),2);
vs2=flip(permute(reshape(M(:,5),[nz0,ny0,nx0]),[3,2,1]),2);

vp=zeros(nx,ny,nz);
vs=vp;
rho=vp;

vp(nx_a_minus+1:nx-nx_a_plus,ny_a_minus+1:ny-ny_a_plus,nz_a_minus+1:nz-nz_a_plus)=vp2;
vs(nx_a_minus+1:nx-nx_a_plus,ny_a_minus+1:ny-ny_a_plus,nz_a_minus+1:nz-nz_a_plus)=vs2;

% copy in x direction vp
vp(1:nx_a_minus,:,:)=repmat(vp(nx_a_minus+1,:,:),[nx_a_minus,1,1]);
vp(nx-nx_a_plus+1:end,:,:)=repmat(vp(nx-nx_a_plus,:,:),[nx_a_plus,1,1]);

% copy in y direction vp
vp(:,1:ny_a_minus,:)=repmat(vp(:,ny_a_minus+1,:),[1,ny_a_minus,1]);
vp(:,ny-ny_a_plus+1:end,:)=repmat(vp(:,ny-ny_a_plus,:),[1,ny_a_plus,1]);

% copy in z direction vp
vp(:,:,1:nz_a_minus)=repmat(vp(:,:,nz_a_minus+1),[1,1,nz_a_minus]);
vp(:,:,nz-nz_a_plus+1:end)=repmat(vp(:,:,nz-nz_a_plus),[1,1,nz_a_plus]);

% copy in x direction vs
vs(1:nx_a_minus,:,:)=repmat(vs(nx_a_minus+1,:,:),[nx_a_minus,1,1]);
vs(nx-nx_a_plus+1:end,:,:)=repmat(vs(nx-nx_a_plus,:,:),[nx_a_plus,1,1]);

% copy in y direction vs
vs(:,1:ny_a_minus,:)=repmat(vs(:,ny_a_minus+1,:),[1,ny_a_minus,1]);
vs(:,ny-ny_a_plus+1:end,:)=repmat(vs(:,ny-ny_a_plus,:),[1,ny_a_plus,1]);

% copy in z direction vs
vs(:,:,1:nz_a_minus)=repmat(vs(:,:,nz_a_minus+1),[1,1,nz_a_minus]);
vs(:,:,nz-nz_a_plus+1:end)=repmat(vs(:,:,nz-nz_a_plus),[1,1,nz_a_plus]);
%% read topography
M=load('./DTM_LiDAR_5m_CF.dat');
%% save velocity field before fluctuated
tit={'WE','SN','Al','vp','vs','vp/vs'};
vtkwrite('campi_vp.vtk','structured_grid',X,Y,Z,'scalars','campi_vp',vp);
vtkwrite('campi_vs.vtk','structured_grid',X,Y,Z,'scalars','campi_vs',vs);
%% Dimensions
dx=100;
dy=100;
dz=100;
%% random media
% vp
kappa=.0002*vp;
a=1000;

nkx=nx;
nky=ny;
nkz=nz;

[v,~]=random_media_generation(nx,ny,nz,dx,dy,dz,a,kappa,nkx,nky,nkz);
vtkwrite('vp_fluctuated.vtk','structured_grid',X,Y,Z,'scalars','fluctuation',v);
vp=vp+v;
figure; imagesc(real(v(1:nx,1:ny,30)));

% vs
kappa=.0002*vs;
a=1000;

[v,~]=random_media_generation(nx,ny,nz,dx,dy,dz,a,kappa,nkx,nky,nkz);

vs=vs+v;
%% set air
vp(:,:,1:nz_a_minus)=0;
vs(:,:,1:nz_a_minus)=0;
%% set topography
X2=unique(M(:,1));
Y2=unique(M(:,2));

nx2=length(X2);
ny2=length(Y2);

M2=reshape(M(:,3),[nx2,ny2]);
[Y3,X3]=meshgrid(Y2,X2);

topo=interp2(Y3,X3,M2,Y(:,:,1),X(:,:,1),'nearest');

topo(isnan(topo))=0;

topo_amplitude=topo;

topo=fix(topo/dz);
% 
C.rho=310*vp.^.25;
C.rho(C.rho==0)=2000;
mu=C.rho.*vs.^2;
lambda=C.rho.*vp.^2-2*mu;
topo_Z=zeros(nx,ny);
topo_Z=topo_amplitude+Z(1,1,(Z(1,1,:)==0));
%{
for i=1:nx
    for j=1:ny
        [~,tt]=max(find(mu(i,j,:)==0));
        if tt<-topo(i,j)+find(Z(1,1,:)==0)
            lambda(i,j,1:(-topo(i,j)+find(Z(1,1,:)==0)))=0;
            mu(i,j,1:(-topo(i,j)+find(Z(1,1,:)==0)))=0;
        else
            lambda(i,j,-topo(i,j)+find(Z(1,1,:)==0):tt)=lambda(i,j,tt+1);
            mu(i,j,-topo(i,j)+find(Z(1,1,:)==0):tt)=mu(i,j,tt+1);
            
        end
    end
end
%}
%% initialize density and stiffness
C.C11=zeros(nx,ny,nz);
C.C12=C.C11;
C.C13=C.C11;
C.C14=C.C11;
C.C15=C.C11;
C.C16=C.C11;
C.C22=C.C11;
C.C23=C.C11;
C.C24=C.C11;
C.C25=C.C11;
C.C26=C.C11;
C.C33=C.C11;
C.C34=C.C11;
C.C35=C.C11;
C.C36=C.C11;
C.C44=C.C11;
C.C45=C.C11;
C.C46=C.C11;
C.C55=C.C11;
C.C56=C.C11;
C.C66=C.C11;
%% transfer stiffness, density and viscosity
C.C11=lambda+2*mu;
C.C12=lambda;
C.C13=lambda;
C.C14=0;
C.C15=0;
C.C16=0;
C.C22=lambda+2*mu;
C.C23=lambda;
C.C24=0;
C.C25=0;
C.C26=0;
C.C33=lambda+2*mu;
C.C34=0;
C.C35=0;
C.C36=0;
C.C44=mu;
C.C45=0;
C.C46=0;
C.C55=mu;
C.C56=0;
C.C66=mu;
C.lambda2=0;
C.mu2=0;
%% receiver location
filename = 'stazioni_151007.txt';
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
stazioni=A.data;
namest=A.textdata;
coor=stazioni(:,1)+stazioni(:,2)/60;
b=stazioni(:,3)+stazioni(:,4)/60;
latlong=[coor,b];
[Er,Nr,Zone]=deg2utm(latlong(:,1),latlong(:,2));
r1=fix((Er'-min(X0))/dx);
r2=fix((Nr'-min(Y0))/dy);
r3=fix((max(Z(:))-stazioni(:,5)'*1000)/dz);
%{
for ir=1:length(r1)
    [~,tt]=max(find(C.C55(r1(ir),r2(ir),:)==0));
    if r3(ir)<=tt
        r3(ir)=tt+1;
    end
end
%}
%% Source and source signals
[Es,Ns,~]=deg2utm(40+49.50/60,14+9.02/60);
s1=fix((Es'-min(X0))/dx);
s2=fix((Ns'-min(Y0))/dy);
s3=fix((max(Z(:))-(-1.53)*1000)/dz);
%% save model setup
campi_model.X=X;
campi_model.Y=Y;
campi_model.Z=Z;
campi_model.nx=nx;
campi_model.ny=ny;
campi_model.nz=nz;
campi_model.dx=dx;
campi_model.dy=dy;
campi_model.dz=dz;
campi_model.C=C;
campi_model.s1=s1;
campi_model.s2=s2;
campi_model.s3=s3;
campi_model.r1=r1;
campi_model.r2=r2;
campi_model.r3=r3;
save('campi_model_random.mat','campi_model');
vtkwrite('campi_model_random.vtk','structured_grid',X,Y,Z,'scalars','C44',C.C44);