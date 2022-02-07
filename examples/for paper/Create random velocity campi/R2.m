function [FF,k3]=R2(k2,a,kx,ky,kz)
[ky2,kx2,kz2]=meshgrid(ky,kx,kz);
FF=zeros(size(kx2));
k=sqrt(kx2.^2+ky2.^2+kz2.^2);

nx=length(kx);
ny=length(ky);
nz=length(kz);

ex=0;
ey=0;
ez=0;
if mod(nx,2)~=0
    nx=nx+1;
    ex=1;
    if size(k2)~=[1,1,1]
    tk2=zeros(nx,ny,nz);
    tk2(1:floor(nx/2),:,:)=k2(1:floor(nx/2),:,:);
    tk2(floor(nx/2)+2:nx,:,:)=k2(floor(nx/2)+1:nx-1,:,:);
    tk2(floor(nx/2)+1,:,:)=k2(floor(nx/2),:,:);
    k2=tk2;
    end
end

if mod(ny,2)~=0
    ny=ny+1;
    ey=1;
    if size(k2)~=[1,1,1]
    tk2=zeros(nx,ny,nz);
    tk2(:,1:floor(ny/2),:)=k2(:,1:floor(ny/2),:);
    tk2(:,floor(ny/2)+2:ny,:)=k2(:,floor(ny/2)+1:ny-1,:);
    tk2(:,floor(ny/2)+1,:)=k2(:,floor(ny/2),:);
    k2=tk2;
    end
end

if mod(nz,2)~=0
    nz=nz+1;
    ez=1;
    if size(k2)~=[1,1,1]
    tk2=zeros(nx,ny,nz);
    tk2(:,:,1:floor(nz/2))=k2(:,:,1:floor(nz/2));
    tk2(:,:,floor(nz/2)+2:nz)=k2(:,:,floor(nz/2)+1:nz-1);
    tk2(:,:,floor(nz/2)+1)=k2(:,:,floor(nz/2));
    k2=tk2;
    end
end

k3=zeros(nx,ny,nz);

k3(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2))=k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2));

k3(floor(nx/2)+1:end,1:floor(ny/2),1:floor(nz/2))=flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),1);
k3(1:floor(nx/2),floor(ny/2)+1:end,1:floor(nz/2))=flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),2);
k3(1:floor(nx/2),1:floor(ny/2),floor(nz/2)+1:end)=flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),3);

k3(1:floor(nx/2),floor(ny/2)+1:end,floor(nz/2)+1:end)=flip(flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),2),3);
k3(floor(nx/2)+1:end,1:floor(ny/2),floor(nz/2)+1:end)=flip(flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),1),3);
k3(floor(nx/2)+1:end,floor(ny/2)+1:end,1:floor(nz/2))=flip(flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),1),2);

k3(floor(nx/2)+1:end,floor(ny/2)+1:end,floor(nz/2)+1:end)=flip(flip(flip(k(1:floor(nx/2),1:floor(ny/2),1:floor(nz/2)),1),2),3);
%%

% von karman

d=3;
N=-.2;
FF=k2.*(a.^-2+k3.^2).^(-d/4-N/2);


% gaussian
%{
FF=k2.*exp(-a .^2.*k3.^2/8);
%}
% exponential
%{
FF=k2.*(a.^-2+k3.^2).^(-(d+1)/4);
%}
%%
if ex==1
    FF(nx/2,:,:)=[];
    k3(nx/2,:,:)=[];
end

if ey==1
    FF(:,ny/2,:)=[];
    k3(:,ny/2,:)=[];
end
if ez==1
    FF(:,:,nz/2)=[];
    k3(:,:,nz/2)=[];
end
end