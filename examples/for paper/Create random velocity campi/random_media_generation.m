function [v,k3]=random_media_generation(nx,ny,nz,dx,dy,dz,a,kappa,nkx,nky,nkz)

rng(1)
W=randn(nx,ny,nz);

% sampling interval
ksx=1/dx;
ksy=1/dy;
ksz=1/dz;


kx=ksx*(1:nkx)/nkx*2*pi;
ky=ksy*(1:nky)/nky*2*pi;
kz=ksz*(1:nkz)/nkz*2*pi;


FW=fftn(W,[nkx,nky,nkz]);
% spectral filter
[FF,k3]=R2(kappa,a,kx,ky,kz);
Fv=FF.*FW;
% inverse transform of velocity in wavenumber domain
v=real(ifftn(Fv,[nkx,nky,nkz]));
end