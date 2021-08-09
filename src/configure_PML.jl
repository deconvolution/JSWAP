"
Configure PML
"
function PML_configuration(nx,ny,nz,dx,dy,dz,lambda,mu,rho,nPML,Rc,lp,PML_active)
    # PML
    vmax=sqrt.((lambda+2*mu) ./rho);
    beta0=(ones(nx,ny,nz) .*vmax .*(nPML+1) .*log(1/Rc)/2/lp/dx);
    beta1=(zeros(nx,ny,nz));
    beta2=copy(beta1);
    beta3=copy(beta1);
    tt=(1:lp)/lp;


    # PML coefficient
    # 6:lp+1
    # nx-lp-5
    beta1=zeros(nx,ny,nz);
    tt=copy(beta1);
    tt[5+2:5+lp+1,:,:]=repeat(reshape((abs.((1:lp) .-lp .-1) ./lp) .^nPML,lp,1,1),1,ny,nz);
    tt[-5+nx-lp:-5+nx-1,:,:]=repeat(reshape(((abs.(nx .-lp .-(nx-lp+1:nx))) ./lp) .^nPML,lp,1,1),1,ny,nz);
    beta1=vmax*(nPML+1)*log(1/Rc)/2/lp/dx.*tt;

    beta2=zeros(nx,ny,nz);
    tt=copy(beta2);
    tt[:,5+2:5+lp+1,:]=repeat(reshape((abs.((1:lp) .-lp .-1) ./lp) .^nPML,1,lp,1),nx,1,nz);
    tt[:,-5+ny-lp:-5+ny-1,:]=repeat(reshape(((abs.(ny .-lp .-(ny-lp+1:ny))) ./lp) .^nPML,1,lp,1),nx,1,nz);
    beta2=vmax*(nPML+1)*log(1/Rc)/2/lp/dy.*tt;

    beta3=zeros(nx,ny,nz);
    tt=copy(beta3);
    tt[:,:,5+2:5+lp+1]=repeat(reshape((abs.((1:lp) .-lp .-1) ./lp) .^nPML,1,1,lp),nx,ny,1);
    tt[:,:,-5+nz-lp:-5+nz-1]=repeat(reshape(((abs.(nz .-lp .-(nz-lp+1:nz))) ./lp) .^nPML,1,1,lp),nx,ny,1);
    beta3=vmax*(nPML+1)*log(1/Rc)/2/lp/dz.*tt;

    for i=1:6
        beta1[i,:,:] .=beta1[7,:,:];
    end

    for i=nx-5:nx
        beta1[i,:,:]=beta1[end-6,:,:];
    end

    for j=1:6
        beta2[:,j,:]=beta2[:,7,:];
    end
    for j=ny-5:ny
        beta2[:,j,:]=beta2[:,end-6,:];
    end

    for k=1:6
        beta3[:,:,k]=beta3[:,:,7];
    end

    for k=nz-5:nz
        beta3[:,:,k]=beta3[:,:,end-6];
    end

    if PML_active[1]==0
        beta1[5+2:5+lp+1,5+lp+2:-5+ny-lp-1,5+lp+2:-5+nz-lp-1] .=0;
    end

    if PML_active[2]==0
        beta1[-5+nz-lp:-5+nz-1,5+lp+2:-5+ny-lp-1,5+lp+2:-5+nz-lp-1] .=0;
    end

    if PML_active[3]==0
        beta2[5+lp+2:-5+nx-lp-1,5+2:5+lp+1,5+lp+2:-5+nz-lp-1] .=0;
    end

    if PML_active[4]==0
        beta2[5+lp+2:-5+nx-lp-1,-5+nz-lp:-5+nz-1,5+lp+2:-5+nz-lp-1] .=0;
    end

    if PML_active[5]==0
        beta3[5+lp+2:-5+nx-lp-1,5+lp+2:-5+ny-lp-1,5+2:5+lp+1] .=0;
    end

    if PML_active[6]==0
        beta3[5+lp+2:-5+nx-lp-1,5+lp+2:-5+ny-lp-1,-5+nz-lp:-5+nz-1] .=0;
    end

    # 3D PML coefficient
    IND=unique(findall(x->x!=0,beta1.*beta2.*beta3));
    IND2=unique(findall(x->x==x,beta1.*beta2+beta2.*beta3+beta3.*beta1));
    IND3=setdiff(IND2,IND);
    beta=beta1+beta2+beta3;
    beta[IND]=beta[IND]/3;
    beta[IND3]=beta[IND3]/2;
    return beta
end
