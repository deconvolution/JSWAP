##
@parallel function compute_sigma_adjoint(dt,dx,dy,dz,lambda,mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p)

    @all(sigmas11)=dt*@all(v1_1_plus)/dx+
    @all(sigmas11)-
    dt*@all(beta).*@all(sigmas11);

    @all(sigmas22)=dt*@all(v2_2_plus)/dy+
    @all(sigmas22)-
    dt*@all(beta).*@all(sigmas22);

    @all(sigmas33)=dt*@all(v3_3_plus)/dz+
    @all(sigmas33)-
    dt*@all(beta).*@all(sigmas33);

    @all(sigmas23)=dt*(@all(v2_3_minus)/dz+@all(v3_2_minus)/dy)+
    @all(sigmas23)-
    dt*@all(beta).*@all(sigmas23);

    @all(sigmas13)=dt*(@all(v1_3_minus)/dz+@all(v3_1_minus)/dx)+
    @all(sigmas13)-
    dt*@all(beta).*@all(sigmas13);

    @all(sigmas12)=dt*(@all(v2_1_minus)/dx+@all(v1_2_minus)/dy)+
    @all(sigmas12)-
    dt*@all(beta).*@all(sigmas12);

    @all(p)=-dt*(@all(v1_1_plus)/dx+@all(v2_2_plus)/dy+@all(v3_3_plus)/dz)+
    @all(p)-
    dt*@all(beta).*@all(p);

    return nothing
end
##
@parallel function compute_auxiliary_in_vadjoint_3D(lambda,mu,
    sigmas11,sigmas22,sigmas33,p,auxiliary_in_vadjoint_3D,
    auxiliary_in_vadjoint_3D2,auxiliary_in_vadjoint_3D3)

    @all(auxiliary_in_vadjoint_3D)=4/3*@all(mu) .*@all(sigmas11)-
    2/3*@all(mu) .*@all(sigmas22)-
    2/3*@all(mu) .*@all(sigmas33)-
    1/3*(3*@all(lambda)+2*@all(mu)) .*@all(p);

    @all(auxiliary_in_vadjoint_3D2)=-2/3*@all(mu) .*@all(sigmas11)+
    4/3*@all(mu) .*@all(sigmas22) -
    2/3*@all(mu) .*@all(sigmas33) -
    1/3*(3*@all(lambda)+2*@all(mu)) .*@all(p);

    @all(auxiliary_in_vadjoint_3D3)=-2/3*@all(mu) .*@all(sigmas11)-
    2/3*@all(mu) .*@all(sigmas22) +
    4/3*@all(mu) .*@all(sigmas33) -
    1/3*(3*@all(lambda)+2*@all(mu)) .*@all(p);

    return nothing
end
##
@parallel function compute_v_adjoint_3D(dt,dx,dy,dz,rho,mu,beta,
    v1,v2,v3,
    sigmas11_1_minus,
    sigmas22_2_minus,
    sigmas33_3_minus,
    sigmas23_2_plus,sigmas23_3_plus,
    sigmas13_1_plus,sigmas13_3_plus,
    sigmas12_1_plus,sigmas12_2_plus,
    p_1_minus,p_2_minus,p_3_minus,
    auxiliary_in_vadjoint_3D_1_minus,
    auxiliary_in_vadjoint_3D2_2_minus,
    auxiliary_in_vadjoint_3D3_3_minus)

    @all(v1)=dt./@all(rho) .*(@all(auxiliary_in_vadjoint_3D_1_minus)/dx+
    @all(mu) .*@all(sigmas12_2_plus)/dy+
    @all(mu) .*@all(sigmas13_3_plus)/dz)+
    @all(v1)-
    dt*@all(beta) .*@all(v1);

    @all(v2)=dt./@all(rho) .*(@all(mu) .*@all(sigmas12_1_plus)/dx+
    @all(auxiliary_in_vadjoint_3D2_2_minus)/dy+
    @all(mu) .*@all(sigmas23_3_plus)/dz)+
    @all(v2)-
    dt*@all(beta) .*@all(v2);

    @all(v3)=dt./@all(rho) .*(@all(mu) .*@all(sigmas13_1_plus)/dx+
    @all(mu) .*@all(sigmas23_2_plus)/dy+
    @all(auxiliary_in_vadjoint_3D3_3_minus))+
    @all(v3)-
    dt*@all(beta) .*@all(v3);

    return nothing
end
##
@timeit ti "iso_3D" function iso_3D_adjoint(dt,dx,dy,dz,nt,
nx,ny,nz,X,Y,Z,r1,r2,r3,s1,s2,s3,src1,src2,src3,srcp,
r1t,r2t,r3t,
s1t,s2t,s3t,
lp,nPML,Rc,pml_active,
C,
plot_interval,
wavefield_interval,
path,
path_pic,
path_model,
path_wavefield);

global data

IND_accuracy=9;

# zero stress condition at the boundaries
C.lambda[1:5,:,:] .=0;
C.lambda[end-4:end,:,:] .=0;
C.lambda[:,1:5,:] .=0;
C.lambda[:,end-4:end,:] .=0;
C.lambda[:,:,1:5] .=0;
C.lambda[:,:,end-4:end] .=0;

C.mu[1:5,:,:] .=0;
C.mu[end-4:end,:,:] .=0;
C.mu[:,1:5,:] .=0;
C.mu[:,end-4:end,:] .=0;
C.mu[:,:,1:5] .=0;
C.mu[:,:,end-4:end] .=0;

# shift coordinate of stiffness

tt=zeros(nx,ny,nz);
#=
tt[1:end-1,1:end-1,1:end-1]=C.lambda[2:end,2:end,2:end];
tt[end,:,:]=C.lambda[end,:,:];
tt[:,end,:]=C.lambda[:,end,:];
tt[:,:,end]=C.lambda[:,:,end];
C.lambda=tt;

tt[1:end-1,1:end-1,1:end-1]=C.rho[2:end,2:end,2:end];
tt[end,:,:]=C.rho[end,:,:];
tt[:,end,:]=C.rho[:,end,:];
tt[:,:,end]=C.rho[:,:,end];
C.rho=tt;
=#
tt=nothing;

d0=Dates.now();
# source number
ns=length(s3);

# create main folder
if isdir(path)==0
    mkdir(path);
end

# create folder for picture
n_picture=1;
n_wavefield=1;
if path_pic!=nothing
    if isdir(path_pic)==0
        mkdir(path_pic);
    end
    # initialize pvd
    pvd=paraview_collection(string(path,"/time_info_adjoint"));
end

# create folder for model
if path_model!=nothing
    if isdir(path_model)==0
        mkdir(path_model)
    end
    vtkfile = vtk_grid(string(path_model,"/material_properties"),X,Y,Z);
    vtkfile["lambda"]=C.lambda;
    vtkfile["mu"]=C.mu;
    vtkfile["rho"]=C.rho;
    vtk_save(vtkfile);
    CSV.write(string(path_model,"/receiver location.csv"),DataFrame([r1t' r2t' r3t'],:auto));
    #CSV.write(string(path_model,"/source location.csv"),DataFrame([s1t' s2t' s3t'],:auto));
end

# create folder for wavefield
if path_wavefield!=nothing
    if isdir(path_wavefield)==0
        mkdir(path_wavefield)
    end
end

# PML
vmax=sqrt.((C.lambda+2*C.mu) ./C.rho);
beta0=(ones(nx,ny,nz) .*vmax .*(nPML+1) .*log(1/Rc)/2/lp/dx);
beta1=(@zeros(nx,ny,nz));
beta2=copy(beta1);
beta3=copy(beta1);
tt=(1:lp)/lp;


# PML coefficient
beta1=@zeros(nx,ny,nz);
tt=copy(beta1);
tt[2:lp+1,:,:]=repeat(reshape((abs.((1:lp) .-lp .-1) ./lp) .^nPML,lp,1,1),1,ny,nz);
tt[nx-lp:nx-1,:,:]=repeat(reshape(((abs.(nx .-lp .-(nx-lp+1:nx))) ./lp) .^nPML,lp,1,1),1,ny,nz);
beta1=vmax*(nPML+1)*log(1/Rc)/2/lp/dx.*tt;

beta2=zeros(nx,ny,nz);
tt=copy(beta2);
tt[:,2:lp+1,:]=repeat(reshape((abs.((1:lp) .-lp .-1) ./lp) .^nPML,1,lp,1),nx,1,nz);
tt[:,ny-lp:ny-1,:]=repeat(reshape(((abs.(ny .-lp .-(ny-lp+1:ny))) ./lp) .^nPML,1,lp,1),nx,1,nz);
beta2=vmax*(nPML+1)*log(1/Rc)/2/lp/dy.*tt;

beta3=zeros(nx,ny,nz);
tt=copy(beta3);
tt[:,:,2:lp+1]=repeat(reshape((abs.((1:lp) .-lp .-1) ./lp) .^nPML,1,1,lp),nx,ny,1);
tt[:,:,nz-lp:nz-1]=repeat(reshape(((abs.(nz .-lp .-(nz-lp+1:nz))) ./lp) .^nPML,1,1,lp),nx,ny,1);
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

if pml_active[1]==0
    beta1[5+2:5+lp+1,5+lp+2:-5+ny-lp-1,5+lp+2:-5+nz-lp-1] .=0;
end

if pml_active[2]==0
    beta1[-5+nz-lp:-5+nz-1,5+lp+2:-5+ny-lp-1,5+lp+2:-5+nz-lp-1] .=0;
end

if pml_active[3]==0
    beta2[5+lp+2:-5+nx-lp-1,5+2:5+lp+1,5+lp+2:-5+nz-lp-1] .=0;
end

if pml_active[4]==0
    beta2[5+lp+2:-5+nx-lp-1,-5+nz-lp:-5+nz-1,5+lp+2:-5+nz-lp-1] .=0;
end

if pml_active[5]==0
    beta3[5+lp+2:-5+nx-lp-1,5+lp+2:-5+ny-lp-1,5+2:5+lp+1] .=0;
end

if pml_active[6]==0
    beta3[5+lp+2:-5+nx-lp-1,5+lp+2:-5+ny-lp-1,-5+nz-lp:-5+nz-1] .=0;
end

# 3D PML coefficient
IND=unique(findall(x->x!=0,beta1.*beta2.*beta3));
IND2=unique(findall(x->x==x,beta1.*beta2+beta2.*beta3+beta3.*beta1));
IND3=setdiff(IND2,IND);
beta=beta1+beta2+beta3;
beta[IND]=beta[IND]/3;
beta[IND3]=beta[IND3]/2;

vmax=beta01=beta02=beta03=tt=beta1=beta2=beta3=IND=IND2=IND3=nothing;

# receiver configuration
R1=@zeros(nt,length(r3));
R2=copy(R1);
R3=copy(R1);
P=@zeros(nt,length(r3));

# wave vector
v1=@zeros(nx,ny,nz);
v2=copy(v1);
v3=copy(v1);

sigmas11=copy(v1);
sigmas22=copy(v1);
sigmas33=copy(v1);
sigmas23=copy(v1);
sigmas13=copy(v1);
sigmas12=copy(v1);
p=copy(v1);

v1_1_plus=copy(v1);
v1_2_minus=copy(v1);
v1_3_minus=copy(v1);
v2_1_minus=copy(v1);
v2_2_plus=copy(v1);
v2_3_minus=copy(v1);
v3_1_minus=copy(v1);
v3_2_minus=copy(v1);
v3_3_plus=copy(v1);

auxiliary_in_vadjoint_3D=copy(v1);
auxiliary_in_vadjoint_3D2=copy(v1);
auxiliary_in_vadjoint_3D3=copy(v1);

dtt1=@zeros(nx-IND_accuracy,ny,nz);
dtt2=@zeros(nx,ny-IND_accuracy,nz);
dtt3=@zeros(nx,ny,nz-IND_accuracy);

sigmas11_1_minus=copy(v1);
sigmas22_2_minus=copy(v1);
sigmas33_3_minus=copy(v1);
sigmas23_2_plus=copy(v1);
sigmas23_3_plus=copy(v1);
sigmas13_1_plus=copy(v1);
sigmas13_3_plus=copy(v1);
sigmas12_1_plus=copy(v1);
sigmas12_2_plus=copy(v1);
p_1_minus=copy(v1);
p_2_minus=copy(v1);
p_3_minus=copy(v1);


auxiliary_in_vadjoint_3D_1_minus=copy(v1);
auxiliary_in_vadjoint_3D2_2_minus=copy(v1);
auxiliary_in_vadjoint_3D3_3_minus=copy(v1);


l=1;
# save wavefield
if path_wavefield!=nothing && wavefield_interval!=0
    if mod(l,wavefield_interval)==0
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/v1_",n_wavefield,".mat"),data);

        write2mat(string(path_wavefield,"/v2_",n_wavefield,".mat"),data);

        write2mat(string(path_wavefield,"/v3_",n_wavefield,".mat"),data);

        write2mat(string(path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
        write2mat(string(path_wavefield,"/sigmas22_",n_wavefield,".mat"),data);
        write2mat(string(path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
        write2mat(string(path_wavefield,"/sigmas23_",n_wavefield,".mat"),data);
        write2mat(string(path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
        write2mat(string(path_wavefield,"/sigmas12_",n_wavefield,".mat"),data);
        write2mat(string(path_wavefield,"/p_",n_wavefield,".mat"),data);
        n_wavefield=n_wavefield+1;
    end
end
#
pro_bar=Progress(nt,1,"adjoint_simulation...",50);
for l=1:nt-1
    @parallel Dx_inn(v1,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_plus(dtt1,v1_1_plus);

    @parallel Dy_inn(v1,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,v1_2_minus);

    @parallel Dz_inn(v1,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,v1_3_minus);

    @parallel Dx_inn(v2,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,v2_1_minus);

    @parallel Dy_inn(v2,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_plus(dtt2,v2_2_plus);

    @parallel Dz_inn(v2,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,v2_3_minus);

    @parallel Dx_inn(v3,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,v3_1_minus);

    @parallel Dy_inn(v3,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,v3_2_minus);

    @parallel Dz_inn(v3,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_plus(dtt3,v3_3_plus);
    @timeit ti "compute_sigma" @parallel compute_sigma_adjoint(dt,dx,dy,dz,C.lambda,C.mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p);

    @timeit ti "minus" @parallel compute_auxiliary_in_vadjoint_3D(C.lambda,C.mu,
    sigmas11,sigmas22,sigmas33,p,auxiliary_in_vadjoint_3D,
    auxiliary_in_vadjoint_3D2,auxiliary_in_vadjoint_3D3);

    @parallel Dx_inn(sigmas11,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,sigmas11_1_minus);

    @parallel Dy_inn(sigmas22,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,sigmas22_2_minus);

    @parallel Dz_inn(sigmas33,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,sigmas33_3_minus);

    @parallel Dy_inn(sigmas23,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_plus(dtt2,sigmas23_2_plus);
    @parallel Dz_inn(sigmas23,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_plus(dtt3,sigmas23_3_plus);

    @parallel Dx_inn(sigmas13,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_plus(dtt1,sigmas13_1_plus);
    @parallel Dz_inn(sigmas13,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_plus(dtt3,sigmas13_3_plus);

    @parallel Dx_inn(sigmas12,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_plus(dtt1,sigmas12_1_plus);
    @parallel Dy_inn(sigmas12,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_plus(dtt2,sigmas12_2_plus);

    @parallel Dx_inn(p,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,p_1_minus);
    @parallel Dy_inn(p,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,p_2_minus);
    @parallel Dz_inn(p,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,p_3_minus);

    @parallel Dx_inn(auxiliary_in_vadjoint_3D,dtt1);
    @parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,auxiliary_in_vadjoint_3D_1_minus);


    @parallel Dy_inn(auxiliary_in_vadjoint_3D2,dtt2);
    @parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,auxiliary_in_vadjoint_3D2_2_minus);

    @parallel Dz_inn(auxiliary_in_vadjoint_3D3,dtt3);
    @parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,auxiliary_in_vadjoint_3D3_3_minus);



    @timeit ti "compute_v" @parallel compute_v_adjoint_3D(dt,dx,dy,dz,C.rho,C.mu,beta,
    v1,v2,v3,
    sigmas11_1_minus,
    sigmas22_2_minus,
    sigmas33_3_minus,
    sigmas23_2_plus,sigmas23_3_plus,
    sigmas13_1_plus,sigmas13_3_plus,
    sigmas12_1_plus,sigmas12_2_plus,
    p_1_minus,p_2_minus,p_3_minus,
    auxiliary_in_vadjoint_3D_1_minus,
    auxiliary_in_vadjoint_3D2_2_minus,
    auxiliary_in_vadjoint_3D3_3_minus);


    if ns==1
        v1[CartesianIndex.(s1,s2,s3)]=v1[CartesianIndex.(s1,s2,s3)]+1 ./C.rho[CartesianIndex.(s1,s2,s3)] .*src1[l];
        v2[CartesianIndex.(s1,s2,s3)]=v2[CartesianIndex.(s1,s2,s3)]+1 ./C.rho[CartesianIndex.(s1,s2,s3)] .*src2[l];
        v3[CartesianIndex.(s1,s2,s3)]=v3[CartesianIndex.(s1,s2,s3)]+1 ./C.rho[CartesianIndex.(s1,s2,s3)] .*src3[l];
        p[CartesianIndex.(s1,s2,s3)]=p[CartesianIndex.(s1,s2,s3)]+srcp[l];
    else
        v1[CartesianIndex.(s1,s2,s3)]=v1[CartesianIndex.(s1,s2,s3)]+1 ./C.rho[CartesianIndex.(s1,s2,s3)] .*src1[l,:]';
        v2[CartesianIndex.(s1,s2,s3)]=v2[CartesianIndex.(s1,s2,s3)]+1 ./C.rho[CartesianIndex.(s1,s2,s3)] .*src2[l,:]';
        v3[CartesianIndex.(s1,s2,s3)]=v3[CartesianIndex.(s1,s2,s3)]+1 ./C.rho[CartesianIndex.(s1,s2,s3)] .*src3[l,:]';
        p[CartesianIndex.(s1,s2,s3)]=p[CartesianIndex.(s1,s2,s3)]+srcp[l,:]';

    end

    # assign recordings
    @timeit ti "receiver" R1[l+1,:]=reshape(v1[CartesianIndex.(r1,r2,r3)],length(r3),);
    @timeit ti "receiver" R2[l+1,:]=reshape(v2[CartesianIndex.(r1,r2,r3)],length(r3),);
    @timeit ti "receiver" R3[l+1,:]=reshape(v3[CartesianIndex.(r1,r2,r3)],length(r3),);
    @timeit ti "receiver" P[l+1,:]=reshape(p[CartesianIndex.(r1,r2,r3)],length(r3),);
    # save wavefield
    if path_wavefield!=nothing && wavefield_interval!=0
        if mod(l,wavefield_interval)==0
            data=v1;
            write2mat(string(path_wavefield,"/v1_",n_wavefield,".mat"),data);
            data=v2;
            write2mat(string(path_wavefield,"/v2_",n_wavefield,".mat"),data);
            data=v3;
            write2mat(string(path_wavefield,"/v3_",n_wavefield,".mat"),data);
            data=sigmas11;
            write2mat(string(path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
            data=sigmas22;
            write2mat(string(path_wavefield,"/sigmas22_",n_wavefield,".mat"),data);
            data=sigmas33;
            write2mat(string(path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
            data=sigmas23;
            write2mat(string(path_wavefield,"/sigmas23_",n_wavefield,".mat"),data);
            data=sigmas13;
            write2mat(string(path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
            data=sigmas12;
            write2mat(string(path_wavefield,"/sigmas12_",n_wavefield,".mat"),data);
            data=p;
            write2mat(string(path_wavefield,"/p_",n_wavefield,".mat"),data);
            n_wavefield=n_wavefield+1;
        end
    end

    # plot
    if path_pic!=nothing && plot_interval!=0
        if mod(l,plot_interval)==0 || l==nt-1
            vtkfile = vtk_grid(string(path_pic,"/wavefield_pic_",n_picture),X,Y,Z);
            vtkfile["v1"]=v1;
            vtkfile["v2"]=v2;
            vtkfile["v3"]=v3;
            vtkfile["p"]=p;
            vtkfile["sigmas33"]=sigmas33;
            vtkfile["lambda"]=C.lambda;
            vtkfile["mu"]=C.mu;
            vtkfile["rho"]=C.rho;
            pvd[dt*(l+1)]=vtkfile;
            n_picture=n_picture+1;
        end
    end

    next!(pro_bar);
end
#=
R1=R1 .*Rm[:,:,1];
R2=R2 .*Rm[:,:,2];
R3=R3 .*Rm[:,:,2];
P=P .*Rm[:,:,3];

data=R1;
write2mat(string(path_rec,"/rec_1.mat"),data);
data=R2;
write2mat(string(path_rec,"/rec_2.mat"),data);
data=R3;
write2mat(string(path_rec,"/rec_3.mat"),data);
data=P;
write2mat(string(path_rec,"/rec_p.mat"),data);
=#
if path_pic!=nothing && plot_interval!=0
    vtk_save(pvd);
end

return v1,v2,v3
end
## correlate
function compute_sensitivity_kernel_iso_3D(path,nt,nx,ny,nz,dt,dx,dy,dz,X,Y,Z,lambda,mu,
    path_forward_wavefield,path_adjoint_wavefield,wavefield_interval)
    if isdir(string(path,"/sensitivity_kernel"))==0
        mkdir(string(path,"/sensitivity_kernel"));
    end

    if isdir(string(path,"/sensitivity_kernel/pic"))==0
        mkdir(string(path,"/sensitivity_kernel/pic"));
    end
    n_picture=1;
    pvd=paraview_collection(string(path,"/sensitivity_kernel/time_info"));

    n_wavefield=convert(Int32,floor(nt/wavefield_interval))-1;
    gradient_lambda=@zeros(nx,ny,nz);
    gradient_mu=@zeros(nx,ny,nz);
    klambda=@zeros(nx,ny,nz);
    kmu=@zeros(nx,ny,nz);

    pro_bar=Progress(n_wavefield,1,"compute_gradient...",50);

    for l2=2:n_wavefield
        data=readmat(string(path_forward_wavefield,"/v1_",l2,".mat"),"data");
        v1f=data;
        data=readmat(string(path_forward_wavefield,"/v2_",l2,".mat"),"data");
        v2f=data;
        data=readmat(string(path_forward_wavefield,"/v3_",l2,".mat"),"data");
        v3f=data;

        # adjoint
        data=readmat(string(path_adjoint_wavefield,"/sigmas11_",n_wavefield-l2+1,".mat"),"data");
        sigmas11a=data;
        data=readmat(string(path_adjoint_wavefield,"/sigmas22_",n_wavefield-l2+1,".mat"),"data");
        sigmas22a=data;
        data=readmat(string(path_adjoint_wavefield,"/sigmas33_",n_wavefield-l2+1,".mat"),"data");
        sigmas33a=data;
        data=readmat(string(path_adjoint_wavefield,"/sigmas23_",n_wavefield-l2+1,".mat"),"data");
        sigmas23a=data;
        data=readmat(string(path_adjoint_wavefield,"/sigmas13_",n_wavefield-l2+1,".mat"),"data");
        sigmas13a=data;
        data=readmat(string(path_adjoint_wavefield,"/sigmas12_",n_wavefield-l2+1,".mat"),"data");
        sigmas12a=data;
        data=readmat(string(path_adjoint_wavefield,"/p_",n_wavefield-l2+1,".mat"),"data");
        pa=data;
        @parallel correlate_forward_adjoint(klambda,kmu,dt,dx,dy,dz,v1f,v2f,v3f,sigmas11a,sigmas22a,sigmas33a,
        sigmas23a,sigmas13a,sigmas12a,pa);
        @parallel compute_gradient(gradient_lambda,gradient_mu,klambda,kmu);

        vtkfile = vtk_grid(string(path,"/sensitivity_kernel/pic/sensitivity",n_picture),X,Y,Z);
        vtkfile["lambda"]=C.lambda;
        vtkfile["mu"]=C.mu;
        vtkfile["rho"]=C.rho;
        vtkfile["klambda"]=klambda;
        vtkfile["kmu"]=kmu;
        vtkfile["gradient_lambda"]=gradient_lambda/n_picture;
        vtkfile["gradient_mu"]=gradient_mu/n_picture;
        pvd[dt*wavefield_interval*l2]=vtkfile;
        vtk_save(vtkfile);
        n_picture=n_picture+1;
        next!(pro_bar);
    end
    vtk_save(pvd);
    return nothing
end

@parallel function correlate_forward_adjoint(klambda,kmu,dt,dx,dy,dz,v1f,v2f,v3f,sigmas11a,sigmas22a,sigmas33a,
    sigmas23a,sigmas13a,sigmas12a,pa)
    @inn(klambda)=@inn(pa) .*(@d_xi(v1f)/dx+@d_yi(v2f)/dy+@d_zi(v3f)/dz);
    @inn(kmu)=@inn(pa).*(2/3*@d_xi(v1f)/dx+2/3*@d_yi(v2f)/dy+2/3*@d_zi(v3f)/dz)-
    @inn(sigmas13a) .*(@d_xi(v3f)/dx+@d_zi(v1f)/dz)-
    @inn(sigmas23a) .*(@d_yi(v3f)/dy+@d_zi(v2f)/dz)-
    @inn(sigmas12a) .*(@d_xi(v2f)/dx+@d_yi(v1f)/dy)+
    @inn(sigmas11a) .*(-4/3*@d_xi(v1f)/dx+2/3*@d_yi(v2f)/dy+2/3*@d_zi(v3f)/dz)+
    @inn(sigmas22a) .*(2/3*@d_xi(v1f)/dx-4/3*@d_yi(v2f)/dy+2/3*@d_zi(v3f)/dz)+
    @inn(sigmas33a) .*(2/3*@d_xi(v1f)/dx+2/3*@d_yi(v2f)/dy-4/3*@d_zi(v3f)/dz);
    return nothing
end

@parallel function compute_gradient(gradient_lambda,gradient_mu,klambda,kmu)
    @inn(gradient_lambda)=@inn(gradient_lambda)+@inn(klambda);
    @inn(gradient_mu)=@inn(gradient_mu)+@inn(kmu);
    return nothing
end
