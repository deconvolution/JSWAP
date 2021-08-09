module JSWAP_CPU_3D
export JSWAP_CPU_3D_isotropic_solver
## Using ParallelStencil
include("./ParallelStencil/ParallelStencil.jl");
using Random,MAT,Plots,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,CSV,
.ParallelStencil,.ParallelStencil.FiniteDifferences3D
## Use CPU for ParallelStencil
const USE_GPU=false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3);
else
    @init_parallel_stencil(Threads, Float64, 3);
end
## timing
ti=TimerOutput();
##
"
Compuites auxiliary variable used to compute sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function compute_ax(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt)

    @all(Ax)=@all(ax);
    @all(Ax2)=@all(ax2);
    @all(Ax3)=@all(ax3);
    @all(Ax4)=@all(ax4);
    @all(Ax5)=@all(ax5);
    @all(Ax6)=@all(ax6);
    @all(Ax7)=@all(ax7);

    @all(ax)=4*@all(mu) .*@all(v1_1_plus)/dx+
    (-2*@all(mu)) .*@all(v2_2_plus)/dy+
    (-2*@all(mu)) .*@all(v3_3_plus)/dz;

    @all(ax2)=(-2*@all(mu)) .*@all(v1_1_plus)/dx+
    (4*@all(mu)) .*@all(v2_2_plus)/dy+
    (-2*@all(mu)) .*@all(v3_3_plus)/dz;

    @all(ax3)=(-2*@all(mu)) .*@all(v1_1_plus)/dx+
    (-2*@all(mu)) .*@all(v2_2_plus)/dy+
    (4*@all(mu)) .*@all(v3_3_plus)/dz;

    @all(ax4)=@all(mu).*@all(v2_1_minus)/dx+
    @all(mu).*@all(v1_2_minus)/dy;

    @all(ax5)=@all(mu) .*@all(v3_1_minus)/dx+
    @all(mu) .*@all(v1_3_minus)/dz;

    @all(ax6)=@all(mu).*@all(v3_2_minus)/dy+
    @all(mu).*@all(v2_3_minus)/dz;

    @all(ax7)=(3*@all(lambda)+2*@all(mu)) .*@all(v1_1_plus)/dx+
    (3*@all(lambda)+2*@all(mu)) .*@all(v2_2_plus)/dy+
    (3*@all(lambda)+2*@all(mu)) .*@all(v3_3_plus)/dz;

    @all(ax_dt)=(@all(ax)-@all(Ax))/dt;
    @all(ax2_dt)=(@all(ax2)-@all(Ax2))/dt;
    @all(ax3_dt)=(@all(ax3)-@all(Ax3))/dt;
    @all(ax4_dt)=(@all(ax4)-@all(Ax4))/dt;
    @all(ax5_dt)=(@all(ax5)-@all(Ax5))/dt;
    @all(ax6_dt)=(@all(ax6)-@all(Ax6))/dt;
    @all(ax7_dt)=(@all(ax7)-@all(Ax7))/dt;
    return nothing
end
##

"
Computes sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function compute_sigma(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt)
    @all(sigmas11)=1/3*dt*(
    @all(ax)+@all(inv_Qa) .*@all(ax_dt))+
    @all(sigmas11)-
    dt*@all(beta).*@all(sigmas11);

    @all(sigmas22)=1/3*dt*(
    @all(ax2)+@all(inv_Qa) .*@all(ax2_dt))+
    @all(sigmas22)-
    dt*@all(beta).*@all(sigmas22);

    @all(sigmas33)=1/3*dt*(
    @all(ax3)+@all(inv_Qa) .*@all(ax3_dt))+
    @all(sigmas33)-
    dt*@all(beta).*@all(sigmas33);

    @all(sigmas12)=dt*(
    @all(ax4)+@all(inv_Qa) .*@all(ax4_dt))+
    @all(sigmas12)-
    dt*@all(beta).*@all(sigmas12);

    @all(sigmas13)=dt*(
    @all(ax5)+@all(inv_Qa) .*@all(ax5_dt))+
    @all(sigmas13)-
    dt*@all(beta).*@all(sigmas13);

    @all(sigmas23)=dt*(
    @all(ax6)+@all(inv_Qa) .*@all(ax6_dt))+
    @all(sigmas23)-
    dt*@all(beta).*@all(sigmas23);

    @all(p)=-1/3*dt*(
    @all(ax7)+@all(inv_Qa) .*@all(ax7_dt))+
    @all(p)-
    dt*@all(beta).*@all(p);

    return nothing
end
##

"
Computes v, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function compute_v(dt,dx,dy,dz,rho,beta,
    v1,v2,v3,
    sigmas11_1_minus,
    sigmas22_2_minus,
    sigmas33_3_minus,
    sigmas23_2_plus,sigmas23_3_plus,
    sigmas13_1_plus,sigmas13_3_plus,
    sigmas12_1_plus,sigmas12_2_plus,
    p_1_minus,p_2_minus,p_3_minus)

    @all(v1)=dt./@all(rho) .*((@all(sigmas11_1_minus)-@all(p_1_minus))/dx+
    @all(sigmas12_2_plus)/dy+
    @all(sigmas13_3_plus)/dz)+
    @all(v1)-
    dt*@all(beta) .*@all(v1);

    @all(v2)=dt./@all(rho) .*(@all(sigmas12_1_plus)/dx+
    (@all(sigmas22_2_minus)-@all(p_2_minus))/dy+
    @all(sigmas23_3_plus)/dz)+
    @all(v2)-
    dt*@all(beta) .*@all(v2);

    @all(v3)=dt./@all(rho) .*(@all(sigmas13_1_plus)/dx+
    @all(sigmas23_2_plus)/dy+
    (@all(sigmas33_3_minus)-@all(p_3_minus))/dz)+
    @all(v3)-
    dt*@all(beta) .*@all(v3);

    return nothing
end
##

"
Finite-difference solver for isotropic viscoelastic media. JSWAP_CPU_3D_isotropic_solver
"

@timeit ti "iso_3D" function JSWAP_CPU_3D_isotropic_solver(dt,dx,dy,dz,nt,
nx,ny,nz,X,Y,Z,r1,r2,r3,s1,s2,s3,src1,src2,src3,srcp,
r1t,r2t,r3t,
Rm,
s1t,s2t,s3t,
lp,nPML,Rc,pml_active,
C,
plot_interval,
wavefield_interval,
path,
path_pic,
path_model,
path_wavefield,
path_rec);

global data

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

IND_accuracy=9;
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
    pvd=paraview_collection(string(path,"/time_info"));
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
    CSV.write(string(path_model,"/source location.csv"),DataFrame([s1t' s2t' s3t'],:auto));
end

# create folder for wavefield
if path_wavefield!=nothing
    if isdir(path_wavefield)==0
        mkdir(path_wavefield)
    end
end

# create folder for rec
if path_rec!=nothing
    if isdir(path_rec)==0
        mkdir(path_rec)
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
# 6:lp+1
# nx-lp-5
beta1=@zeros(nx,ny,nz);
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

dtt1=@zeros(nx-IND_accuracy,ny,nz);
dtt2=@zeros(nx,ny-IND_accuracy,nz);
dtt3=@zeros(nx,ny,nz-IND_accuracy);

ax=copy(v1);
ax2=copy(v1);
ax3=copy(v1);
ax4=copy(v1);
ax5=copy(v1);
ax6=copy(v1);
ax7=copy(v1);
Ax=copy(v1);
Ax2=copy(v1);
Ax3=copy(v1);
Ax4=copy(v1);
Ax5=copy(v1);
Ax6=copy(v1);
Ax7=copy(v1);
ax_dt=copy(v1);
ax2_dt=copy(v1);
ax3_dt=copy(v1);
ax4_dt=copy(v1);
ax5_dt=copy(v1);
ax6_dt=copy(v1);
ax7_dt=copy(v1);

l=1;
# save wavefield
if path_wavefield!=nothing && wavefield_interval!=0
    if mod(l,wavefield_interval)==0
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/v1_",n_wavefield,".mat"),data);
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/v3_",n_wavefield,".mat"),data);
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
        data=zeros(nx,ny,nz);
        write2mat(string(path_wavefield,"/p_",n_wavefield,".mat"),data);
        n_wavefield=n_wavefield+1;
    end
end

pro_bar=Progress(nt,1,"forward_simulation...",50);
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

    @timeit ti "compute_sigma" @parallel compute_ax(dt,dx,dy,dz,C.inv_Qa,C.lambda,C.mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt);

    @timeit ti "compute_sigma" @parallel compute_sigma(dt,dx,dy,dz,C.inv_Qa,C.lambda,C.mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt);

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

    @timeit ti "compute_v" @parallel compute_v(dt,dx,dy,dz,C.rho,beta,
    v1,v2,v3,
    sigmas11_1_minus,
    sigmas22_2_minus,
    sigmas33_3_minus,
    sigmas23_2_plus,sigmas23_3_plus,
    sigmas13_1_plus,sigmas13_3_plus,
    sigmas12_1_plus,sigmas12_2_plus,
    p_1_minus,p_2_minus,p_3_minus);

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

if path_pic!=nothing && plot_interval!=0
    vtk_save(pvd);
end

return v1,v2,v3,R1,R2,R3,P
end
##
"
d/dx, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dx_inn(in,out)
    @inn(out)=@d_xi_12(in);
    return nothing
end
"
d/dy, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dy_inn(in,out)
    @inn(out)=@d_yi_12(in);
    return nothing
end
"
d/dz, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dz_inn(in,out)
    @inn(out)=@d_zi_12(in);
    return nothing
end
"
shift grid in x direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[6:end-4,iy,iz]=in[:,iy,iz]
"
@parallel_indices (iy,iz) function u_1_plus(in,out)
out[6:end-4,iy,iz]=in[:,iy,iz];
return nothing
end

"
shift grid in x direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[5:end-5,iy,iz]=in[:,iy,iz]
"
@parallel_indices (iy,iz) function u_1_minus(in,out)
out[5:end-5,iy,iz]=in[:,iy,iz];
return nothing
end

"
shift grid in y direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,6:end-4,iz]=in[ix,:,iz]
"
@parallel_indices (ix,iz) function u_2_plus(in,out)
out[ix,6:end-4,iz]=in[ix,:,iz];
return nothing
end

"
shift grid in y direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,5:end-5,iz]=in[ix,:,iz]
"
@parallel_indices (ix,iz) function u_2_minus(in,out)
out[ix,5:end-5,iz]=in[ix,:,iz];
return nothing
end

"
shift grid in z direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,iy,6:end-4]=in[ix,iy,:]
"
@parallel_indices (ix,iy) function u_3_plus(in,out)
out[ix,iy,6:end-4]=in[ix,iy,:];
return nothing
end

"
shift grid in z direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,iy,5:end-5]=in[ix,iy,:]
"
@parallel_indices (ix,iy) function u_3_minus(in,out)
out[ix,iy,5:end-5]=in[ix,iy,:];
return nothing
end

end
