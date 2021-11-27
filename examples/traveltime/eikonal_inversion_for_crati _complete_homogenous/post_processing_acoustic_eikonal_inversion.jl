## read data
include(path_to_input[1]);
v=input2.v;
nx=input2.nx;
ny=input2.ny;
nz=input2.nz;
X=input2.X;
Y=input2.Y;
Z=input2.Z;
h=input2.dx;
s1=input2.s1;
s2=input2.s2;
s3=input2.s3;
r1=input2.r1;
r2=input2.r2;
r3=input2.r3;
s1t=input2.s1t;
s2t=input2.s2t;
s3t=input2.s3t;
r1t=input2.r1t;
r2t=input2.r2t;
r3t=input2.r3t;
n_iteration=input2.n_iteration;
path=input2.path;
R_true=input2.R_true;

s_E=zeros(n_iteration,1);
s_fu=copy(s_E);
alp=0;
## create struct data
mutable struct data2
    X
    Y
    Z
    nx
    ny
    nz
    dx
    dy
    dz
    vp
end
data=data2(0,0,0,0,0,0,0,0,0,0);
## inversion
n_decrease_fu=0;
max_gradient=input2.max_gradient;
fu=input2.fu;
for l=1:n_iteration
    global v,n_decrease_fu,alp,max_gradient,fu
    include(path_to_input[1]);
    # adjust max_gradient and fu
    if mod(l,100)==0 && l!=1
        max_gradient=max_gradient/2;
        fu=fu/2;
    end
    DV=zeros(input2.nx,input2.ny,input2.nz);
    E=zeros(size(input2.s1,1),1);

    M=[0:30:size(input2.s1,1);size(input2.s1,1)];
    for m=1:size(M,1)-1
        Threads.@threads for I=(M[m]+1):M[m+1]
            T,R_cal=JSWAP.eikonal.acoustic_eikonal_forward(nx=nx,
            ny=ny,
            nz=nz,
            h=h,
            v=v,
            s1=s1[I][1],
            s2=s2[I][1],
            s3=s3[I][1],
            T0=0,
            s1t=s1t[I][1],
            s2t=s2t[I][1],
            s3t=s3t[I][1],
            r1=r1[I],
            r2=r2[I],
            r3=r3[I],
            r1t=r1t[I],
            r2t=r2t[I],
            r3t=r3t[I],
            X=X,
            Y=Y,
            Z=Z,
            path=path,
            write_t=0);

            lambda=JSWAP.eikonal.acoustic_eikonal_adjoint(nx=nx,
            ny=ny,
            nz=nz,
            h=h,
            T=T,
            r1=r1[I],
            r2=r2[I],
            r3=r3[I],
            s1=s1[I][1],
            s2=s2[I][1],
            s3=s3[I][1],
            R_cal=R_cal,
            R_true=R_true[I]);
            E[I]=JSWAP.norm(R_cal-R_true[I],2);
            DV[:,:,:]=DV[:,:,:]+lambda ./v .^3;
        end
    end

    s_E[l]=sum(E);
    #=
    if l>=2
        if s_E[l]>s_E[l-1]
            n_decrease_fu=n_decrease_fu+1;
        end
    end
    if n_decrease_fu>=4
        fu=round(Int64,fu-4);
        if fu<=1
            fu=1;
        end
    end
    =#
    s_fu[l]=fu;
    JSWAP.CSV.write(string("./inversion_progress/E_",l,".csv"),
    JSWAP.DataFrame([reshape(s_E,length(s_E),) reshape(s_fu,length(s_fu),)],:auto));
    mat"""
    $DV=imgaussfilt3($DV,$fu);
    """
    if l==1
        alp=max_gradient/maximum(abs.(DV));
    end
    v=v-alp*DV;
    vtkfile=JSWAP.vtk_grid(string("./inversion_progress/v_",l),input2.X,
    input2.Y,input2.Z);
    vtkfile["v"]=v;
    vtkfile["DV"]=DV;
    JSWAP.vtk_save(vtkfile);
    ## write velocity
    data.X=input2.X;
    data.Y=input2.Y;
    data.Z=input2.Z;
    data.vp=v;
    data.nx=input2.nx;
    data.ny=input2.ny;
    data.nz=input2.nz;
    data.dx=input2.dx;
    data.dy=input2.dy;
    data.dz=input2.dz;
    file=JSWAP.matopen(string("./inversion_progress/", "output_vp_",l,".mat"), "w");
    write(file,"data",data);
    close(file);

    println("iteration = ",l);
end

##
# compute new source location
#=
for J=1:size(input2.s1,1)
global Tt;
Tt=zeros(nx,ny,nz);
M=[0:30:size(input2.s1,1);size(input2.s1,1)];
for m=1:size(M,1)
Threads.@threads for I=M[m]+1:M[m+1]
global Tt;
T,R_cal=JSWAP.eikonal.acoustic_eikonal_forward(nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
s1=r1[J][I],
s2=r2[J][I],
s3=r3[J][I],
T0=-R_true[J][I],
r1t=r1t[J][I],
r2t=r2t[J][I],
r3t=r3t[J][I],
r1=r1[I],
r2=r2[I],
r3=r3[I],
r1t=r1t[I],
r2t=r2t[I],
r3t=r3t[I],
X=X,
Y=Y,
Z=Z,
path=path,
write_t=0);
Tt=Tt+abs.(T);
end
end
tt=findmax(Tt);
tt2=findmax(abs.(Z[1,1,:]));

s1[J][1]=tt[2][1];
s2[J][1]=tt[2][2];
s3[J][1]=tt[2][3];

s1t[J][1]=s1[J][1]*h;
s2t[J][1]=s2[J][1]*h;
s3t[J][1]=(s3[J][1]-tt2[2])*h;

JSWAP.CSV.write(string("./inversion_progress/source_",file_name[J],".csv"),
JSWAP.DataFrame([s1t[J],s2t[J],s3t[J]],:auto));
=#
