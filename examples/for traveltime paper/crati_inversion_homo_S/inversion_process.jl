## import packages
using JSWAP,MATLAB
## inversion paramemters
n_iteration=50;
max_gradient=100;
fu=8;

R_true=Vector{Vector{Float64}}();
s1=Vector{Vector{Int64}}();
s2=Vector{Vector{Int64}}();
s3=Vector{Vector{Int64}}();
r1=Vector{Vector{Int64}}();
r2=Vector{Vector{Int64}}();
r3=Vector{Vector{Int64}}();
s1t=Vector{Vector{Float64}}();
s2t=Vector{Vector{Float64}}();
s3t=Vector{Vector{Float64}}();
r1t=Vector{Vector{Float64}}();
r2t=Vector{Vector{Float64}}();
r3t=Vector{Vector{Float64}}();
##
tt=readmat("./m.mat","data");
nx=round(Int64,tt["nx"]);
ny=round(Int64,tt["ny"]);
nz=round(Int64,tt["nz"]);
X=tt["X"];
Y=tt["Y"];
Z=tt["Z"];
h=tt["dx"];
dx=h;
dy=h;
dz=h;
v=tt["vs"];

tt=readdir("./crati_traveltime_input/");
file_name=tt;
for I=1:size(tt,1)
    global R_true,s1,s2,s3,r1,r2,r3;
    tt2=JSWAP.readmat(string("./crati_traveltime_input/",tt[I]),"data");
    if size(tt2["Rs"],1)!=0
        R_true=push!(R_true,tt2["Rs"][:,4]);
        s1=push!(s1,round.(Int64,tt2["S"][:,1]));
        s2=push!(s2,round.(Int64,tt2["S"][:,2]));
        s3=push!(s3,round.(Int64,tt2["S"][:,3]));
        r1=push!(r1,round.(Int64,tt2["Rs"][:,1]));
        r2=push!(r2,round.(Int64,tt2["Rs"][:,2]));
        r3=push!(r3,round.(Int64,tt2["Rs"][:,3]));
    end
end
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
zero_Z=findmin(abs.(Z[1,1,:]));
# receiver true location x
r1t=r1*dx;
# receiver true location y
r2t=r2*dx;
# receiver true location z
s3t=copy(s3);
for i=1:size(s3t,1)
    s3t[i]=h*(s3[i] .-zero_Z[2]);
end
# source true location x
s1t=s1*dx;
# source true location y
s2t=s2*dx;
# source true location z
r3t=copy(r3);
for i=1:size(r3t,1)
    r3t[i]=h*(r3[i] .-zero_Z[2]);
end

tt=zeros(1,size(s1t,1));
tt2=zeros(1,size(s2t,1));
tt3=zeros(1,size(s3t,1));
for i=1:size(s1t,1)
    tt[i]=s1t[i][1];
    tt2[i]=s2t[i][1];
    tt3[i]=s3t[i][1];
end
JSWAP.CSV.write(string("./source_overview2.csv"),
JSWAP.DataFrame([tt' tt2' tt3'],:auto));
## plot
write_t=0;
##
s_E=zeros(n_iteration,1);
s_fu=zeros(n_iteration,1);
s_max_gradient=zeros(n_iteration,1);
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
td=0;
for l=1:n_iteration
    global v,n_decrease_fu,alp,max_gradient,fu,td;


    DV=zeros(nx,ny,nz);
    E=zeros(size(s1,1),1);

    M=[0:50:size(s1,1);size(s1,1)];

    for m=1:size(M,1)-1
        Threads.@threads for I=(M[m]+1):M[m+1]
            global R_true,td;
            input_s1=zeros(Int64,1,1);
            input_s2=zeros(Int64,1,1);
            input_s3=zeros(Int64,1,1);
            input_s1[:] .=s1[I][1];
            input_s2[:] .=s2[I][1];
            input_s3[:] .=s3[I][1];

            T,R_cal=JSWAP.eikonal.acoustic_eikonal_forward(nx=nx,
            ny=ny,
            nz=nz,
            h=h,
            v=v,
            s1=input_s1,
            s2=input_s2,
            s3=input_s3,
            T0=0,
            s1t=s1t[I][1],
            s2t=s2t[I][1],
            s3t=s3t[I][1],
            r1=r1[I]',
            r2=r2[I]',
            r3=r3[I]',
            r1t=r1t[I]',
            r2t=r2t[I]',
            r3t=r3t[I]',
            X=X,
            Y=Y,
            Z=Z,
            path=string("./inversion_configuration/",I,"/"),
            write_t=0);

            lambda=JSWAP.eikonal.acoustic_eikonal_adjoint(nx=nx,
            ny=ny,
            nz=nz,
            h=h,
            T=T,
            r1=reshape(r1[I],1,size(r1[I],1)),
            r2=reshape(r2[I],1,size(r2[I],1)),
            r3=reshape(r3[I],1,size(r3[I],1)),
            s1=s1[I][1],
            s2=s2[I][1],
            s3=s3[I][1],
            R_cal=R_cal,
            R_true=R_true[I]',
            N=ones(size(R_cal))*(-1));

            E[I]=JSWAP.norm(R_cal-R_true[I]',2);
            DV[:,:,:]=DV[:,:,:]+lambda ./v .^3;
        end
    end
    s_E[l]=sum(E);

    s_fu[l]=fu;

    mat"""
    $DV=imgaussfilt3($DV,$fu);
    """
    if l>=2
        if s_E[l]>s_E[l-1]
            td=td+1;
        end
    end
    #=
    if td==5
        fu=fu-1;
        max_gradient=max_gradient*.9;
        if fu<=1
            fu=1;
        end
        if max_gradient<=100 && l<=100
            max_gradient=100;
        end
        td=0;
    end
    =#
    if mod(l,10)==0
        fu=fu-1;
        max_gradient=max_gradient*.8;
        if fu<=1
            fu=1;
        end
        if max_gradient<=100 && l<=100
            max_gradient=100;
        end
    end

    s_max_gradient[l]=max_gradient;
    JSWAP.CSV.write(string("./inversion_progress/E_",l,".csv"),
    JSWAP.DataFrame([reshape(s_E,length(s_E),) reshape(s_fu,length(s_fu),) reshape(s_max_gradient,length(s_max_gradient),)],:auto));

    v=v-max_gradient/maximum(abs.(DV))*DV;
    vtkfile=JSWAP.vtk_grid(string("./inversion_progress/v_",l),X,
    Y,Z);
    vtkfile["v"]=v;
    vtkfile["DV"]=DV;
    JSWAP.vtk_save(vtkfile);
    ## write velocity
    data.X=X;
    data.Y=Y;
    data.Z=Z;
    data.vp=v;
    data.nx=nx;
    data.ny=ny;
    data.nz=nz;
    data.dx=dx;
    data.dy=dy;
    data.dz=dz;
    file=JSWAP.matopen(string("./inversion_progress/", "output_vp_",l,".mat"), "w");
    write(file,"data",data);
    close(file);
    println("iteration = ",l);
end
