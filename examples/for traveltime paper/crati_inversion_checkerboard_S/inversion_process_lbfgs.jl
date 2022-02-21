## import packages
using JSWAP,MATLAB,Statistics
## inversion paramemters
n_iteration=40;
max_gradient=200;
fu=2;

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
v=tt["vp"];
v[:] .=5346;

tt=readdir("./crati_traveltime_checkerboard_input/");
file_name=tt;
for I=1:size(tt,1)
    global R_true,s1,s2,s3,r1,r2,r3;
    tt2=JSWAP.readmat(string("./crati_traveltime_checkerboard_input/",tt[I]),"data");
    R_true=push!(R_true,tt2["Rp"][:,4]);
    s1=push!(s1,round.(Int64,tt2["S"][:,1]));
    s2=push!(s2,round.(Int64,tt2["S"][:,2]));
    s3=push!(s3,round.(Int64,tt2["S"][:,3]));
    r1=push!(r1,round.(Int64,tt2["Rp"][:,1]));
    r2=push!(r2,round.(Int64,tt2["Rp"][:,2]));
    r3=push!(r3,round.(Int64,tt2["Rp"][:,3]));
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
## initialization for l-bfgs
m=8;
# size of the problem
alpha=zeros(m,1);
s=zeros(nx,ny,nz,m);
y=copy(s);
tt_v=zeros(nx,ny,nz);
v_old=zeros(nx,ny,nz);
D=copy(v_old);
G=zeros(nx,ny,nz);
G_old=zeros(nx,ny,nz);
p=zeros(nx,ny,nz);
s_E=zeros(n_iteration,1);
DV=zeros(nx,ny,nz);

for l=1:n_iteration
    global v,max_graient,fu,G,s_E;
    if mod(l,10)==0
        fu=fu-1;
        if fu<=1
            fu=1;
        end
    end

    # batch for parallelization
    M=[0:30:size(s1,1);size(s1,1)];

    DV[:] .=0;
    if l==1
        E=0;
        for o=1:size(M,1)-1
            Threads.@threads for I=(M[o]+1):M[o+1]
                global R_true;
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

                E=E+JSWAP.norm(R_cal-R_true[I]',2);
                DV[:,:,:]=DV[:,:,:]+lambda ./v .^3;
            end
        end
        s_E[l]=E;
        s_fu[l]=fu;
        D[:,:,:]=DV;
        mat"""
        $G=imgaussfilt3($D,$fu);
        """
        G_old[:,:,:]=G[:,:,:];
        p[:,:,:]=G[:,:,:];
    else
        ## lbfgs
        E=0;
        for o=1:size(M,1)-1
            Threads.@threads for I=(M[o]+1):M[o+1]
                global R_true;
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

                E=E+JSWAP.norm(R_cal-R_true[I]',2);
                DV[:,:,:]=DV[:,:,:]+lambda ./v .^3;
            end
        end
        G_old[:,:,:]=G[:,:,:];
        D[:,:,:]=DV;
        mat"""
        $G=imgaussfilt3($D,$fu);
        """

        s_E[l]=E;
        s_fu[l]=fu;

        p[:,:,:]=G[:,:,:];

        for j=m:-1:max(m-l+2,1)
            alpha[j]=sum(s[:,:,:,j] .*p,dims=[1,2,3])[1]/sum(y[:,:,:,j] .*s[:,:,:,j],dims=[1,2,3])[1];
            p[:,:,:]=p-alpha[j]*y[:,:,:,j];
        end
        p[:,:,:]=sum(s[:,:,:,end] .*y[:,:,:,end],dims=[1,2,3])[1]/sum(y[:,:,:,end] .^2,dims=[1,2,3])[1]*p;
        for j=max(m-l+2,1):m
            beta=sum(y[:,:,:,j] .*p,dims=[1,2,3])[1]/sum(s[:,:,:,j] .*y[:,:,:,end],dims=[1,2,3])[1];
            p[:,:,:]=p+s[:,:,:,j]*(alpha[j]-beta);
        end

    end
    ## evaluate error again for line search
    tt_error=Inf;
    n=1;
    max_gradient2=max_gradient*2;
    while tt_error>s_E[l]+10^-4/maximum(abs.(p))*max_gradient2*sum(p.*G)
        max_gradient2=max_gradient2/2;
        tt_v[:,:,:]=v-p/maximum(abs.(p))*max_gradient2;
        E=0;
        for o=1:size(M,1)-1
            Threads.@threads for I=(M[o]+1):M[o+1]
                global R_true;
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
                v=tt_v,
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

                E=E+JSWAP.norm(R_cal-R_true[I]',2);
                tt_error=copy(E);
            end
        end

        n=n+1;
        if n>=4
            break;
        end
    end
    s_E[l]=E;

    s_max_gradient[l]=max_gradient2;
    v_old[:,:,:]=v[:,:,:];
    v[:,:,:]=tt_v[:,:,:];
    ## compute gradient again
    DV[:] .=0;
    for o=1:size(M,1)-1
        Threads.@threads for I=(M[o]+1):M[o+1]
            global R_true;
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

            DV[:,:,:]=DV[:,:,:]+lambda ./v .^3;
        end
    end

    D[:,:,:]=DV;

    mat"""
    $G=imgaussfilt3($D,$fu);
    """

    s[:,:,:,1:end-1]=s[:,:,:,2:end];
    y[:,:,:,1:end-1]=y[:,:,:,2:end];

    y[:,:,:,end]=G-G_old;
    s[:,:,:,end]=v-v_old;

    ## write
    JSWAP.CSV.write(string("./inversion_progress/E_",l,".csv"),
    JSWAP.DataFrame([reshape(s_E,length(s_E),) reshape(s_fu,length(s_fu),) reshape(s_max_gradient,length(s_max_gradient),)],:auto));

    vtkfile=JSWAP.vtk_grid(string("./inversion_progress/v_",l),X,
    Y,Z);
    vtkfile["v"]=v;
    vtkfile["G"]=G;
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
