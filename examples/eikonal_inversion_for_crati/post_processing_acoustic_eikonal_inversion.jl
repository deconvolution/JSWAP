## read data
include(path_to_input[1]);
v=input2.v;
nx=input2.nx;
ny=input2.ny;
nz=input2.nz;
X=input2.X;
Y=input2.Y;
Z=input2.Z;
h=input2.h;
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
##
n_decrease_fu=0;
for l=1:n_iteration
    include(path_to_input[1]);
    max_gradient=input2.max_gradient;
    fu=input2.fu;
    global v,n_decrease_fu
    dv=zeros(input2.nx,input2.ny,input2.nz,size(input2.s1,1));
    E=zeros(length(path_to_input),1);
    Threads.@threads for I=1:length(path_to_input)
        T,R_cal=JSWAP.eikonal.acoustic_eikonal_forward(nx,ny,nz,
        h,v,s1[I],s2[I],s3[I],s1t[I],s2t[I],s3t[I],
        r1[I],r2[I],r3[I],
        r1t[I],r2t[I],r3t[I],X,Y,Z,path,0);
        lambda=JSWAP.eikonal.acoustic_eikonal_adjoint(nx,ny,nz,
        h,T,r1[I],r2[I],r3[I],s1[I],s2[I],s3[I],R_cal,R_true[I]);
        E[I]=JSWAP.norm(R_cal-R_true[I],2);
        dv[:,:,:,I]=lambda ./v .^3;
    end
    s_E[l]=sum(E);

    DV=reshape(sum(dv,dims=4),round(Int64,nx),round(Int64,ny),round(Int64,nz));
    if l>=2
        if s_E[l]>s_E[l-1]
            n_decrease_fu=n_decrease_fu+1;
        end
    end
    if n_decrease_fu>=4
        fu=round(Int64,fu/2);
    end
    s_fu[l]=fu;
    JSWAP.CSV.write(string("./inversion_progress/E_",l,".csv"),
    JSWAP.DataFrame([reshape(s_E,length(s_E),) reshape(s_fu,length(s_fu),)],:auto));
    mat"""
    $DV=imgaussfilt3($DV,$fu);
    """

    v=v-max_gradient/maximum(abs.(DV))*DV;
    vtkfile=JSWAP.vtk_grid(string("./inversion_progress/v_",l),input2.X,
    input2.Y,input2.Z);
    vtkfile["v"]=input2.v;
    vtkfile["DV"]=DV;
    JSWAP.vtk_save(vtkfile);
    print("\n",l)
end
