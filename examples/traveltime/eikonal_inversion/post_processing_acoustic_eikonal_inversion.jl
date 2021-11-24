## read data
I=1;
include(path_to_input[I]);
v=input2.v;
nx=input2.nx;
ny=input2.ny;
nz=input2.nz;
X=input2.X;
Y=input2.Y;
Z=input2.Z;
h=input2.h;
max_gradient=input2.max_gradient;
fu=input2.fu;
n_iteration=input2.n_iteration;
s_E=zeros(n_iteration,1);
R_true=Vector{Vector{Float64}}();
s1=Vector{Vector{Int32}}();
s2=Vector{Vector{Int32}}();
s3=Vector{Vector{Int32}}();
r1=Vector{Vector{Int32}}();
r2=Vector{Vector{Int32}}();
r3=Vector{Vector{Int32}}();
s1t=Vector{Vector{Float64}}();
s2t=Vector{Vector{Float64}}();
s3t=Vector{Vector{Float64}}();
r1t=Vector{Vector{Float64}}();
r2t=Vector{Vector{Float64}}();
r3t=Vector{Vector{Float64}}();
path=Vector{String}();

for I=1:length(path_to_input)
    global R_true,s1,s2,s3,r1,r2,r3,s1t,s2t,s3t,r1t,r2t,r3t,path
    include(path_to_input[I]);
    R_true=push!(R_true,reshape(input2.R_true,length(input2.R_true),));
    s1=push!(s1,reshape(input2.s1,length(input2.s1),));
    s2=push!(s2,reshape(input2.s2,length(input2.s2),));
    s3=push!(s3,reshape(input2.s3,length(input2.s3),));
    r1=push!(r1,reshape(input2.r1,length(input2.r1),));
    r2=push!(r2,reshape(input2.r2,length(input2.r2),));
    r3=push!(r3,reshape(input2.r3,length(input2.r3),));
    s1t=push!(s1t,reshape(input2.s1t,length(input2.s1t),));
    s2t=push!(s2t,reshape(input2.s2t,length(input2.s2t),));
    s3t=push!(s3t,reshape(input2.s3t,length(input2.s3t),));
    r1t=push!(r1t,reshape(input2.r1t,length(input2.r1t),));
    r2t=push!(r2t,reshape(input2.r2t,length(input2.r2t),));
    r3t=push!(r3t,reshape(input2.r3t,length(input2.r3t),));
    path=push!(path,input2.path);
end
##
for l=1:n_iteration
    global v
    dv=zeros(input2.nx,input2.ny,input2.nz,length(path_to_input));
    E=zeros(length(path_to_input),1);
    Threads.@threads for I=1:length(path_to_input)
        T,R_cal=JSWAP.eikonal.acoustic_eikonal_forward(nx,ny,nz,
        h,v,s1[I],s2[I],s3[I],s1t[I],s2t[I],s3t[I],
        r1[I],r2[I],r3[I],
        r1t[I],r2t[I],r3t[I],X,Y,Z,path[I],0);
        lambda=JSWAP.eikonal.acoustic_eikonal_adjoint(nx,ny,nz,
        h,T,r1[I],r2[I],r3[I],s1[I],s2[I],s3[I],R_cal,R_true[I]);
        E[I]=JSWAP.norm(R_cal-R_true[I],2);
        dv[:,:,:,I]=lambda ./v .^3;
    end
    s_E[l]=sum(E);
    DV=reshape(sum(dv,dims=4),input2.nx,input2.ny,input2.nz);
    mat"""
    $DV=imgaussfilt($DV,$fu);
    """
    v=v-max_gradient/maximum(abs.(DV))*DV;
    vtkfile=JSWAP.vtk_grid(string("./inversion_progress/v_",l),input2.X,
    input2.Y,input2.Z);
    vtkfile["v"]=input2.v;
    vtkfile["DV"]=DV;
    JSWAP.vtk_save(vtkfile);
    print("\n",l)
end
