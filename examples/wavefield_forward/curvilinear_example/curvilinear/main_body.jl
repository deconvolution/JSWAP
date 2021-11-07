## Initialization of JSWAP
include("./JSWAP_initialization.jl");
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
## Specify where the input file is
path_to_input=["./input_template.jl"];
## Run solvers
for I=1:length(path_to_input)
    include(path_to_input[I]);
    JSWAP.CPU_3D.JSWAP_CPU_3D_forward_isotropic_curvilinear_solver(input2);
end
##
using JSWAP
using JSWAP.CPU_3D.ParallelStencil,JSWAP.CPU_3D.ParallelStencil.FiniteDifferences3D
a=zeros(11,15,11);
a[4,8,4]=3;
b=zeros(11,12,11);
c=zeros(size(a));

@parallel JSWAP.CPU_3D.Dy_c(a,b)
@parallel (3:11-2,3:11-2) JSWAP.CPU_3D.uc_2_plus(b,c)
##
a=zeros(11,11,15);
a[4,4,8]=3;
b=zeros(11,11,12);
c=zeros(size(a));

@parallel JSWAP.CPU_3D.Dz_c(a,b)
@parallel (3:11-2,3:11-2) JSWAP.CPU_3D.uc_3_plus(b,c)
