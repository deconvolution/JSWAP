## Initialization of JSWAP
include("JSWAP_initialization.jl");
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
## Specify where the input file is
path_to_input=["./input_test.jl"];
## Run solvers
I=1;
include(path_to_input[I]);
v1,v2,v3,R1,R2,R3,P=JSWAP.CPU_3D.isotropic_forward_solver(input2);
