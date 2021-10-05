## Initialization of JSWAP
include("./JSWAP_initialization.jl");
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
## Specify where the input file is
path_to_input=["./input_normal.jl"];
## Run solvers
for I=1:length(path_to_input)
    include(path_to_input[I]);
    JSWAP.CPU_3D.forward_solver(input2);
end
